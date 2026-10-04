use crate::enzyme::{group_digests, Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{validate_mods, validate_var_mods, ModificationSpecificity, VarModEntry};
use crate::peptide::Peptide;
use dashmap::DashSet;
use fnv::{FnvBuildHasher, FnvHashSet};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Arc;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    pub semi_enzymatic: Option<bool>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme: Enzyme::new(
                &en.cleave_at.clone().unwrap_or_else(|| "KR".into()),
                // An omitted `restrict` means no restriction for a custom `cleave_at`, but
                // when `cleave_at` is omitted too the enzyme is trypsin, whose rule is "not
                // before P" (DOCS.md). `{"missed_cleavages": 2}` used to drop that rule.
                &en.restrict.unwrap_or_else(|| match en.cleave_at {
                    None => "P".into(),
                    Some(_) => "".into(),
                }),
                en.c_terminal.unwrap_or(true),
                en.semi_enzymatic.unwrap_or(false),
            ),
        }
    }
}

#[derive(Deserialize, Default)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,

    pub enzyme: Option<EnzymeBuilder>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Which kind of fragment ions to generate (a, b, c, x, y, z)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids
    pub static_mods: Option<HashMap<String, f32>>,
    /// Variable modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`15.9949`) or an object with `mass` and
    /// optional `max_count` fields (`{"mass": 15.9949, "max_count": 1}`).
    pub variable_mods: Option<HashMap<String, Vec<VarModEntry>>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Hard cap on the total peptide variants generated per input peptide,
    /// including its unmodified form. Values below 1 are normalized to 1.
    /// Variants with fewer PTMs are preferred (generated first).
    pub max_combinations: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Pre-filter the database with a minimal amount of memory at the cost of speed
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods),
            variable_mods: validate_var_mods(self.variable_mods),
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.expect("A fasta file must be provided!"),
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_low_memory: self.prefilter_low_memory.unwrap_or(true),
        }
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    pub static_mods: HashMap<ModificationSpecificity, f32>,
    pub variable_mods: HashMap<ModificationSpecificity, Vec<VarModEntry>>,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    pub prefilter_low_memory: bool,
}

impl Parameters {
    /// Flatten variable modifications into a stable order. This matters when
    /// `max_combinations` truncates variants: equivalent configurations must
    /// retain the same variants regardless of randomized `HashMap` iteration.
    fn variable_modifications(&self) -> Vec<(ModificationSpecificity, f32, Option<usize>)> {
        let mut mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().enumerate().map(|(entry_order, entry)| {
                    (*specificity, entry_order, entry.mass(), entry.max_count())
                })
            })
            .collect::<Vec<_>>();
        mods.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        mods.into_iter()
            .map(|(specificity, _, mass, max_count)| (specificity, mass, max_count))
            .collect()
    }

    pub fn auto_calculate_prefilter_chunk_size(&mut self, fasta: &Fasta) {
        const MAX_PEPS_PER_CHUNK: usize = 2usize.pow(23);
        self.prefilter_chunk_size = match self.prefilter_chunk_size {
            0 => {
                let enzyme = self.enzyme.clone().into();
                let total_unmodified_pep_count: usize = fasta.digest(&enzyme).len();
                let combination_factor = if self.max_variable_mods >= usize::BITS as usize {
                    usize::MAX
                } else {
                    1usize << self.max_variable_mods
                };
                let mut mod_count_estimate =
                    (self.variable_mods.len() + 1).saturating_mul(combination_factor);
                if let Some(max_combinations) = self.max_combinations {
                    mod_count_estimate = mod_count_estimate.min(max_combinations);
                }
                let chunk_count = mod_count_estimate.saturating_mul(total_unmodified_pep_count)
                    / MAX_PEPS_PER_CHUNK;
                fasta
                    .targets
                    .len()
                    .checked_div(chunk_count)
                    .unwrap_or(fasta.targets.len())
            }
            x => x,
        };
    }

    pub fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme.clone().into();
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest(&enzyme);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_digests(digests);
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );

        let mods = self.variable_modifications();

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .for_each(|digest| {
                targets.insert(digest.reference.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .map(Peptide::try_from)
            .filter_map(Result::ok)
            .flat_map_iter(|peptide| {
                peptide
                    .apply(
                        &mods,
                        &self.static_mods,
                        self.max_variable_mods,
                        self.max_combinations,
                    )
                    .into_iter()
                    .filter(|peptide| {
                        peptide.monoisotopic >= self.peptide_min_mass
                            && peptide.monoisotopic <= self.peptide_max_mass
                    })
                    .flat_map(|peptide| {
                        if self.generate_decoys {
                            vec![peptide.reverse(), peptide].into_iter()
                        } else {
                            vec![peptide].into_iter()
                        }
                    })
                    .filter(|peptide| !peptide.decoy || !targets.contains(&(peptide.sequence[..])))
            })
            .collect::<Vec<_>>();

        // `digest` already removed decoys colliding with targets
        Self::sort_and_dedup(&mut target_decoys);

        target_decoys
    }

    /// Merge peptides from separately digested FASTA chunks (prefilter): drop decoys that
    /// collide with a target, then sort and deduplicate.
    pub fn reorder_peptides(target_decoys: &mut Vec<Peptide>) {
        // A decoy with the sequence of a target must go, as in `digest`. With a
        // prefiltered database the chunks are digested separately, so a decoy from one
        // chunk can equal a target from another; merging them below used to append the
        // decoy's (unrelated) proteins to the target and destroy its uniqueness.
        let targets: FnvHashSet<Arc<[u8]>> = target_decoys
            .iter()
            .filter(|p| !p.decoy)
            .map(|p| p.sequence.clone())
            .collect();
        target_decoys.retain(|p| !p.decoy || !targets.contains(&p.sequence));
        drop(targets);
        Self::sort_and_dedup(target_decoys);
    }

    fn sort_and_dedup(target_decoys: &mut Vec<Peptide>) {
        log::trace!("sorting and deduplicating peptides");
        let init_size = target_decoys.len();
        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.monoisotopic == keep.monoisotopic
                && remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

        let num_dropped = init_size - target_decoys.len();
        log::trace!(
            "dropped {} t/d pairs, remaining {}",
            num_dropped,
            target_decoys.len(),
        );
    }

    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        let target_decoys = self.digest(&fasta);
        self.build_from_peptides(target_decoys)
    }

    /// Ions of `peptide` that go into the fragment index (b1, b2, y1, y2... are
    /// excluded according to `min_ion_index`)
    fn index_ions<'a>(
        &'a self,
        peptide: &'a Peptide,
    ) -> impl Iterator<Item = crate::ion_series::Ion> + 'a {
        self.ion_kinds
            .iter()
            .flat_map(move |kind| IonSeries::new(peptide, *kind).enumerate())
            .filter(move |(ion_idx, ion)| match ion.kind {
                // Don't store b1, b2, y1, y2 ions for preliminary scoring
                Kind::A | Kind::B | Kind::C => (ion_idx + 1) > self.min_ion_index,
                Kind::X | Kind::Y | Kind::Z => {
                    peptide.sequence.len().saturating_sub(1) - ion_idx > self.min_ion_index
                }
            })
            .map(|(_, ion)| ion)
    }

    /// All theoretical fragments of `peptides`, in [`binned_fragments`] order: by m/z
    /// ([`f32::total_cmp`]), fragments of equal m/z by peptide index.
    fn sorted_fragments(&self, peptides: &[Peptide]) -> Vec<Theoretical> {
        self.sorted_fragments_with(
            peptides,
            4096,
            MAX_BINS,
            (rayon::current_num_threads() * 4).clamp(1, 256),
        )
    }

    fn sorted_fragments_with(
        &self,
        peptides: &[Peptide],
        chunk_size: usize,
        max_bins: usize,
        max_groups: usize,
    ) -> Vec<Theoretical> {
        // The fragments of a peptide are lighter than the peptide plus the terminal group
        // of an a/c/x/z ion, so practically every fragment falls into a regular bin of
        // [GUESS_LOW, heaviest peptide + GUESS_MARGIN] without a separate pass over all
        // ions for their m/z range. Anything outside (e.g. after a negative modification
        // mass) lands in an overflow bin and still ends up in order.
        let heaviest = peptides
            .par_iter()
            .map(|p| p.monoisotopic)
            .reduce(|| f32::NEG_INFINITY, f32::max);
        binned_fragments(
            peptides.len(),
            |ix| peptides[ix].sequence.len(),
            |ix| {
                self.index_ions(&peptides[ix])
                    .map(|ion| ion.monoisotopic_mass)
            },
            MzBins::new(GUESS_LOW, heaviest + GUESS_MARGIN, max_bins),
            chunk_size,
            max_groups,
        )
    }

    pub fn build_from_peptides(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        log::trace!("generating fragments");

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by
        // [`SpectrumProcessor`] or during scoring - all theoretical
        // fragments are monoisotopic/uncharged
        // All of our theoretical fragments, sorted by m/z from low to high
        let mut fragments = self.sorted_fragments(&target_decoys);
        log::trace!("finalizing index");

        // Now, we bucket all of our theoretical fragments, and within each bucket
        // sort by precursor m/z - and save the minimum *fragment* m/z in a separate
        // vector so that we can perform an efficient binary search to reduce
        // the number of in silico fragments we evaluate
        //
        // Imagine our theoretical fragments look like this
        //
        // Fragment        A      B       C       D       E       F       G       H
        // Fragment m/z [ 1.0    1.2     1.3     2.5     2.5     2.6     3.5     4.0 ]
        // Parent m/z   [ 500    439     291     800     142     515     517     232 ]
        //
        // If we apply a bucket size of 4 we will end up with the following:
        //
        // Fragment        C      B       A       D       E       H       F       G
        // Fragment m/z [ 1.3    1.2     1.0     2.5     2.5     4.0     3.5     2.6 ]
        // Parent m/z   [ 291    439     500     800     142     232     515     517 ]
        //              |___________________________|   |____________________________|
        //               Bucket 1: min m/z 1.0          Bucket 2: min m/z 2.5
        //
        // * Example query: Fragment m/z 1.3 - 1.9 & Precursor m/z: 450 - 900
        // 1) Perform a binary search to narrow down our window to Bucket 1 only
        //      * Bucket 2 has a min m/z outside of our query range - nothing here can match
        //
        // Fragment        C      B       A       D
        // Fragment m/z [ 1.3    1.2     1.0     2.5
        // Parent m/z   [ 291    439     500     800
        //                            |_____________|
        //                                    ^
        //                                    |
        // Window with matching precursors ___|

        // and within Bucket 1, we can perform another binary search to find fragments
        // matching our desired precursor m/z tolerance

        // Each bucket (page) is sorted by peptide index; fragments of the same peptide
        // keep their m/z order, so the index depends on the fragments alone (see
        // [`binned_fragments`]).
        let peptide_bits = index_bits(target_decoys.len());
        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map_init(SortScratch::default, |scratch, chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                sort_by_peptide(chunk, peptide_bits, scratch);
                min
            })
            .collect::<Vec<_>>();

        let potential_mods = self
            .variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries.iter().map(|entry| (*specificity, entry.mass()))
            })
            .collect::<Vec<(ModificationSpecificity, f32)>>();

        let page_skip = page_skip(&fragments, self.bucket_size);

        IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            page_skip,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            decoy_tag: self.decoy_tag,
        }
    }
}

/// Ask the kernel to back `buf` with transparent huge pages before it is touched.
///
/// The fragment index is a multi-GB array searched with random, cache-missing binary
/// searches; with 4 KiB pages nearly every probe also misses the TLB. Most
/// distributions ship THP in `madvise` mode, so without this hint the index stays on
/// small pages. PXD041421: index build -20%, search -2..-9% (same effect as
/// `GLIBC_TUNABLES=glibc.malloc.hugetlb=1`). A no-op elsewhere or if THP is disabled.
fn advise_huge_pages<T>(buf: &mut [T]) {
    #[cfg(target_os = "linux")]
    {
        // SAFETY: sysconf has no preconditions
        let page = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
        let page = if page > 0 { page as usize } else { 4096 };
        let start = buf.as_mut_ptr() as usize;
        let end = start + std::mem::size_of_val(buf);
        let aligned = (start + page - 1) & !(page - 1);
        if end > aligned {
            // SAFETY: the range lies inside `buf`'s allocation; MADV_HUGEPAGE only
            // changes how the kernel backs these pages, never their contents.
            unsafe {
                libc::madvise(
                    aligned as *mut libc::c_void,
                    end - aligned,
                    libc::MADV_HUGEPAGE,
                );
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = buf;
}

/// Upper bound on the number of regular m/z bins of [`binned_fragments`]
const MAX_BINS: usize = 8192;
/// Fragment masses expected in `GUESS_LOW..=heaviest peptide + GUESS_MARGIN` (Da);
/// others go to the overflow bins
const GUESS_LOW: f32 = 50.0;
const GUESS_MARGIN: f32 = 200.0;
/// Bins up to this size are sorted by comparison instead of by counting
const SMALL_BIN: usize = 256;

/// Bins of the fragment m/z for [`binned_fragments`], in [`mz_key`] order: bin 0 holds
/// every key below that of `lowest`, the last bin every key above that of `highest`
/// (negative masses, infinities and NaN included), and the regular bins in between
/// split the keys of `lowest..=highest` by their top bits.
#[derive(Clone, Copy, Debug, PartialEq)]
struct MzBins {
    lo: u32,
    hi: u32,
    /// Number of low key bits that a regular bin does not fix
    shift: u32,
    /// Number of bins, the two overflow bins included
    len: usize,
}

impl MzBins {
    /// At most `max_bins` regular bins (at least one) over the keys of `lowest..=highest`
    fn new(lowest: f32, highest: f32, max_bins: usize) -> Self {
        let lo = mz_key(lowest);
        let hi = mz_key(highest).max(lo);
        let mut shift = 0;
        // terminates: at shift 31 the difference is at most 1
        while ((hi >> shift) - (lo >> shift)) as usize >= max_bins.max(2) {
            shift += 1;
        }
        MzBins {
            lo,
            hi,
            shift,
            len: ((hi >> shift) - (lo >> shift)) as usize + 3,
        }
    }

    #[inline(always)]
    fn of(&self, mz: f32) -> usize {
        let k = mz_key(mz);
        if k < self.lo {
            0
        } else if k > self.hi {
            self.len - 1
        } else {
            ((k >> self.shift) - (self.lo >> self.shift)) as usize + 1
        }
    }

    fn is_overflow(&self, bin: usize) -> bool {
        bin == 0 || bin + 1 == self.len
    }
}

/// Reusable buffers of the bin and page sorts
#[derive(Default)]
struct SortScratch {
    items: Vec<Theoretical>,
    counts: Vec<u32>,
}

/// The fragments of items (peptides) `0..n` in (m/z, item index) order: by m/z
/// ([`f32::total_cmp`]), and fragments of equal m/z by item index. `ions(ix)` yields
/// the fragment m/z of item `ix` and must yield the same values every time it is
/// called; `weight(ix)` estimates their number.
///
/// This replaces a global sort of all fragments (`par_sort_unstable`, ~3*10^8 of them
/// for a large search space, which moved the 2.4 GB array ~log2(n) times through
/// memory and stopped scaling at ~16 threads) by a two-pass build:
///
/// 1. The items are split into at most `max_groups` contiguous groups of whole
///    `chunk_size` chunks with about equal weight. Each group counts its fragments per
///    m/z bin ([`MzBins`]).
/// 2. Each group generates its fragments a second time and writes every one straight
///    into its own range of its bin in the final allocation. Within a bin, the ranges
///    follow group order and every group writes in item order, so a bin holds its
///    fragments in item order.
/// 3. Each bin is put in order by the key bits that it does not fix, with a stable
///    counting sort (fragments of equal m/z keep item order), or, for small and
///    overflow bins, by comparison of (m/z key, item).
///
/// The result is a function of the multiset of (item, m/z) pairs alone: two fragments
/// that agree in both are the same bytes, so neither the thread count, the grouping
/// nor the order in which `ions` yields the m/z of one item can show. A global
/// unstable sort orders fragments of equal m/z arbitrarily, so the index can differ
/// from that sort's in the order of such fragments and therefore in which of two pages
/// a fragment at a page boundary lands; the search does not depend on either (every
/// lookup visits all pages whose m/z range meets its window and tests every fragment
/// in them against the same bounds, and the minimum m/z of every page is the same).
///
/// Memory: the final array (no second copy of the fragments), one count or write
/// position per group and bin, and a scratch buffer per worker of the size of the
/// largest bin it sorts.
fn binned_fragments<W, F, I>(
    n: usize,
    weight: W,
    ions: F,
    bins: MzBins,
    chunk_size: usize,
    max_groups: usize,
) -> Vec<Theoretical>
where
    W: Fn(usize) -> usize + Sync,
    F: Fn(usize) -> I + Sync,
    I: Iterator<Item = f32>,
{
    assert!(
        n <= u32::MAX as usize + 1,
        "too many peptides for a 32-bit index"
    );
    if n == 0 {
        return Vec::new();
    }
    let nbins = bins.len;
    let groups = balanced_groups(n, chunk_size, max_groups, weight);

    // pass 1: the number of fragments per group and bin (one row per group)
    let mut table = vec![0usize; groups.len() * nbins];
    table
        .par_chunks_mut(nbins)
        .zip(groups.par_iter())
        .for_each(|(counts, range)| {
            for ix in range.clone() {
                for mz in ions(ix) {
                    counts[bins.of(mz)] += 1;
                }
            }
        });
    // Counts -> first write position of every (group, bin): bins in m/z order, and
    // within a bin the groups in item order. `bounds[b]..bounds[b + 1]` is bin b.
    let mut next = vec![0usize; nbins];
    for counts in table.chunks(nbins) {
        next.iter_mut().zip(counts).for_each(|(n, c)| *n += c);
    }
    let mut bounds = Vec::with_capacity(nbins + 1);
    let mut total = 0usize;
    for n in next.iter_mut() {
        bounds.push(total);
        total += std::mem::replace(n, total);
    }
    bounds.push(total);
    for counts in table.chunks_mut(nbins) {
        for (c, n) in counts.iter_mut().zip(next.iter_mut()) {
            let count = std::mem::replace(c, *n);
            *n += count;
        }
    }
    drop(next);
    if total == 0 {
        return Vec::new();
    }

    // pass 2: every fragment straight into its group's range of its bin
    let mut fragments: Vec<Theoretical> = Vec::with_capacity(total);
    advise_huge_pages(fragments.spare_capacity_mut());
    struct Out(*mut Theoretical);
    // SAFETY: shared only for writes to disjoint, bounds-checked elements (see below)
    unsafe impl Sync for Out {}
    let out = Out(fragments.as_mut_ptr());
    let (out, table, bounds) = (&out, &table, &bounds);
    let filled = groups.par_iter().enumerate().all(|(g, range)| {
        // (next write position, end) of group g's range in every bin: the range ends
        // where that of the next group in the same bin starts
        let start = &table[g * nbins..(g + 1) * nbins];
        let mut slots: Vec<(usize, usize)> = match table.get((g + 1) * nbins..(g + 2) * nbins) {
            Some(end) => start.iter().copied().zip(end.iter().copied()).collect(),
            None => start
                .iter()
                .copied()
                .zip(bounds[1..].iter().copied())
                .collect(),
        };
        for ix in range.clone() {
            let peptide_index = PeptideIx(ix as u32);
            for mz in ions(ix) {
                let slot = &mut slots[bins.of(mz)];
                assert!(slot.0 < slot.1, "fragment count changed between passes");
                // SAFETY: `slot.0 < total` lies in group g's own range of this bin; the
                // ranges of all groups and bins are disjoint and tile `0..total`, which
                // is within the capacity
                unsafe {
                    out.0.add(slot.0).write(Theoretical {
                        peptide_index,
                        fragment_mz: mz,
                    })
                };
                slot.0 += 1;
            }
        }
        slots.iter().all(|(next, end)| next == end)
    });
    assert!(filled, "fragment count changed between passes");
    // SAFETY: every range was filled completely (asserted above) and the ranges tile
    // `0..total`, so all `total` elements are initialised
    unsafe { fragments.set_len(total) };

    // pass 3: order every bin (each arrived in item order)
    let mut rest = &mut fragments[..];
    let mut todo = Vec::with_capacity(nbins);
    for b in 0..nbins {
        let (bin, tail) = std::mem::take(&mut rest).split_at_mut(bounds[b + 1] - bounds[b]);
        rest = tail;
        if bin.len() > 1 {
            todo.push((b, bin));
        }
    }
    todo.into_par_iter()
        .for_each_init(SortScratch::default, |scratch, (b, bin)| {
            let low_bits = (!bins.is_overflow(b)).then_some(bins.shift);
            sort_bin(bin, low_bits, scratch);
        });
    fragments
}

/// At most `max_groups` contiguous ranges of whole `chunk_size` chunks of `0..n` (n > 0),
/// with about equal total `weight`, covering `0..n`
fn balanced_groups<W>(
    n: usize,
    chunk_size: usize,
    max_groups: usize,
    weight: W,
) -> Vec<std::ops::Range<usize>>
where
    W: Fn(usize) -> usize + Sync,
{
    let chunk_size = chunk_size.max(1);
    let weights: Vec<usize> = (0..n.div_ceil(chunk_size))
        .into_par_iter()
        .map(|c| {
            (c * chunk_size..((c + 1) * chunk_size).min(n))
                .map(&weight)
                .sum()
        })
        .collect();
    let ngroups = max_groups.clamp(1, weights.len());
    let per_group = weights.iter().sum::<usize>().div_ceil(ngroups).max(1);
    let mut groups = Vec::with_capacity(ngroups);
    let (mut start, mut acc) = (0, 0);
    for (c, w) in weights.iter().enumerate() {
        acc += w;
        if acc >= per_group || c + 1 == weights.len() {
            groups.push(start * chunk_size..((c + 1) * chunk_size).min(n));
            start = c + 1;
            acc = 0;
        }
    }
    groups
}

/// Total order of [`binned_fragments`]: m/z key, then peptide index
#[inline(always)]
fn mz_order(frag: &Theoretical) -> u64 {
    (mz_key(frag.fragment_mz) as u64) << 32 | frag.peptide_index.0 as u64
}

/// Puts a bin of [`binned_fragments`], which arrives in peptide order, in (m/z, peptide
/// index) order. `low_bits`: the number of low key bits that the bin does not fix, or
/// `None` for an overflow bin.
fn sort_bin(bin: &mut [Theoretical], low_bits: Option<u32>, scratch: &mut SortScratch) {
    match low_bits {
        // one m/z per bin: peptide order is the order
        Some(0) => {}
        // stable counting sort on the free key bits: equal m/z keep peptide order
        Some(bits) if bits <= 16 && bin.len() > SMALL_BIN && bin.len() <= u32::MAX as usize => {
            let mask = (1u32 << bits) - 1;
            let counts = &mut scratch.counts;
            counts.clear();
            counts.resize(1 << bits, 0);
            for frag in bin.iter() {
                counts[(mz_key(frag.fragment_mz) & mask) as usize] += 1;
            }
            let mut sum = 0;
            for c in counts.iter_mut() {
                sum += std::mem::replace(c, sum);
            }
            let items = &mut scratch.items;
            items.clear();
            items.extend_from_slice(bin);
            for frag in items.iter() {
                let c = &mut counts[(mz_key(frag.fragment_mz) & mask) as usize];
                bin[*c as usize] = *frag;
                *c += 1;
            }
        }
        // the full key decides; the arrival order does not matter
        _ => bin.sort_unstable_by_key(mz_order),
    }
}

/// Number of bits of the largest index below `n`
fn index_bits(n: usize) -> u32 {
    usize::BITS - n.saturating_sub(1).leading_zeros()
}

/// Stable LSD radix sort of `page` by peptide index (all indices below `2^bits`), in at
/// most 12-bit digits. A page holds `bucket_size` fragments in m/z order with random
/// peptide indices; a comparison sort of each page cost ~20 ns per fragment.
fn sort_by_peptide(page: &mut [Theoretical], bits: u32, scratch: &mut SortScratch) {
    if page.len() < 2 || bits == 0 {
        return;
    }
    if page.len() <= 64 || page.len() > u32::MAX as usize {
        page.sort_by_key(|frag| frag.peptide_index);
        return;
    }
    let passes = bits.min(32).div_ceil(12);
    let width = bits.min(32).div_ceil(passes);
    let mask = (1u32 << width) - 1;
    let counts = &mut scratch.counts;
    counts.clear();
    counts.resize(1 << width, 0);
    let items = &mut scratch.items;
    items.clear();
    items.extend_from_slice(page);
    let (mut src, mut dst): (&mut [Theoretical], &mut [Theoretical]) = (&mut items[..], page);
    for pass in 0..passes {
        let shift = pass * width;
        counts.iter_mut().for_each(|c| *c = 0);
        for frag in src.iter() {
            counts[((frag.peptide_index.0 >> shift) & mask) as usize] += 1;
        }
        let mut sum = 0;
        for c in counts.iter_mut() {
            sum += std::mem::replace(c, sum);
        }
        for frag in src.iter() {
            let c = &mut counts[((frag.peptide_index.0 >> shift) & mask) as usize];
            dst[*c as usize] = *frag;
            *c += 1;
        }
        std::mem::swap(&mut src, &mut dst);
    }
    // after an even number of passes the result is in `items`
    if passes % 2 == 0 {
        dst.copy_from_slice(src);
    }
}

/// Unsigned key whose order is that of [`f32::total_cmp`]
#[inline(always)]
fn mz_key(x: f32) -> u32 {
    let b = x.to_bits();
    if b >> 31 == 1 {
        !b
    } else {
        b | 0x8000_0000
    }
}

/// Stride of the per-page skip table
const SKIP: usize = 64;

/// Peptide index of every `SKIP`-th fragment of each page, `bucket_size / SKIP` entries
/// per page (padded with `u32::MAX` for the last, shorter page)
fn page_skip(fragments: &[Theoretical], bucket_size: usize) -> Vec<u32> {
    let per_page = bucket_size.div_ceil(SKIP);
    fragments
        .par_chunks(bucket_size)
        .flat_map_iter(|page| {
            (0..per_page).map(move |k| {
                page.get(k * SKIP)
                    .map_or(u32::MAX, |frag| frag.peptide_index.0)
            })
        })
        .collect()
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
}

#[derive(Default)]
pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    pub fragments: Vec<Theoretical>,
    pub ion_kinds: Vec<Kind>,
    pub min_value: Vec<f32>,
    /// Peptide index of every `SKIP`-th fragment of each page (see [`page_skip`])
    pub page_skip: Vec<u32>,
    /// Keep a list of potential (AA, mass) modifications for RT prediction
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
}

impl IndexedDatabase {
    /// Absolute fragment indices `[start, end)` within `page` that cover every fragment
    /// whose peptide index lies in `lo..=hi` (plus at most one fragment below `lo`, as
    /// with [`binary_search_slice`]; callers filter exact bounds anyway).
    ///
    /// Pages are sorted by peptide index. A plain binary search over an 8192-entry page
    /// costs ~13 dependent probes spread over 64 KB, nearly all cache and TLB misses; this
    /// is the dominant cost of the search. The skip table (4 bytes per `SKIP` fragments,
    /// ~18 MB for human tryptic) first narrows the search to a few `SKIP`-sized blocks.
    pub fn page_range(&self, page: usize, lo: usize, hi: usize) -> (usize, usize) {
        let (s, e) = self.page_blocks(page, lo, hi);
        let (l, r) = binary_search_slice(
            &self.fragments[s..e],
            |frag, bound| (frag.peptide_index.0 as usize).cmp(bound),
            lo,
            hi,
        );
        (s + l, s + r)
    }

    /// Skip-table entries of `page`
    pub fn page_skip(&self, page: usize) -> &[u32] {
        let per_page = self.bucket_size.div_ceil(SKIP);
        self.page_skip
            .get(page * per_page..(page + 1) * per_page)
            .unwrap_or_default()
    }

    /// First step of [`Self::page_range`]: the `SKIP`-sized blocks of `page` that can
    /// hold peptide indices `lo..=hi`
    pub fn page_blocks(&self, page: usize, lo: usize, hi: usize) -> (usize, usize) {
        let start = page * self.bucket_size;
        let end = ((page + 1) * self.bucket_size).min(self.fragments.len());
        let per_page = self.bucket_size.div_ceil(SKIP);
        let (s, e) = match self.page_skip.get(page * per_page..(page + 1) * per_page) {
            Some(skip) => {
                // blocks before `first` hold only indices < lo; blocks from `last` on hold
                // only indices > hi
                let first = skip
                    .partition_point(|&v| (v as usize) < lo)
                    .saturating_sub(1);
                let last = skip.partition_point(|&v| (v as usize) <= hi);
                (
                    (start + first * SKIP).min(end),
                    (start + last * SKIP).min(end),
                )
            }
            None => (start, end),
        };
        (s, e)
    }

    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
    pub fn query(
        &self,
        precursor_mass: f32,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
    ) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.monoisotopic.total_cmp(bounds),
            precursor_lo,
            precursor_hi,
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            fragment_tol,
            pre_idx_lo,
            pre_idx_hi,
        }
    }

    pub fn size(&self) -> usize {
        self.fragments.len()
    }

    pub fn buckets(&self) -> &[f32] {
        &self.min_value
    }

    pub fn serialize(&self) {
        use std::io::Write;
        let mut wtr = std::io::BufWriter::new(std::fs::File::create("fragments.bin").unwrap());
        for fragment in &self.fragments {
            let _ = wtr.write(&fragment.fragment_mz.to_le_bytes()).unwrap();
            let _ = wtr.write(&fragment.peptide_index.0.to_le_bytes()).unwrap();
        }
        wtr.flush().unwrap();

        let mut wtr = std::io::BufWriter::new(std::fs::File::create("peptides.csv").unwrap());
        writeln!(wtr, "peptide,proteins,monoisotopic,decoy").unwrap();
        for fragment in &self.peptides {
            writeln!(
                wtr,
                "{},{},{},{}",
                fragment,
                fragment.proteins(&self.decoy_tag, self.generate_decoys),
                fragment.monoisotopic,
                fragment.decoy
            )
            .unwrap();
        }
        wtr.flush().unwrap();
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        &self.peptides[index.0 as usize]
    }
}

pub struct IndexedQuery<'d> {
    db: &'d IndexedDatabase,
    precursor_mass: f32,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, mass: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices that can be used with `self.db.fragments`
        let (left_idx, right_idx) = binary_search_slice(
            &self.db.min_value,
            |min, bounds| min.total_cmp(bounds),
            fragment_lo,
            fragment_hi,
        );

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of peptide_index (precursor m/z)
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.fragments.len());

            // Narrow down into our region of interest: the slice of matching precursor mzs
            let slice = &&self.db.fragments[left_idx..right_idx];
            let (inner_left, inner_right) = {
                let (l, r) = self.db.page_range(page, self.pre_idx_lo, self.pre_idx_hi);
                (l - left_idx, r - left_idx)
            };

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                // This looks somewhat complicated, but it's a consequence of
                // how the `binary_search_slice` function works - it will return
                // the set of indices that maximally cover the desired range - the exact
                // `left` and `right` indices may be valid, or just outside of the range.
                // Anything interior of `left` and `right` is guaranteed to be within the
                // precursor tolerance, so we just need to check the edge cases
                //
                // Previously, a direct lookup to check the mass of the current fragment was
                // performed, but the pointer indirection + float comparison can slow down
                // open searches by as much as 2x!!
                // e.g. used to be `self.db[frag.peptide_index].monoisotopic >= precursor_lo`
                (frag.peptide_index.0 > self.pre_idx_lo as u32
                    || (frag.peptide_index.0 == self.pre_idx_lo as u32
                        && self.db[frag.peptide_index].monoisotopic >= precursor_lo))
                    && (frag.peptide_index.0 < self.pre_idx_hi as u32
                        || (frag.peptide_index.0 == self.pre_idx_hi as u32
                            && self.db[frag.peptide_index].monoisotopic <= precursor_hi))
                    && frag.fragment_mz >= fragment_lo
                    && frag.fragment_mz <= fragment_hi
            })
        })
    }
}

/// Ask the CPU to start loading `data` into cache (no-op off x86-64)
#[inline(always)]
pub fn prefetch<T>(data: &[T]) {
    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
        let p = data.as_ptr() as *const i8;
        for offset in (0..std::mem::size_of_val(data)).step_by(64) {
            // SAFETY: prefetching is a hint and never faults; the address is in `data`
            unsafe { _mm_prefetch::<_MM_HINT_T0>(p.add(offset)) };
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = data;
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] > high || right == slice.len()`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = slice
        .partition_point(|a| key(a, &low) == Ordering::Less)
        .saturating_sub(1);

    let right_idx =
        slice[left_idx..].partition_point(|a| key(a, &high) != Ordering::Greater) + left_idx;

    (left_idx, right_idx)
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use quickcheck_macros::quickcheck;

    use super::*;

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }

    #[test]
    fn structured_variable_mod_config_round_trips() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "variable_mods": {
                "M": [15.9949],
                "K": [
                    {"mass": 42.0106, "max_count": 1},
                    {"mass": 14.0157}
                ]
            },
            "max_variable_mods": 2,
            "max_combinations": 0
        }))
        .unwrap();

        let params = builder.make_parameters();
        assert_eq!(params.max_variable_mods, 2);
        assert_eq!(params.max_combinations, Some(1));

        let mods = params.variable_modifications();
        assert_eq!(mods.len(), 3);
        assert_eq!(mods[0].0, ModificationSpecificity::Residue(b'K'));
        assert!((mods[0].1 - 42.0106).abs() < 1e-4);
        assert_eq!(mods[0].2, Some(1));
        assert_eq!(mods[1].0, ModificationSpecificity::Residue(b'K'));
        assert!((mods[1].1 - 14.0157).abs() < 1e-4);
        assert_eq!(mods[1].2, None);
        assert_eq!(mods[2].0, ModificationSpecificity::Residue(b'M'));
        assert!((mods[2].1 - 15.9949).abs() < 1e-4);
        assert_eq!(mods[2].2, None);

        let serialized = serde_json::to_value(params).unwrap();
        let k_entries = &serialized["variable_mods"]["K"];
        assert!(k_entries[0].is_object());
        assert_eq!(k_entries[0]["max_count"], 1);
        assert!(k_entries[1].is_object());
        assert!(k_entries[1].get("max_count").is_none());
        assert!(serialized["variable_mods"]["M"][0].is_number());
    }

    #[test]
    fn decoy_colliding_with_target_is_dropped_not_merged() {
        let peptide = |protein: &str, decoy| {
            Peptide::try_from(crate::enzyme::Digest {
                decoy,
                sequence: "AELDWGK".into(),
                protein: protein.into(),
                ..Default::default()
            })
            .unwrap()
        };
        // e.g. the target from one prefilter chunk and a decoy from another
        let mut peptides = vec![peptide("P2", true), peptide("P1", false)];
        Parameters::reorder_peptides(&mut peptides);
        assert_eq!(peptides.len(), 1);
        assert!(!peptides[0].decoy);
        assert_eq!(peptides[0].proteins, vec![Arc::from("P1")]);
    }

    #[test]
    fn page_range_covers_the_same_fragments_as_a_full_page_search() {
        for bucket_size in [32, 64, 128, 1000, 8192] {
            let mut builder = Builder {
                bucket_size: Some(bucket_size),
                enzyme: Some(EnzymeBuilder {
                    missed_cleavages: Some(2),
                    min_len: Some(5),
                    ..Default::default()
                }),
                ..Default::default()
            };
            builder.update_fasta("unused".into());
            let fasta = Fasta::parse(
                include_str!("../../../tests/Q99536.fasta").into(),
                "rev_",
                true,
            );
            let db = builder.make_parameters().build(fasta);
            let n = db.peptides.len();
            let pages = db.fragments.len().div_ceil(db.bucket_size);
            for lo in (0..n).step_by(7) {
                for hi in [lo, lo + 1, lo + 5, lo + 40, n + 3] {
                    for page in 0..pages {
                        let start = page * db.bucket_size;
                        let end = (start + db.bucket_size).min(db.fragments.len());
                        let page_frags = &db.fragments[start..end];
                        let expected = page_frags
                            .iter()
                            .filter(|f| (lo..=hi).contains(&(f.peptide_index.0 as usize)))
                            .count();
                        let (l, r) = db.page_range(page, lo, hi);
                        assert!(start <= l && l <= r && r <= end);
                        let got = db.fragments[l..r]
                            .iter()
                            .filter(|f| (lo..=hi).contains(&(f.peptide_index.0 as usize)))
                            .count();
                        assert_eq!(
                            got, expected,
                            "bucket {bucket_size} page {page} {lo}..={hi}"
                        );
                    }
                }
            }
        }
    }

    /// The fragments of `peptides` in generation (peptide) order, the order in which the
    /// index build used to write them before its global sort
    fn fragments_in_peptide_order(params: &Parameters, peptides: &[Peptide]) -> Vec<Theoretical> {
        peptides
            .iter()
            .enumerate()
            .flat_map(|(idx, peptide)| {
                params.index_ions(peptide).map(move |ion| Theoretical {
                    peptide_index: PeptideIx(idx as u32),
                    fragment_mz: ion.monoisotopic_mass,
                })
            })
            .collect()
    }

    /// The binned build's order, by definition: a stable sort by m/z of the fragments in
    /// peptide order
    fn stable_mz_sort(mut fragments: Vec<Theoretical>) -> Vec<Theoretical> {
        fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));
        fragments
    }

    /// Bit patterns: `PartialEq` of `f32` equates -0.0 and 0.0 and never NaN
    fn bits(fragments: &[Theoretical]) -> Vec<(u32, u32)> {
        fragments
            .iter()
            .map(|f| (f.peptide_index.0, f.fragment_mz.to_bits()))
            .collect()
    }

    fn q99536_parameters(bucket_size: usize) -> Parameters {
        let mut builder = Builder {
            bucket_size: Some(bucket_size),
            enzyme: Some(EnzymeBuilder {
                missed_cleavages: Some(2),
                min_len: Some(5),
                ..Default::default()
            }),
            variable_mods: Some(
                [("M".to_string(), vec![VarModEntry::Mass(15.9949)])]
                    .into_iter()
                    .collect(),
            ),
            ..Default::default()
        };
        builder.update_fasta("unused".into());
        builder.make_parameters()
    }

    fn q99536() -> Fasta {
        Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        )
    }

    #[test]
    fn binned_fragments_are_the_stable_mz_sort() {
        let params = q99536_parameters(8192);
        let peptides = params.digest(&q99536());
        let expected = bits(&stable_mz_sort(fragments_in_peptide_order(
            &params, &peptides,
        )));
        assert!(expected.len() > 1000);
        let ties = expected.windows(2).filter(|w| w[0].1 == w[1].1).count();
        assert!(ties > 0, "the test needs equal fragment m/z");
        // 2 bins: wide key range per bin (comparison sort); MAX_BINS: small bins
        for max_bins in [2, 3, 64, 1000, MAX_BINS] {
            for chunk_size in [1, 7, 4096] {
                for max_groups in [1, 3, 256] {
                    let got =
                        params.sorted_fragments_with(&peptides, chunk_size, max_bins, max_groups);
                    assert_eq!(bits(&got), expected, "{max_bins} {chunk_size} {max_groups}");
                }
            }
        }
        // the default grouping depends on the thread count, the result must not
        for threads in [1, 2, 7] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let got = pool.install(|| params.sorted_fragments(&peptides));
            assert_eq!(bits(&got), expected, "{threads} threads");
        }
        assert!(params.sorted_fragments(&[]).is_empty());

        // fragments outside the guessed range (a1 ions below GUESS_LOW; masses above a
        // peptide mass that is too low) are sorted in the overflow bins
        let mut wide = params.clone();
        wide.ion_kinds = vec![Kind::A, Kind::B, Kind::Y];
        wide.min_ion_index = 0;
        let expected = bits(&stable_mz_sort(fragments_in_peptide_order(
            &wide, &peptides,
        )));
        assert!(f32::from_bits(expected[0].1) < GUESS_LOW);
        for max_bins in [2, 64, MAX_BINS] {
            let got = wide.sorted_fragments_with(&peptides, 7, max_bins, 5);
            assert_eq!(bits(&got), expected);
        }
        let mut light = peptides.clone();
        light.iter_mut().for_each(|p| p.monoisotopic = 100.0);
        let expected = bits(&stable_mz_sort(fragments_in_peptide_order(&params, &light)));
        assert_eq!(
            bits(&params.sorted_fragments_with(&light, 7, MAX_BINS, 5)),
            expected
        );
    }

    /// Random fragments for [`binned_fragments`]: many equal m/z, m/z concentrated in a
    /// few bins (counting sort), special values (signed zeros, infinities, NaN,
    /// subnormals, negative masses) and random bit patterns, with random bins and groups
    #[derive(Clone, Debug)]
    struct RandomFragments {
        items: Vec<Vec<f32>>,
        lowest: f32,
        highest: f32,
        max_bins: usize,
        chunk_size: usize,
        max_groups: usize,
    }

    impl quickcheck::Arbitrary for RandomFragments {
        fn arbitrary(g: &mut quickcheck::Gen) -> Self {
            let center = *g.choose(&[60.0f32, 499.9, 512.0, 1999.0]).unwrap();
            let special = [
                f32::NAN,
                -f32::NAN,
                0.0,
                -0.0,
                f32::INFINITY,
                f32::NEG_INFINITY,
                -5.0,
                1e-40,
                f32::MIN_POSITIVE,
                f32::MAX,
                49.99,
                50.0,
            ];
            let mz = |g: &mut quickcheck::Gen| match u8::arbitrary(g) % 8 {
                0 => *g.choose(&special).unwrap(),
                1 => f32::from_bits(u32::arbitrary(g)),
                2 => center + (u8::arbitrary(g) % 8) as f32 * 1e-3,
                _ => center + (u16::arbitrary(g) % 2000) as f32 * 1e-3,
            };
            let n = usize::arbitrary(g) % 400;
            let items = (0..n)
                .map(|_| (0..usize::arbitrary(g) % 16).map(|_| mz(g)).collect())
                .collect();
            RandomFragments {
                items,
                lowest: *g.choose(&[50.0, 0.0, -1.0, 500.0, f32::NAN]).unwrap(),
                highest: *g
                    .choose(&[f32::NEG_INFINITY, 100.0, 600.0, 2500.0, 1e9, f32::INFINITY])
                    .unwrap(),
                max_bins: *g.choose(&[0, 1, 2, 3, 5, 64, 1000, 8192, 65536]).unwrap(),
                chunk_size: 1 + usize::arbitrary(g) % 9,
                max_groups: usize::arbitrary(g) % 10,
            }
        }
    }

    #[quickcheck]
    fn binned_fragments_match_a_stable_mz_sort(input: RandomFragments) -> bool {
        let RandomFragments {
            items,
            lowest,
            highest,
            max_bins,
            chunk_size,
            max_groups,
        } = input;
        let expected: Vec<Theoretical> = items
            .iter()
            .enumerate()
            .flat_map(|(ix, mzs)| {
                mzs.iter().map(move |&mz| Theoretical {
                    peptide_index: PeptideIx(ix as u32),
                    fragment_mz: mz,
                })
            })
            .collect();
        let expected = bits(&stable_mz_sort(expected));
        let bins = MzBins::new(lowest, highest, max_bins);
        let weight = |ix: usize| items[ix].len();
        let got = binned_fragments(
            items.len(),
            weight,
            |ix| items[ix].iter().copied(),
            bins,
            chunk_size,
            max_groups,
        );
        // the order in which the m/z of one item arrive does not matter
        let reversed = binned_fragments(
            items.len(),
            weight,
            |ix| items[ix].iter().rev().copied(),
            bins,
            chunk_size,
            max_groups,
        );
        bits(&got) == expected && bits(&reversed) == expected
    }

    #[quickcheck]
    fn mz_bins_follow_the_key_order(
        a: u32,
        b: u32,
        max_bins: u16,
        lowest: f32,
        highest: f32,
    ) -> bool {
        let bins = MzBins::new(lowest, highest, max_bins as usize);
        let (a, b) = (f32::from_bits(a), f32::from_bits(b));
        let (ba, bb) = (bins.of(a), bins.of(b));
        let monotone = match mz_key(a).cmp(&mz_key(b)) {
            Ordering::Less => ba <= bb,
            Ordering::Equal => ba == bb,
            Ordering::Greater => ba >= bb,
        };
        // a regular bin fixes the key bits above `shift`
        let regular =
            ba != bb || bins.is_overflow(ba) || mz_key(a) >> bins.shift == mz_key(b) >> bins.shift;
        monotone && regular && ba < bins.len && bins.len <= (max_bins as usize).max(2) + 2
    }

    #[test]
    fn bin_sort_by_counting_equals_the_full_order() {
        // xorshift: deterministic pseudo-random input
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut scratch = SortScratch::default();
        for low_bits in [0, 1, 4, 13, 16, 17, 24] {
            for len in [2, 3, SMALL_BIN, SMALL_BIN + 1, 1000, 5000] {
                // all keys of the bin share the bits above `low_bits`; few distinct low
                // bits and peptides, so that both repeat
                let high = mz_key(700.0) >> low_bits << low_bits;
                let distinct = (next() % 64 + 1) as u32;
                let mut bin: Vec<Theoretical> = (0..len)
                    .map(|_| {
                        let low = if low_bits == 0 {
                            0
                        } else {
                            (next() as u32 % distinct) & ((1u32 << low_bits) - 1).max(1)
                        };
                        let key = high | low;
                        // inverse of `mz_key` for positive values
                        let mz = f32::from_bits(key & 0x7fff_ffff);
                        Theoretical {
                            peptide_index: PeptideIx((next() % 50) as u32),
                            fragment_mz: mz,
                        }
                    })
                    .collect();
                // a bin arrives in peptide order
                bin.sort_by_key(|f| f.peptide_index);
                let mut expected = bin.clone();
                expected.sort_unstable_by_key(mz_order);
                sort_bin(&mut bin, Some(low_bits), &mut scratch);
                assert_eq!(bits(&bin), bits(&expected), "{low_bits} bits, {len}");
                sort_bin(&mut expected, None, &mut scratch);
                assert_eq!(bits(&bin), bits(&expected), "{low_bits} bits, {len}");
            }
        }
    }

    #[test]
    fn page_radix_sort_is_a_stable_sort_by_peptide() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut scratch = SortScratch::default();
        for bits_ in [1, 5, 12, 13, 24, 25, 32] {
            for len in [2, 3, 64, 65, 1000, 8192] {
                let max = if bits_ == 32 {
                    u32::MAX
                } else {
                    (1u32 << bits_) - 1
                };
                // repeated peptide indices (several fragments of one peptide in a page),
                // and the largest index
                let mut page: Vec<Theoretical> = (0..len)
                    .map(|i| Theoretical {
                        peptide_index: PeptideIx(match next() % 4 {
                            0 => max,
                            1 => (next() % 8) as u32 & max,
                            _ => (next() as u32) & max,
                        }),
                        fragment_mz: i as f32,
                    })
                    .collect();
                let mut expected = page.clone();
                expected.sort_by_key(|f| f.peptide_index);
                sort_by_peptide(&mut page, bits_, &mut scratch);
                assert_eq!(bits(&page), bits(&expected), "bits {bits_} len {len}");
            }
        }
        assert_eq!(index_bits(0), 0);
        assert_eq!(index_bits(1), 0);
        assert_eq!(index_bits(2), 1);
        assert_eq!(index_bits(4096), 12);
        assert_eq!(index_bits(4097), 13);
    }

    /// The layout of a global unstable sort, which leaves fragments of equal m/z, and
    /// fragments of one peptide within a page, in arbitrary order: here the opposite of
    /// the binned build's (peptide index descending, m/z descending)
    fn opposite_tie_layout(params: &Parameters, peptides: &[Peptide]) -> IndexedDatabase {
        let mut fragments = fragments_in_peptide_order(params, peptides);
        fragments.sort_by(|a, b| {
            a.fragment_mz
                .total_cmp(&b.fragment_mz)
                .then(b.peptide_index.cmp(&a.peptide_index))
        });
        let min_value = fragments
            .chunks_mut(params.bucket_size)
            .map(|page| {
                let min = page[0].fragment_mz;
                page.sort_by(|a, b| {
                    a.peptide_index
                        .cmp(&b.peptide_index)
                        .then(b.fragment_mz.total_cmp(&a.fragment_mz))
                });
                min
            })
            .collect();
        IndexedDatabase {
            page_skip: page_skip(&fragments, params.bucket_size),
            peptides: peptides.to_vec(),
            fragments,
            min_value,
            bucket_size: params.bucket_size,
            ..Default::default()
        }
    }

    #[test]
    fn search_does_not_depend_on_the_order_of_ties() {
        // small pages, so that page boundaries fall between fragments of equal m/z
        for bucket_size in [1, 2, 4, 8, 64, 8192] {
            let params = q99536_parameters(bucket_size);
            assert_eq!(params.bucket_size, bucket_size);
            let fasta = q99536();
            let peptides = params.digest(&fasta);
            let db = params.clone().build(fasta);
            let other = opposite_tie_layout(&params, &peptides);
            assert_eq!(db.peptides.len(), other.peptides.len());
            assert_eq!(db.fragments.len(), other.fragments.len());
            // the minimum m/z of a page is an order statistic
            assert_eq!(
                db.min_value.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                other
                    .min_value
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>()
            );
            // pages hold the same m/z values, possibly of other peptides
            for (a, b) in db
                .fragments
                .chunks(bucket_size)
                .zip(other.fragments.chunks(bucket_size))
            {
                let mut a: Vec<u32> = a.iter().map(|f| mz_key(f.fragment_mz)).collect();
                let mut b: Vec<u32> = b.iter().map(|f| mz_key(f.fragment_mz)).collect();
                a.sort_unstable();
                b.sort_unstable();
                assert_eq!(a, b);
            }
            let differs = bits(&db.fragments) != bits(&other.fragments);
            assert!(differs || bucket_size == 8192, "the layouts should differ");

            // every query finds the same fragments; query m/z of existing fragments (ties)
            let tolerances = [
                (Tolerance::Ppm(-20.0, 20.0), Tolerance::Ppm(-10.0, 10.0)),
                (Tolerance::Da(-0.5, 0.5), Tolerance::Da(0.0, 0.0)),
                (Tolerance::Da(-150.0, 150.0), Tolerance::Da(-0.05, 0.05)),
            ];
            for (pi, peptide) in peptides.iter().enumerate().step_by(5) {
                for (precursor_tol, fragment_tol) in tolerances {
                    let qa = db.query(peptide.monoisotopic, precursor_tol, fragment_tol);
                    let qb = other.query(peptide.monoisotopic, precursor_tol, fragment_tol);
                    for frag in db.fragments.iter().skip(pi % 11).step_by(97) {
                        let mut a = bits(
                            &qa.page_search(frag.fragment_mz)
                                .copied()
                                .collect::<Vec<_>>(),
                        );
                        let mut b = bits(
                            &qb.page_search(frag.fragment_mz)
                                .copied()
                                .collect::<Vec<_>>(),
                        );
                        a.sort_unstable();
                        b.sort_unstable();
                        assert_eq!(
                            a, b,
                            "bucket {bucket_size} peptide {pi} m/z {}",
                            frag.fragment_mz
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn mz_key_orders_like_total_cmp() {
        let values = [
            f32::NEG_INFINITY,
            -2.5,
            -0.0,
            0.0,
            1e-30,
            150.0,
            150.00002,
            9000.0,
            f32::INFINITY,
            f32::NAN,
            -f32::NAN,
        ];
        for a in values {
            for b in values {
                assert_eq!(mz_key(a).cmp(&mz_key(b)), a.total_cmp(&b), "{a} {b}");
            }
        }
    }

    #[quickcheck]
    fn mz_key_orders_random_bits_like_total_cmp(a: u32, b: u32) -> bool {
        let (x, y) = (f32::from_bits(a), f32::from_bits(b));
        mz_key(x).cmp(&mz_key(y)) == x.total_cmp(&y)
    }

    #[test]
    fn digestion() {
        let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

        let fasta = Fasta::parse(fasta.into(), "rev_", false);

        // Make sure that FASTA parsed OK
        assert_eq!(
            fasta.targets,
            vec![
                (
                    Arc::from("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::from("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let params = Parameters {
            bucket_size: 128,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            static_mods: HashMap::default(),
            variable_mods: [(
                ModificationSpecificity::ProteinN(None),
                vec![VarModEntry::Mass(42.0)],
            )]
            .into_iter()
            .collect(),
            max_variable_mods: 2,
            max_combinations: None,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
        };

        let peptides = params.digest(&fasta);

        let expected = [
            "EQALLK",
            "LEQSMR",
            "AQLTQLK",
            "MEWKLEQSMR",
            "[+42]-MEWKLEQSMR",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

        let sequences = peptides.iter().map(|p| p.to_string()).collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }
}
