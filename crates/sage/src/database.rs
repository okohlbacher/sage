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

    /// All theoretical fragments, in peptide order, written into one exact-size
    /// allocation. Collecting a parallel `flat_map` instead goes through a
    /// `LinkedList<Vec<_>>` and a concatenation copy, so two copies of the largest
    /// structure in Sage are alive at once (6.7 GB instead of 3.0 GB peak before
    /// sorting for human tryptic, 282 M fragments).
    #[cfg(test)]
    fn fragments(&self, peptides: &[Peptide]) -> Vec<Theoretical> {
        self.fragments_in_blocks(peptides, 4096)
    }

    #[cfg(test)]
    fn fragments_in_blocks(&self, peptides: &[Peptide], block_size: usize) -> Vec<Theoretical> {
        let totals: Vec<usize> = peptides
            .par_chunks(block_size)
            .map(|chunk| chunk.iter().map(|p| self.index_ions(p).count()).sum())
            .collect();
        let total: usize = totals.iter().sum();

        let mut fragments: Vec<Theoretical> = Vec::with_capacity(total);
        advise_huge_pages(fragments.spare_capacity_mut());
        let mut rest = &mut fragments.spare_capacity_mut()[..total];
        let mut blocks = Vec::with_capacity(totals.len());
        for &n in &totals {
            let (head, tail) = std::mem::take(&mut rest).split_at_mut(n);
            blocks.push(head);
            rest = tail;
        }
        blocks
            .into_par_iter()
            .zip(peptides.par_chunks(block_size))
            .enumerate()
            .for_each(|(block, (out, chunk))| {
                let mut i = 0;
                for (offset, peptide) in chunk.iter().enumerate() {
                    for ion in self.index_ions(peptide) {
                        out[i].write(Theoretical {
                            peptide_index: PeptideIx((block * block_size + offset) as u32),
                            fragment_mz: ion.monoisotopic_mass,
                        });
                        i += 1;
                    }
                }
                assert_eq!(i, out.len(), "fragment count changed between passes");
            });
        // SAFETY: the blocks tile `0..total` and each wrote exactly `out.len()`
        // elements (asserted above), so all `total` elements are initialised.
        unsafe { fragments.set_len(total) };
        fragments
    }

    /// All theoretical fragments sorted by m/z (`f32::total_cmp`; equal m/z in peptide
    /// order), written into one exact-size allocation without a comparison sort.
    ///
    /// The fragments are counted per m/z bin first, then generated a second time and
    /// written straight into their bin (peptide order within a bin), and each bin is
    /// put in order with a counting sort on the key bits that the bin does not fix.
    /// A global `par_sort_unstable` of the ~3*10^8 fragments of a large search space
    /// moved the 2.4 GB array ~log2(n) times through memory. The order of fragments
    /// with equal m/z can differ from that sort's; the search does not depend on it
    /// (every lookup visits all fragments within its m/z and precursor bounds).
    fn sorted_fragments(&self, peptides: &[Peptide]) -> Vec<Theoretical> {
        self.sorted_fragments_with(peptides, 4096, MAX_BINS)
    }

    fn sorted_fragments_with(
        &self,
        peptides: &[Peptide],
        chunk_size: usize,
        max_bins: usize,
    ) -> Vec<Theoretical> {
        // pass 1: fragment count and m/z key range of every chunk of peptides
        let stats: Vec<(usize, u32, u32)> = peptides
            .par_chunks(chunk_size)
            .map(|chunk| {
                let (mut n, mut lo, mut hi) = (0usize, u32::MAX, 0u32);
                for peptide in chunk {
                    for ion in self.index_ions(peptide) {
                        let k = mz_key(ion.monoisotopic_mass);
                        n += 1;
                        lo = lo.min(k);
                        hi = hi.max(k);
                    }
                }
                (n, lo, hi)
            })
            .collect();
        let total: usize = stats.iter().map(|s| s.0).sum();
        if total == 0 {
            return Vec::new();
        }
        let lo = stats.iter().map(|s| s.1).min().unwrap_or(0);
        let hi = stats.iter().map(|s| s.2).max().unwrap_or(0);
        let mut shift = 0;
        while ((hi >> shift) - (lo >> shift)) as usize >= max_bins.max(2) {
            shift += 1;
        }
        let base = lo >> shift;
        let nbins = ((hi >> shift) - base) as usize + 1;
        let bin = |mz: f32| ((mz_key(mz) >> shift) - base) as usize;

        // groups of consecutive chunks with about the same number of fragments
        let ngroups = (rayon::current_num_threads() * 4)
            .clamp(1, 256)
            .min(stats.len());
        let per_group = total.div_ceil(ngroups);
        let mut groups: Vec<std::ops::Range<usize>> = Vec::with_capacity(ngroups);
        let (mut start, mut acc) = (0, 0);
        for (ix, s) in stats.iter().enumerate() {
            acc += s.0;
            if acc >= per_group || ix + 1 == stats.len() {
                groups.push(start * chunk_size..((ix + 1) * chunk_size).min(peptides.len()));
                start = ix + 1;
                acc = 0;
            }
        }

        // pass 2: per-group histograms over the bins
        let hist: Vec<Vec<usize>> = groups
            .par_iter()
            .map(|range| {
                let mut h = vec![0usize; nbins];
                for peptide in &peptides[range.clone()] {
                    for ion in self.index_ions(peptide) {
                        h[bin(ion.monoisotopic_mass)] += 1;
                    }
                }
                h
            })
            .collect();
        // bin-major, then group order: within a bin, fragments stay in peptide order.
        // `slots[g][b]` = (next write position, end) of group g's range in bin b.
        let mut bounds = Vec::with_capacity(nbins + 1);
        let mut slots = vec![vec![(0usize, 0usize); nbins]; groups.len()];
        let mut pos = 0;
        for b in 0..nbins {
            bounds.push(pos);
            for (g, h) in hist.iter().enumerate() {
                slots[g][b] = (pos, pos + h[b]);
                pos += h[b];
            }
        }
        bounds.push(pos);
        drop(hist);
        assert_eq!(pos, total, "fragment count changed between passes");

        // pass 3: write every fragment into its bin
        let mut fragments: Vec<Theoretical> = Vec::with_capacity(total);
        advise_huge_pages(fragments.spare_capacity_mut());
        struct Out(*mut Theoretical);
        // SAFETY: the groups write disjoint, bounds-checked index ranges (their slots)
        unsafe impl Sync for Out {}
        let out = Out(fragments.as_mut_ptr());
        let out = &out;
        groups
            .par_iter()
            .zip(slots.par_iter_mut())
            .for_each(|(range, slots)| {
                for (ix, peptide) in peptides[range.clone()].iter().enumerate() {
                    let peptide_index = PeptideIx((range.start + ix) as u32);
                    for ion in self.index_ions(peptide) {
                        let slot = &mut slots[bin(ion.monoisotopic_mass)];
                        assert!(slot.0 < slot.1, "fragment count changed between passes");
                        // SAFETY: `slot.0` lies in this group's own range of the bin, and
                        // the ranges of all groups and bins tile `0..total`
                        unsafe {
                            out.0.add(slot.0).write(Theoretical {
                                peptide_index,
                                fragment_mz: ion.monoisotopic_mass,
                            })
                        };
                        slot.0 += 1;
                    }
                }
            });
        assert!(
            slots.iter().flatten().all(|(next, end)| next == end),
            "fragment count changed between passes"
        );
        // SAFETY: every slot range was filled completely (asserted above) and the slot
        // ranges tile `0..total`, so all `total` elements are initialised
        unsafe { fragments.set_len(total) };

        // pass 4: order each bin by the key bits below `shift` (stable: equal m/z keep
        // peptide order)
        let mut bins = Vec::with_capacity(nbins);
        let mut rest = &mut fragments[..];
        for b in 0..nbins {
            let (head, tail) = std::mem::take(&mut rest).split_at_mut(bounds[b + 1] - bounds[b]);
            bins.push(head);
            rest = tail;
        }
        let low_mask = if shift == 0 { 0 } else { (1u32 << shift) - 1 };
        bins.into_par_iter().for_each(|slice| {
            if slice.len() < 2 || shift == 0 {
                return;
            }
            if shift <= 16 {
                let mut count = vec![0usize; 1 << shift];
                for f in slice.iter() {
                    count[(mz_key(f.fragment_mz) & low_mask) as usize] += 1;
                }
                let mut sum = 0;
                for c in count.iter_mut() {
                    let n = *c;
                    *c = sum;
                    sum += n;
                }
                let scratch = slice.to_vec();
                for f in scratch {
                    let k = (mz_key(f.fragment_mz) & low_mask) as usize;
                    slice[count[k]] = f;
                    count[k] += 1;
                }
            } else {
                slice.sort_by_key(|f| mz_key(f.fragment_mz));
            }
        });
        fragments
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

        let peptide_bits = u32::BITS - (target_decoys.len().max(1) as u32 - 1).leading_zeros();
        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map_init(Vec::new, |scratch, chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                sort_by_peptide(chunk, scratch, peptide_bits);
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

/// Upper bound on the number of m/z bins of [`Parameters::sorted_fragments`]
const MAX_BINS: usize = 8192;

/// Stable LSD radix sort of `page` by peptide index (indices below `2^bits`), in at
/// most 12-bit digits. A page holds `bucket_size` fragments in m/z order with random
/// peptide indices; a comparison sort of each page cost ~20 ns per fragment.
fn sort_by_peptide(page: &mut [Theoretical], scratch: &mut Vec<Theoretical>, bits: u32) {
    if page.len() < 2 || bits == 0 {
        return;
    }
    if page.len() <= 64 {
        page.sort_by_key(|frag| frag.peptide_index);
        return;
    }
    let passes = bits.div_ceil(12);
    let width = bits.div_ceil(passes);
    let mask = (1u32 << width) - 1;
    let mut count = vec![0usize; 1 << width];
    scratch.clear();
    scratch.extend_from_slice(page);
    let (mut src, mut dst): (&mut [Theoretical], &mut [Theoretical]) = (&mut scratch[..], page);
    for pass in 0..passes {
        let shift = pass * width;
        count.iter_mut().for_each(|c| *c = 0);
        for frag in src.iter() {
            count[((frag.peptide_index.0 >> shift) & mask) as usize] += 1;
        }
        let mut sum = 0;
        for c in count.iter_mut() {
            let n = *c;
            *c = sum;
            sum += n;
        }
        for frag in src.iter() {
            let d = ((frag.peptide_index.0 >> shift) & mask) as usize;
            dst[count[d]] = *frag;
            count[d] += 1;
        }
        std::mem::swap(&mut src, &mut dst);
    }
    // after an even number of passes the result is in `scratch`
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

    #[test]
    fn exact_size_fragments_match_collect() {
        let mut builder = Builder {
            enzyme: Some(EnzymeBuilder {
                missed_cleavages: Some(2),
                min_len: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        };
        builder.update_fasta("unused".into());
        let params = builder.make_parameters();
        let fasta = Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        );
        let peptides = params.digest(&fasta);
        // the previous implementation
        let expected = peptides
            .iter()
            .enumerate()
            .flat_map(|(idx, peptide)| {
                params.index_ions(peptide).map(move |ion| Theoretical {
                    peptide_index: PeptideIx(idx as u32),
                    fragment_mz: ion.monoisotopic_mass,
                })
            })
            .collect::<Vec<_>>();
        assert!(expected.len() > 1000);
        // small blocks so that the block tiling is exercised
        for block_size in [1, 7, 4096] {
            assert_eq!(params.fragments_in_blocks(&peptides, block_size), expected);
        }
    }

    #[test]
    fn binned_fragments_are_the_stable_mz_sort() {
        let mut builder = Builder {
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
        let params = builder.make_parameters();
        let fasta = Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        );
        let peptides = params.digest(&fasta);
        let mut expected = params.fragments(&peptides);
        // equal m/z: generation (peptide) order
        expected.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));
        assert!(expected.len() > 1000);
        let ties = expected
            .windows(2)
            .filter(|w| w[0].fragment_mz == w[1].fragment_mz)
            .count();
        assert!(ties > 0, "the test needs equal fragment m/z");
        // 2 bins: wide key range per bin (comparison-sort fallback)
        for max_bins in [2, 3, 64, 1000, MAX_BINS] {
            for chunk_size in [1, 7, 4096] {
                let got = params.sorted_fragments_with(&peptides, chunk_size, max_bins);
                assert_eq!(got, expected, "max_bins {max_bins} chunk {chunk_size}");
            }
        }
        assert!(params.sorted_fragments(&[]).is_empty());

        // pages: radix sort by peptide index = stable sort
        let mut scratch = Vec::new();
        for bits in [1, 5, 12, 13, 24, 25, 32] {
            for len in [2, 3, 64, 65, 1000, 8192] {
                let mut page: Vec<Theoretical> = (0..len)
                    .map(|i| Theoretical {
                        peptide_index: PeptideIx(
                            ((i as u64 * 2654435761) % (1u64 << bits)) as u32
                                & (u32::MAX >> (32 - bits)),
                        ),
                        fragment_mz: i as f32,
                    })
                    .collect();
                let mut expected = page.clone();
                expected.sort_by_key(|f| f.peptide_index);
                sort_by_peptide(&mut page, &mut scratch, bits);
                assert_eq!(page, expected, "bits {bits} len {len}");
            }
        }

        // keys order like `total_cmp`, including signed zeros, infinities and NaN
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
        ];
        for a in values {
            for b in values {
                assert_eq!(mz_key(a).cmp(&mz_key(b)), a.total_cmp(&b), "{a} {b}");
            }
        }
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
