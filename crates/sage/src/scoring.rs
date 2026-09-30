use crate::database::{binary_search_slice, prefetch, IndexedDatabase, PeptideIx};
use crate::heap::bounded_min_heapify;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::spectrum::{Precursor, ProcessedSpectrum};
use serde::{Deserialize, Serialize};
use std::ops::AddAssign;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum ScoreType {
    SageHyperScore,
    OpenMSHyperScore,
}

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug, PartialEq)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    longest_b: usize,
    longest_y: usize,
    hyperscore: f64,
    ppm_difference: f32,
    precursor_charge: u8,
    isotope_error: i8,
}

impl Eq for Score {}

// Must agree with `Ord`: `bounded_min_heapify` compares with `<`/`>`, i.e. `PartialOrd`.
// A derived `PartialOrd` compared `peptide` (the first field) instead of hyperscore, so the
// low-memory prefilter kept the candidates with the highest peptide index, not the best.
impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hyperscore
            .partial_cmp(&other.hyperscore)
            .unwrap_or(std::cmp::Ordering::Less)
    }
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
}

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

#[derive(Serialize, Clone, Debug, Default)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    // psm_id help to match with matched fragments table.
    pub psm_id: usize,
    pub peptide_len: usize,
    /// Spectrum id
    pub spec_id: String,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
    /// Ion mobility
    pub ims: f32,
    /// Predicted ion mobility, if enabled
    pub predicted_ims: f32,
    /// Difference between predicted & observed ion mobility
    pub delta_ims_model: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_next: f64,
    /// Difference between hyperscore of this candidate, and the best candidate
    pub delta_best: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,
    pub protein_group_q: f32,

    pub ms2_intensity: f32,

    pub protein_groups: Option<String>,
    pub num_protein_groups: u32,

    pub fragments: Option<Fragments>,
}

/// Matching Fragment details
#[derive(Serialize, Default, Clone, Debug)]
pub struct Fragments {
    /// Observed fragment charge state.
    #[serde(skip_serializing)]
    pub charges: Vec<i32>,
    pub kinds: Vec<Kind>,
    pub fragment_ordinals: Vec<i32>,
    pub intensities: Vec<f32>,
    pub mz_calculated: Vec<f32>,
    pub mz_experimental: Vec<f32>,
}

static PSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

fn increment_psm_counter() -> usize {
    PSM_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl ScoreType {
    pub fn score(&self, matched_b: u16, matched_y: u16, summed_b: f32, summed_y: f32) -> f64 {
        let score = match self {
            // Calculate the X!Tandem hyperscore
            Self::SageHyperScore => {
                let i = (summed_b + 1.0) as f64 * (summed_y + 1.0) as f64;

                i.ln() + lnfact(matched_b) + lnfact(matched_y)
            }
            // Calculate the OpenMS flavour hyperscore
            Self::OpenMSHyperScore => {
                let summed_intensity = summed_b + summed_y;

                summed_intensity.ln_1p() as f64 + lnfact(matched_b) + lnfact(matched_y)
            }
        };
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

impl Score {
    /// Calculate the hyperscore for a given PSM choosing between implementations based on `score_type`
    fn hyperscore(&self, score_type: ScoreType) -> f64 {
        score_type.score(self.matched_b, self.matched_y, self.summed_b, self.summed_y)
    }
}

pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub report_psms: usize,

    // Rather than use a fixed precursor tolerance, dynamically alter
    // the precursor tolerance window based on MS2 isolation window and charge
    pub wide_window: bool,
    pub annotate_matches: bool,
    pub score_type: ScoreType,
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

impl<'db> Scorer<'db> {
    /// Perform a quick first-pass scoring, where we consider a peptide "identified"
    /// if it meets the following criterion:
    ///  * prefilter_low_memory = true: in the top `report_psms` hits for a spectrum
    ///  * prefilter_low_memory = false: has at least `min_matched_peaks` fragment ion matches
    /// * `keep`: A vector of atomic bools is used to maintain an identification list across scans
    pub fn quick_score(
        &self,
        query: &ProcessedSpectrum,
        prefilter_low_memory: bool,
        keep: &[AtomicBool],
    ) {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });
        let hits = self.initial_hits(query, precursor);

        if prefilter_low_memory {
            let mut score_vector = hits
                .preliminary
                .iter()
                .filter_map(|pre| {
                    if pre.peptide == PeptideIx::default() {
                        return None;
                    }
                    let (score, _) = self.score_candidate(query, pre);
                    if (score.matched_b + score.matched_y) < self.min_matched_peaks {
                        return None;
                    }
                    Some(score)
                })
                .collect::<Vec<_>>();

            let k = self.report_psms.min(score_vector.len());
            bounded_min_heapify(&mut score_vector, k);
            for score in &score_vector[..k] {
                keep[score.peptide.0 as usize].store(true, Ordering::Relaxed);
            }
        } else {
            for pre in &hits.preliminary {
                if pre.peptide != PeptideIx::default() {
                    keep[pre.peptide.0 as usize].store(true, Ordering::Relaxed);
                }
            }
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }

    /// Perform a k-select and truncation of an [`InitialHits`] list.
    ///
    /// Determine how many candidates to actually calculate hyperscore for.
    /// Hyperscore is relatively computationally expensive, so we don't want
    /// to calculate it for every possible candidate (100s - 10,000s depending on search)
    /// when we are only going to report a few PSMs. But we also want to calculate
    /// it for enough candidates that we don't accidentally miss the best hit!
    ///
    /// Given that hyperscore is dominated by the number of matched peaks, it seems
    /// reasonable to assume that the highest hyperscore will belong to one of the
    /// top 50 candidates sorted by # of matched peaks.
    fn trim_hits(&self, hits: &mut InitialHits) {
        let k = 50.clamp(
            (self.report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);
    }

    /// Preliminary Score, return # of matched peaks per candidate
    /// Returned hits are guaranteed to be the top-K hits (see above comment)
    /// from among all potential candidates, but the returned vector is not
    /// in sorted order.
    fn matched_peaks_with_isotope(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> InitialHits {
        let candidates = self.db.query(
            precursor_mass - isotope_error as f32 * NEUTRON,
            precursor_tol,
            self.fragment_tol,
        );

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut hits = InitialHits {
            matched_peaks: 0,
            scored_candidates: 0,
            preliminary: vec![PreScore::default(); potential],
        };

        for peak_mass in query.masses.iter() {
            for charge in 1..max_fragment_charge {
                let mass = peak_mass * charge as f32;
                for frag in candidates.page_search(mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let sc = &mut hits.preliminary[idx];
                    if sc.matched == 0 {
                        hits.scored_candidates += 1;
                        sc.precursor_charge = precursor_charge;
                        sc.peptide = frag.peptide_index;
                        sc.isotope_error = isotope_error;
                    }

                    sc.matched += 1;
                    hits.matched_peaks += 1;
                }
            }
        }
        if hits.matched_peaks == 0 {
            // no candidate matched: don't hand the dense (all-empty) array to the caller
            hits.preliminary = Vec::new();
            return hits;
        }

        self.trim_hits(&mut hits);
        hits
    }

    /// Preliminary scoring for all isotope errors at once.
    ///
    /// Gives the same result as calling [`Self::matched_peaks_with_isotope`] for each
    /// isotope error and concatenating, but binary-searches each fragment page once for
    /// the union of the isotope windows instead of once per window. Within a page,
    /// fragments are sorted by peptide index (i.e. precursor mass) and the windows are
    /// ~1 Da apart, so their sub-ranges are a few entries from each other: one search
    /// plus a short scan replaces one cache-missing binary search per window. This is
    /// the dominant cost of preliminary scoring.
    fn matched_peaks_isotope_windows(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        struct Window {
            isotope: i8,
            idx_lo: usize,
            idx_hi: usize,
            mass_lo: f32,
            mass_hi: f32,
            hits: InitialHits,
        }

        let db = self.db;
        let mut windows = (self.min_isotope_err..=self.max_isotope_err)
            .map(|isotope| {
                let mass = precursor_mass - isotope as f32 * NEUTRON;
                let query = db.query(mass, precursor_tol, self.fragment_tol);
                let (mass_lo, mass_hi) = precursor_tol.bounds(mass);
                Window {
                    isotope,
                    idx_lo: query.pre_idx_lo,
                    idx_hi: query.pre_idx_hi,
                    mass_lo,
                    mass_hi,
                    hits: InitialHits::default(),
                }
            })
            .collect::<Vec<_>>();
        // All windows' dense candidate arrays are alive at once here, versus one at a time
        // in the per-window path. For very wide (Da) precursor windows that can be GBs
        // per thread, so fall back to one window at a time above ~4 M candidates.
        // Decided before allocating anything.
        const MAX_FUSED_CANDIDATES: usize = 1 << 22;
        if windows
            .iter()
            .map(|w| w.idx_hi - w.idx_lo + 1)
            .sum::<usize>()
            > MAX_FUSED_CANDIDATES
        {
            drop(windows);
            return (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        query,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                    );
                    hits
                },
            );
        }
        for w in windows.iter_mut() {
            w.hits.preliminary = vec![PreScore::default(); w.idx_hi - w.idx_lo + 1];
        }
        let union_lo = windows.iter().map(|w| w.idx_lo).min().unwrap_or_default();
        let union_hi = windows.iter().map(|w| w.idx_hi).max().unwrap_or_default();

        // Each (peak, fragment charge, page) lookup costs a few dependent cache misses
        // (skip table, then fragments). Run the lookups as a pipeline, `AHEAD` apart:
        // prefetch a page's skip entries, then narrow by them and prefetch the fragment
        // blocks, then search and scan -- so the misses of several lookups overlap.
        const AHEAD: usize = 8;
        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        let mut lookups: Vec<(usize, f32, f32)> = Vec::new();
        for peak_mass in query.masses.iter() {
            for charge in 1..max_fragment_charge {
                let (fragment_lo, fragment_hi) =
                    self.fragment_tol.bounds(peak_mass * charge as f32);
                let (page_lo, page_hi) = binary_search_slice(
                    &db.min_value,
                    |min, bounds| min.total_cmp(bounds),
                    fragment_lo,
                    fragment_hi,
                );
                lookups.extend((page_lo..page_hi).map(|page| (page, fragment_lo, fragment_hi)));
            }
        }
        let mut blocks = vec![(0usize, 0usize); lookups.len()];
        for i in 0..lookups.len() + 2 * AHEAD {
            if let Some(&(page, _, _)) = lookups.get(i) {
                prefetch(db.page_skip(page));
            }
            if let Some(j) = i.checked_sub(AHEAD).filter(|&j| j < lookups.len()) {
                let (s, e) = db.page_blocks(lookups[j].0, union_lo, union_hi);
                prefetch(&db.fragments[s..e]);
                blocks[j] = (s, e);
            }
            let Some(k) = i.checked_sub(2 * AHEAD).filter(|&k| k < lookups.len()) else {
                continue;
            };
            let (_, fragment_lo, fragment_hi) = lookups[k];
            let (s, e) = blocks[k];
            let (l, r) = binary_search_slice(
                &db.fragments[s..e],
                |frag, bound| (frag.peptide_index.0 as usize).cmp(bound),
                union_lo,
                union_hi,
            );
            for frag in &db.fragments[s + l..s + r] {
                // same (positive) predicate as `page_search`, so NaN never matches
                if !(frag.fragment_mz >= fragment_lo && frag.fragment_mz <= fragment_hi) {
                    continue;
                }
                let ix = frag.peptide_index.0 as usize;
                // windows may overlap (wide Da tolerances): check each of them
                for w in windows.iter_mut() {
                    // same edge handling as `IndexedQuery::page_search`
                    let inside = (ix > w.idx_lo
                        || (ix == w.idx_lo && db.peptides[ix].monoisotopic >= w.mass_lo))
                        && (ix < w.idx_hi
                            || (ix == w.idx_hi && db.peptides[ix].monoisotopic <= w.mass_hi));
                    if !inside {
                        continue;
                    }
                    let sc = &mut w.hits.preliminary[ix - w.idx_lo];
                    if sc.matched == 0 {
                        w.hits.scored_candidates += 1;
                        sc.precursor_charge = precursor_charge;
                        sc.peptide = frag.peptide_index;
                        sc.isotope_error = w.isotope;
                    }
                    sc.matched += 1;
                    w.hits.matched_peaks += 1;
                }
            }
        }

        windows
            .into_iter()
            .fold(InitialHits::default(), |mut acc, mut w| {
                if w.hits.matched_peaks > 0 {
                    self.trim_hits(&mut w.hits);
                } else {
                    w.hits.preliminary = Vec::new();
                }
                acc += w.hits;
                acc
            })
    }

    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        // also for a single isotope window (the default): its lookups are pipelined too
        let mut hits = self.matched_peaks_isotope_windows(
            query,
            precursor_mass,
            precursor_charge,
            precursor_tol,
        );
        if self.min_isotope_err != self.max_isotope_err {
            // one window is trimmed already; trimming again could reorder ties
            self.trim_hits(&mut hits);
        }
        hits
    }

    /// Charges listed for `precursor` when a spectrum carries several precursors at the
    /// same m/z that differ only in charge (the MGF reader's `CHARGE=2+ and 3+`)
    fn listed_charges(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> Option<Vec<u8>> {
        if self.override_precursor_charge || query.precursors.len() < 2 {
            return None;
        }
        let mut charges = query
            .precursors
            .iter()
            .filter(|p| p.mz == precursor.mz)
            .filter_map(|p| p.charge)
            .collect::<Vec<_>>();
        charges.sort_unstable();
        charges.dedup();
        (charges.len() > 1).then_some(charges)
    }

    fn initial_hits(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        // Search in wide-window/DIA mode
        let mut hits = if self.wide_window {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits +=
                        self.matched_peaks(query, precursor_mass, precursor_charge, precursor_tol);
                    hits
                },
            )
        } else if let Some(charges) = self.listed_charges(query, precursor) {
            // Several candidate charges were listed for this precursor (MGF "2+ and 3+"):
            // search exactly those. Only the first precursor used to be searched.
            charges
                .into_iter()
                .fold(InitialHits::default(), |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        self.precursor_tol,
                    );
                    hits
                })
        } else if precursor.charge.is_some() && !self.override_precursor_charge {
            let charge = precursor.charge.unwrap();
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            self.matched_peaks(query, precursor_mass, charge, self.precursor_tol)
        } else {
            // Not all selected ion precursors have charge states annotated (or user has set
            // `override_precursor_charge`)
            // assume it could be z=2, z=3, z=4 and search all three
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        self.precursor_tol,
                    );
                    hits
                },
            )
        };
        self.trim_hits(&mut hits);
        hits
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let hits = self.initial_hits(query, precursor);
        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, precursor, &hits, self.report_psms, &mut features);
        features
    }

    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    fn build_features(
        &self,
        query: &ProcessedSpectrum,
        precursor: &Precursor,
        hits: &InitialHits,
        report_psms: usize,
        features: &mut Vec<Feature>,
    ) {
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .map(|pre| self.score_candidate(query, pre))
            .filter(|s| (s.0.matched_b + s.0.matched_y) >= self.min_matched_peaks)
            .collect::<Vec<_>>();

        // Hyperscore is our primary score function for PSMs
        score_vector.sort_by(|a, b| b.0.hyperscore.total_cmp(&a.0.hyperscore));

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].0;
            let fragments: Option<Fragments> = score_vector[idx].1.take();
            let psm_id = increment_psm_counter();

            let peptide = &self.db[score.peptide];
            let precursor_mass = mz * score.precursor_charge as f32;

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.0.hyperscore)
                .unwrap_or_default();

            let best = score_vector
                .first()
                .map(|score| score.0.hyperscore)
                .expect("we know that index 0 is valid");

            // Poisson distribution log10 probability mass function
            // Computed directly in log space to avoid overflow from lambda.powi(k)
            // log10(PMF) = (k*ln(lambda) - lambda - lnfact(k)) / ln(10)
            let k = score.matched_b + score.matched_y;
            let log10_poisson =
                (k as f64 * lambda.ln() - lambda - lnfact(k)) / std::f64::consts::LN_10;

            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error) * 2E6
                / (precursor_mass - isotope_error + peptide.monoisotopic);

            // let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            features.push(Feature {
                // Identifiers
                psm_id,
                peptide_idx: score.peptide,
                spec_id: query.id.clone(),
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: query
                    .precursors
                    .first()
                    .unwrap()
                    .inverse_ion_mobility
                    .unwrap_or(0.0),
                delta_mass,
                isotope_error,
                average_ppm: score.ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.hyperscore - next,
                delta_best: best - score.hyperscore,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                poisson: if log10_poisson.is_finite() {
                    log10_poisson
                } else {
                    f64::NEG_INFINITY
                },
                longest_b: score.longest_b as u32,
                longest_y: score.longest_y as u32,
                longest_y_pct: score.longest_y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                predicted_ims: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                delta_ims_model: 0.999,
                ms2_intensity: score.summed_b + score.summed_y,

                //Fragments
                protein_groups: None,
                num_protein_groups: 0,
                fragments,
                protein_group_q: 1.0,
            })
        }
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum, psm: &Feature) {
        let peptide = &self.db[psm.peptide_idx];
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind));

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, psm.charge);

        // Remove MS2 peaks matched by previous match
        let mut to_remove = Vec::new();
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                if let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                    &query.masses,
                    &query.intensities,
                    frag.monoisotopic_mass / charge as f32,
                    self.fragment_tol,
                    None,
                ) {
                    to_remove.push((
                        query.masses[peak_idx],
                        query.intensities[peak_idx],
                        query.charges[peak_idx],
                    ));
                }
            }
        }

        let mut masses = Vec::with_capacity(query.masses.len());
        let mut intensities = Vec::with_capacity(query.intensities.len());
        let mut charges = Vec::with_capacity(query.charges.len());
        let mut mobilities = Vec::with_capacity(query.mobilities.len());

        for idx in 0..query.masses.len() {
            let peak = (
                query.masses[idx],
                query.intensities[idx],
                query.charges[idx],
            );
            if !to_remove.contains(&peak) {
                masses.push(query.masses[idx]);
                intensities.push(query.intensities[idx]);
                charges.push(query.charges[idx]);
                if !query.mobilities.is_empty() {
                    mobilities.push(query.mobilities[idx]);
                }
            }
        }

        query.masses = masses;
        query.intensities = intensities;
        query.charges = charges;
        query.mobilities = mobilities;
        query.total_ion_current = query.intensities.iter().sum::<f32>();
    }

    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let mut query = query.clone();
        let hits = self.initial_hits(&query, precursor);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, precursor, &hits, 1, &mut candidates);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Calculate full hyperscore for a given PSM
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
    ) -> (Score, Option<Fragments>) {
        let mut score = Score {
            peptide: pre_score.peptide,
            precursor_charge: pre_score.precursor_charge,
            isotope_error: pre_score.isotope_error,
            ..Default::default()
        };
        let peptide = &self.db[score.peptide];
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate());

        let mut b_run = Run::default();
        let mut y_run = Run::default();

        let mut fragments_details = Fragments::default();

        for (idx, frag) in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                let mz = frag.monoisotopic_mass / charge as f32;

                if let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                    &query.masses,
                    &query.intensities,
                    mz,
                    self.fragment_tol,
                    None,
                ) {
                    let peak_mass = query.masses[peak_idx];
                    let peak_intensity = query.intensities[peak_idx];
                    let fragment_charge = query.charges[peak_idx].max(charge);

                    score.ppm_difference +=
                        peak_intensity * (mz - peak_mass).abs() * 2E6 / (mz + peak_mass);

                    let exp_mz = query.peak_mz(peak_idx);
                    let calc_mz = frag.monoisotopic_mass / fragment_charge as f32 + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak_intensity;
                            b_run.matched(idx);
                        }
                        Kind::X | Kind::Y | Kind::Z => {
                            score.matched_y += 1;
                            score.summed_y += peak_intensity;
                            y_run.matched(idx);
                        }
                    }

                    if self.annotate_matches {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => idx as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z => {
                                peptide.sequence.len().saturating_sub(1) as i32 - idx as i32
                            }
                        };
                        fragments_details.kinds.push(frag.kind);
                        fragments_details.charges.push(fragment_charge as i32);
                        fragments_details.mz_experimental.push(exp_mz);
                        fragments_details.mz_calculated.push(calc_mz);
                        fragments_details.fragment_ordinals.push(idx);
                        fragments_details.intensities.push(peak_intensity);
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(self.score_type);
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;

        if self.annotate_matches {
            (score, Some(fragments_details))
        } else {
            // drop(fragments_details);
            (score, None)
        }
    }
}

/// Maintain information about the longest continous ion ladder for a series
#[derive(Default)]
struct Run {
    start: usize,
    length: usize,
    last: usize,
    pub longest: usize,
}

impl Run {
    pub fn matched(&mut self, index: usize) {
        if self.last == index {
            return;
        } else if self.start + self.length == index {
            self.length += 1;
            self.longest = self.longest.max(self.length);
        } else {
            self.start = index;
            self.length = 1;
            self.longest = self.longest.max(self.length);
        }
        self.last = index;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::database::{Builder, EnzymeBuilder};
    use crate::spectrum::{RawSpectrum, Representation, SpectrumProcessor};

    /// The old preliminary search: one fragment-index search per isotope window
    fn per_window_reference(
        scorer: &Scorer,
        query: &ProcessedSpectrum,
        mass: f32,
        charge: u8,
        tol: Tolerance,
    ) -> InitialHits {
        let mut hits = (scorer.min_isotope_err..=scorer.max_isotope_err).fold(
            InitialHits::default(),
            |mut hits, isotope| {
                hits += scorer.matched_peaks_with_isotope(query, mass, charge, tol, isotope);
                hits
            },
        );
        scorer.trim_hits(&mut hits);
        hits
    }

    #[test]
    fn isotope_windows_match_per_window_search() {
        let fasta = crate::fasta::Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        );
        let mut builder = Builder {
            // many small pages, so that page boundaries are exercised
            bucket_size: Some(64),
            enzyme: Some(EnzymeBuilder {
                missed_cleavages: Some(2),
                min_len: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        };
        builder.update_fasta("unused".into());
        let db = builder.make_parameters().build(fasta);
        assert!(db.min_value.len() > 20, "test needs several pages");

        // Spectrum: singly charged b/y ions of three target peptides, precursor taken
        // from the first one, shifted by +1 isotope so that a non-zero window matches
        let targets = db
            .peptides
            .iter()
            .enumerate()
            .filter(|(_, p)| !p.decoy)
            .map(|(ix, _)| ix as u32)
            .step_by(7)
            .take(3)
            .collect::<Vec<_>>();
        let mut mz = db
            .fragments
            .iter()
            .filter(|f| targets.contains(&f.peptide_index.0))
            .map(|f| f.fragment_mz + PROTON)
            .collect::<Vec<_>>();
        mz.sort_by(f32::total_cmp);
        let precursor_mass = db.peptides[targets[0] as usize].monoisotopic + NEUTRON;
        let raw = RawSpectrum {
            ms_level: 2,
            representation: Representation::Centroid,
            precursors: vec![Precursor {
                mz: precursor_mass / 2.0 + PROTON,
                charge: Some(2),
                ..Default::default()
            }],
            intensity: vec![100.0; mz.len()],
            mz,
            ..Default::default()
        };
        let mut query = SpectrumProcessor::new(150, false, 0.0).process(raw);
        // a non-finite peak must not match anything in either path
        query.masses.push(f32::NAN);
        query.intensities.push(1.0);
        query.charges.push(1);

        for (precursor_tol, fragment_tol) in [
            (Tolerance::Ppm(-10.0, 10.0), Tolerance::Ppm(-20.0, 20.0)),
            // windows 1 Da apart with +/-1.5 Da tolerance overlap
            (Tolerance::Da(-1.5, 1.5), Tolerance::Ppm(-20.0, 20.0)),
            (Tolerance::Ppm(-50.0, 20.0), Tolerance::Da(-0.02, 0.02)),
        ] {
            for (min_iso, max_iso, charge, frag_charge) in [
                (-1, 3, 2, Some(1)),
                (0, 2, 2, Some(1)),
                (0, 1, 2, Some(1)),
                (1, 1, 2, Some(1)),
                // fragment-charge folding with higher precursor charges
                (0, 2, 3, Some(2)),
                (-1, 3, 4, None),
                (0, 1, 4, Some(3)),
            ] {
                let scorer = Scorer {
                    db: &db,
                    precursor_tol,
                    fragment_tol,
                    min_matched_peaks: 2,
                    min_isotope_err: min_iso,
                    max_isotope_err: max_iso,
                    min_precursor_charge: 2,
                    max_precursor_charge: 3,
                    override_precursor_charge: false,
                    max_fragment_charge: frag_charge,
                    chimera: false,
                    report_psms: 1,
                    wide_window: false,
                    annotate_matches: false,
                    score_type: ScoreType::SageHyperScore,
                };
                let expected =
                    per_window_reference(&scorer, &query, precursor_mass, charge, precursor_tol);
                let mut fused = scorer.matched_peaks_isotope_windows(
                    &query,
                    precursor_mass,
                    charge,
                    precursor_tol,
                );
                scorer.trim_hits(&mut fused);
                assert!(
                    expected.matched_peaks > 0,
                    "test spectrum must match something"
                );
                assert_eq!(fused.matched_peaks, expected.matched_peaks);
                assert_eq!(fused.scored_candidates, expected.scored_candidates);
                assert_eq!(fused.preliminary, expected.preliminary);

                // what the scorer actually runs (single window: `isotope_errors [n, n]`)
                let mut actual =
                    scorer.matched_peaks(&query, precursor_mass, charge, precursor_tol);
                if min_iso == max_iso {
                    // same order too: ties are ranked by preliminary order later
                    assert_eq!(actual.preliminary, expected.preliminary);
                }
                let mut expected = expected.preliminary;
                actual.preliminary.sort();
                expected.sort();
                assert_eq!(
                    actual.preliminary, expected,
                    "isotopes {min_iso}..={max_iso}"
                );

                // no fragment matches: no dense all-empty candidate array is returned
                let mut empty = query.clone();
                empty.masses.clear();
                let none = scorer.matched_peaks(&empty, precursor_mass, charge, precursor_tol);
                assert_eq!(none.matched_peaks, 0);
                assert!(none.preliminary.is_empty());
            }
        }
    }

    #[test]
    fn all_listed_precursor_charges_are_searched() {
        let fasta = crate::fasta::Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        );
        let mut builder = Builder::default();
        builder.update_fasta("unused".into());
        let db = builder.make_parameters().build(fasta);
        let target = db.peptides.iter().position(|p| !p.decoy).unwrap() as u32;
        let mut mz = db
            .fragments
            .iter()
            .filter(|f| f.peptide_index.0 == target)
            .map(|f| f.fragment_mz + PROTON)
            .collect::<Vec<_>>();
        mz.sort_by(f32::total_cmp);
        // the peptide is 3+, the spectrum lists "2+ and 3+"
        let precursor_mz = db.peptides[target as usize].monoisotopic / 3.0 + PROTON;
        let precursor = |charge| Precursor {
            mz: precursor_mz,
            charge: Some(charge),
            ..Default::default()
        };
        let raw = RawSpectrum {
            ms_level: 2,
            representation: Representation::Centroid,
            precursors: vec![precursor(2), precursor(3)],
            intensity: vec![100.0; mz.len()],
            mz,
            ..Default::default()
        };
        let query = SpectrumProcessor::new(150, false, 0.0).process(raw);
        let scorer = Scorer {
            db: &db,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 2,
            min_isotope_err: 0,
            max_isotope_err: 0,
            // the configured range does not contain 3: only the listed charges count
            min_precursor_charge: 2,
            max_precursor_charge: 2,
            override_precursor_charge: false,
            max_fragment_charge: Some(1),
            chimera: false,
            report_psms: 1,
            wide_window: false,
            annotate_matches: false,
            score_type: ScoreType::SageHyperScore,
        };
        let psms = scorer.score(&query);
        assert_eq!(psms.len(), 1);
        assert_eq!(psms[0].peptide_idx.0, target);
        assert_eq!(psms[0].charge, 3);
    }

    #[test]
    fn score_heap_keeps_best_hyperscores() {
        let score = |peptide, hyperscore| Score {
            peptide: PeptideIx(peptide),
            hyperscore,
            ..Default::default()
        };
        // best hyperscore at the lowest peptide index
        let mut scores = vec![
            score(0, 50.0),
            score(1, 10.0),
            score(2, 5.0),
            score(3, 20.0),
        ];
        bounded_min_heapify(&mut scores, 2);
        let mut kept = scores[..2].iter().map(|s| s.peptide.0).collect::<Vec<_>>();
        kept.sort();
        assert_eq!(kept, vec![0, 3]);
    }

    #[test]
    fn longest_series() {
        let mut run = Run::default();

        run.matched(1);
        run.matched(2);
        run.matched(3);
        run.matched(3);
        run.matched(3);

        assert_eq!(run.length, 3);
        assert_eq!(run.longest, 3);

        run.matched(5);
        run.matched(5);
        assert_eq!(run.length, 1);
        assert_eq!(run.longest, 3);
        run.matched(6);
        assert_eq!(run.length, 2);
    }

    #[test]
    fn test_max_fragment_charge() {
        assert_eq!(max_fragment_charge(None, 1), 2);
        assert_eq!(max_fragment_charge(None, 2), 2);
        assert_eq!(max_fragment_charge(None, 3), 3);
        assert_eq!(max_fragment_charge(None, 4), 4);
        assert_eq!(max_fragment_charge(Some(1), 2), 2);
        assert_eq!(max_fragment_charge(Some(1), 3), 2);
        assert_eq!(max_fragment_charge(Some(2), 4), 3);
        assert_eq!(max_fragment_charge(Some(4), 1), 2);
    }
}
