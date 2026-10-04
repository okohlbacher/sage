//! Self-trained, cross-fitted fragment-ion likelihood model (optional PIN/TSV features).
//!
//! For every PSM, the model asks how well the observed fragment evidence fits what confident
//! PSMs of the same run look like. Each theoretical b/y ion (at fragment charge 1, and also 2
//! for precursor charge >= 3) gets a *context* (ion series, precursor-charge bucket, fragment
//! charge, position along the peptide, residues at the cleavage site, whether the
//! complementary ion was found) and an *outcome*: not found, or found at one of 7 intensity
//! rank bins x 3 mass error bins. Two tables of outcome probabilities per context are learned
//! from the run's own confident PSMs: a *signal* table (the PSM's peptide) and a *noise* table
//! (the same peptide reversed except its C-terminal residue, matched against the same
//! spectrum). The features are
//!
//! * `ion_llr`: the sum over all ions of log P(outcome | signal) - log P(outcome | noise);
//! * `ion_explained`: the expected-to-be-found share of ions that were found, weighted by the
//!   signal model's probability of finding each ion.
//!
//! Evidence is every deisotoped peak of the spectrum (before the `max_peaks` cut), kept by the
//! [`crate::spectrum::SpectrumProcessor`] as an [`IonEvidence`]. Training PSMs are the rank-1
//! targets at 1% FDR of a label-free score (the HyperScore on intensities normalised to the
//! spectrum's most intense peak, so that it is comparable between spectra), with one
//! target-decoy competition for precursor charge <= 2 and one for >= 3.
//!
//! Cross-fitting: the spectra of a file are split in two folds by the parity of their index
//! (among the file's spectra with PSMs, in file order). Each fold trains its own model, and
//! every PSM is scored by the model of the *other* fold, so no PSM is scored by a model it
//! helped to train. A file whose folds have fewer than `min_training_psms` training PSMs each
//! (or no decoy PSM at all) gets neutral features (0).
//!
//! The model follows ProSE's `FragmentIonLikelihoodModel` ("rich" contexts, back-off smoothing
//! with a pseudo-count of 20). All counting is integer and every PSM is scored sequentially, so
//! the result does not depend on the number of threads.

use crate::database::IndexedDatabase;
use crate::mass::{monoisotopic, Tolerance};
use crate::scoring::Feature;
use crate::spectrum::ProcessedSpectrum;
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};

/// Ion series: b, y
const SER: usize = 2;
/// Precursor charge buckets: <= 2, 3, >= 4
const PB: usize = 3;
/// Fragment charges: 1, 2
const FC: usize = 2;
/// Position bins along the peptide
const POS: usize = 10;
/// Cleavage site: other, N-terminal to proline, C-terminal to D/E
const SITES: usize = 3;
/// Complementary ion found at charge 1: no, yes
const COMP: usize = 2;
/// Intensity rank bins: 1-2, 3-5, 6-10, 11-20, 21-40, 41-80, > 80
const RB: usize = 7;
/// Mass error bins: < 1/4, < 1/2, <= 1 of the tolerance
const EB: usize = 3;
/// Outcomes: RB x EB found, plus not found
const NOUT: usize = RB * EB + 1;
const ABSENT: usize = NOUT - 1;
const NCTX: usize = SER * PB * FC * POS * SITES * COMP;
/// Peaks ranked individually by intensity; the rest share the last rank bin
const RANKED: usize = 80;
/// Pseudo-count of the back-off smoothing
const PSEUDO: f64 = 20.0;
/// Training PSMs: q-value of the label-free target-decoy competition
const TRAIN_FDR: f64 = 0.01;
/// Minimum number of training PSMs in each fold
pub const MIN_TRAINING_PSMS: usize = 100;

/// Fragment evidence of an MS2 spectrum for the ion model: every deisotoped peak (before the
/// `max_peaks` cut) as a neutral mass (the same mass Sage matches ions against), sorted
/// ascending, with its intensity rank bin.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct IonEvidence {
    pub masses: Vec<f32>,
    pub rank_bins: Vec<u8>,
}

fn rank_bin(rank: usize) -> u8 {
    match rank {
        0..=2 => 0,
        3..=5 => 1,
        6..=10 => 2,
        11..=20 => 3,
        21..=40 => 4,
        41..=80 => 5,
        _ => 6,
    }
}

impl IonEvidence {
    /// `masses` and `intensities` of the peaks, in any order. Equal masses keep their input
    /// order; equal intensities rank the lower mass first.
    pub fn new(masses: &[f32], intensities: &[f32]) -> Self {
        debug_assert_eq!(masses.len(), intensities.len());
        let n = masses.len();
        let mut order = (0..n).collect::<Vec<_>>();
        // stable: equal masses keep their input order
        order.sort_by(|&a, &b| masses[a].total_cmp(&masses[b]));
        let sorted_masses = order.iter().map(|&i| masses[i]).collect::<Vec<_>>();
        let sorted_intensities = order.iter().map(|&i| intensities[i]).collect::<Vec<_>>();

        let mut rank_bins = vec![rank_bin(RANKED + 1); n];
        let by_intensity = |a: &usize, b: &usize| {
            sorted_intensities[*b]
                .total_cmp(&sorted_intensities[*a])
                .then(a.cmp(b))
        };
        let mut ranked = (0..n).collect::<Vec<_>>();
        if n > RANKED {
            ranked.select_nth_unstable_by(RANKED - 1, by_intensity);
            ranked.truncate(RANKED);
        }
        ranked.sort_unstable_by(by_intensity);
        for (rank, &k) in ranked.iter().enumerate() {
            rank_bins[k] = rank_bin(rank + 1);
        }
        Self {
            masses: sorted_masses,
            rank_bins,
        }
    }

    pub fn len(&self) -> usize {
        self.masses.len()
    }

    pub fn is_empty(&self) -> bool {
        self.masses.is_empty()
    }

    /// Index of the peak nearest to `x` (the lower one on a tie), if within `width`
    fn nearest(&self, x: f32, width: f64) -> Option<usize> {
        let ms = &self.masses;
        if ms.is_empty() {
            return None;
        }
        let i = ms.partition_point(|m| *m < x);
        let k = if i == 0 {
            0
        } else if i == ms.len() {
            ms.len() - 1
        } else if (ms[i] - x).abs() < (ms[i - 1] - x).abs() {
            i
        } else {
            i - 1
        };
        ((ms[k] - x).abs() as f64 <= width).then_some(k)
    }
}

/// b and y ion masses of a residue sequence, with the arithmetic of [`crate::ion_series::IonSeries`]:
/// `b[i]` is b(i+1) and `y[i]` is y(n-1-i), for i in 0..n-1
fn ladder(
    residues: impl Iterator<Item = (u8, f32)>,
    nterm: f32,
    monoisotopic_mass: f32,
    b: &mut Vec<f32>,
    y: &mut Vec<f32>,
) {
    b.clear();
    y.clear();
    let mut cb = nterm;
    let mut cy = monoisotopic_mass - nterm;
    let mut residues = residues.peekable();
    while let Some((r, m)) = residues.next() {
        if residues.peek().is_none() {
            break;
        }
        cb += monoisotopic(r) + m;
        cy += -(monoisotopic(r) + m);
        b.push(cb);
        y.push(cy);
    }
}

/// Fragment-ion matching of one peptide sequence against one spectrum's evidence
struct Matcher<'a> {
    evidence: &'a IonEvidence,
    tolerance: Tolerance,
}

impl Matcher<'_> {
    /// Half-width of the matching window around `x`
    fn width(&self, x: f32) -> f64 {
        match self.tolerance {
            Tolerance::Ppm(lo, hi) => x as f64 * ((hi - lo) as f64 / 2.0) * 1e-6,
            Tolerance::Pct(lo, hi) => x as f64 * ((hi - lo) as f64 / 2.0) * 1e-2,
            Tolerance::Da(lo, hi) => (hi - lo) as f64 / 2.0,
        }
    }

    fn outcome(&self, x: f32) -> usize {
        let width = self.width(x);
        match self.evidence.nearest(x, width) {
            None => ABSENT,
            Some(k) => {
                let rel = (self.evidence.masses[k] - x).abs() as f64 / width;
                let eb = if rel < 0.25 {
                    0
                } else if rel < 0.5 {
                    1
                } else {
                    2
                };
                self.evidence.rank_bins[k] as usize * EB + eb
            }
        }
    }

    /// (context, outcome) of every ion and fragment charge of `sequence` (ion masses `b`, `y`
    /// from [`ladder`]) at precursor charge `z`, in the order b1, y(n-1), b2, y(n-2), ...
    fn contexts(&self, sequence: &[u8], b: &[f32], y: &[f32], z: u8, out: &mut Vec<(u16, u8)>) {
        out.clear();
        let n = sequence.len();
        if n < 2 {
            return;
        }
        let bucket = match z {
            0..=2 => 0,
            3 => 1,
            _ => 2,
        };
        let max_charge = (z as usize).saturating_sub(1).clamp(1, 2);
        // ion i of the walk: series s (0 b, 1 y), ordinal o
        let ions = (0..n - 1).flat_map(|i| [(0usize, i + 1, b[i]), (1usize, n - 1 - i, y[i])]);
        // found at fragment charge 1, by series and ordinal
        let mut singly = [vec![false; n], vec![false; n]];
        let mut outcomes = Vec::with_capacity(2 * (n - 1) * max_charge);
        for (s, o, mass) in ions.clone() {
            for charge in 1..=max_charge {
                let outcome = self.outcome(mass / charge as f32);
                if outcome != ABSENT && charge == 1 {
                    singly[s][o] = true;
                }
                outcomes.push(outcome);
            }
        }
        let mut outcomes = outcomes.into_iter();
        for (s, o, _) in ions {
            let pos = (POS * o / n).min(POS - 1);
            // residues N- and C-terminal of the cleavage
            let nside = if s == 0 { o - 1 } else { n - o - 1 };
            let (nres, cres) = (sequence[nside], sequence[nside + 1]);
            let site = if cres == b'P' {
                1
            } else if nres == b'D' || nres == b'E' {
                2
            } else {
                0
            };
            let comp = singly[1 - s][n - o] as usize;
            for charge in 1..=max_charge {
                let ctx = ((((s * PB + bucket) * FC + (charge - 1)) * POS + pos) * SITES + site)
                    * COMP
                    + comp;
                let outcome = outcomes.next().expect("one outcome per ion and charge");
                out.push((ctx as u16, outcome as u8));
            }
        }
    }
}

/// Log-probabilities of the outcomes per context, smoothed by back-off: context ->
/// (series, bucket, fragment charge, site, complement) -> (series, fragment charge) -> global
fn smooth(counts: &[u32]) -> Vec<f64> {
    debug_assert_eq!(counts.len(), NCTX * NOUT);
    let l1n = SER * FC;
    let l2n = SER * PB * FC * SITES * COMP;
    let mut global = [0f64; NOUT];
    let mut l1 = vec![0f64; l1n * NOUT];
    let mut l2 = vec![0f64; l2n * NOUT];
    let mut parent = vec![0usize; l2n];
    let level2 = |ctx: usize| {
        let rest = ctx / (COMP * SITES * POS);
        rest * (SITES * COMP) + ctx % (COMP * SITES)
    };
    for ctx in 0..NCTX {
        let rest = ctx / (COMP * SITES * POS);
        let charge = rest % FC;
        let series = rest / (FC * PB);
        let a1 = series * FC + charge;
        let a2 = level2(ctx);
        parent[a2] = a1;
        for o in 0..NOUT {
            let c = counts[ctx * NOUT + o] as f64;
            global[o] += c;
            l1[a1 * NOUT + o] += c;
            l2[a2 * NOUT + o] += c;
        }
    }
    let global_total = global.iter().sum::<f64>();
    let gp = global
        .iter()
        .map(|g| (g + 1.0) / (global_total + NOUT as f64))
        .collect::<Vec<_>>();
    let back_off = |table: &[f64], row: usize, prior: &[f64]| -> Vec<f64> {
        let row = &table[row * NOUT..(row + 1) * NOUT];
        let total = row.iter().sum::<f64>();
        row.iter()
            .zip(prior)
            .map(|(c, p)| (c + PSEUDO * p) / (total + PSEUDO))
            .collect()
    };
    let l1p = (0..l1n).map(|a| back_off(&l1, a, &gp)).collect::<Vec<_>>();
    let l2p = (0..l2n)
        .map(|a| back_off(&l2, a, &l1p[parent[a]]))
        .collect::<Vec<_>>();
    let mut logp = vec![0f64; NCTX * NOUT];
    for ctx in 0..NCTX {
        let row = &counts[ctx * NOUT..(ctx + 1) * NOUT];
        let total = row.iter().map(|&c| c as f64).sum::<f64>();
        let prior = &l2p[level2(ctx)];
        for o in 0..NOUT {
            logp[ctx * NOUT + o] = ((row[o] as f64 + PSEUDO * prior[o]) / (total + PSEUDO)).ln();
        }
    }
    logp
}

/// A trained model: per (context, outcome) log P(signal) - log P(noise), and per context the
/// signal probability that the ion is found
struct Model {
    llr: Vec<f64>,
    present: Vec<f64>,
}

impl Model {
    fn new(signal: &[u32], noise: &[u32]) -> Self {
        let signal = smooth(signal);
        let noise = smooth(noise);
        let llr = signal.iter().zip(&noise).map(|(s, n)| s - n).collect();
        let present = (0..NCTX)
            .map(|ctx| 1.0 - signal[ctx * NOUT + ABSENT].exp())
            .collect();
        Self { llr, present }
    }

    /// (`ion_llr`, `ion_explained`) of the (context, outcome) list of a PSM
    fn score(&self, contexts: &[(u16, u8)]) -> (f64, f64) {
        let mut llr = 0.0;
        let mut total = 0.0;
        let mut found = 0.0;
        for &(ctx, outcome) in contexts {
            let (ctx, outcome) = (ctx as usize, outcome as usize);
            llr += self.llr[ctx * NOUT + outcome];
            let p = self.present[ctx];
            total += p;
            if outcome != ABSENT {
                found += p;
            }
        }
        (llr, if total > 0.0 { found / total } else { 0.0 })
    }
}

/// Reusable buffers for matching one peptide
#[derive(Default)]
struct Scratch {
    b: Vec<f32>,
    y: Vec<f32>,
    reversed: Vec<u8>,
    contexts: Vec<(u16, u8)>,
}

/// The (context, outcome) list of `feature`'s peptide (or its reversed noise version) against
/// `evidence`; false for a noise peptide equal to the peptide
fn psm_contexts(
    db: &IndexedDatabase,
    feature: &Feature,
    evidence: &IonEvidence,
    tolerance: Tolerance,
    noise: bool,
    scratch: &mut Scratch,
) -> bool {
    let peptide = &db[feature.peptide_idx];
    let sequence = &peptide.sequence[..];
    let n = sequence.len();
    let nterm = peptide.nterm.unwrap_or_default();
    let matcher = Matcher {
        evidence,
        tolerance,
    };
    if !noise {
        let residues = sequence
            .iter()
            .copied()
            .zip(peptide.modifications.iter().copied());
        ladder(
            residues,
            nterm,
            peptide.monoisotopic,
            &mut scratch.b,
            &mut scratch.y,
        );
        matcher.contexts(
            sequence,
            &scratch.b,
            &scratch.y,
            feature.charge,
            &mut scratch.contexts,
        );
        return true;
    }
    // noise: the sequence reversed except its C-terminal residue
    if n < 3 {
        return false;
    }
    let order = (0..n - 1).rev().chain(std::iter::once(n - 1));
    scratch.reversed.clear();
    scratch.reversed.extend(order.clone().map(|i| sequence[i]));
    if scratch.reversed[..] == *sequence {
        return false;
    }
    let residues = order.map(|i| (sequence[i], peptide.modifications[i]));
    ladder(
        residues,
        nterm,
        peptide.monoisotopic,
        &mut scratch.b,
        &mut scratch.y,
    );
    let reversed = std::mem::take(&mut scratch.reversed);
    matcher.contexts(
        &reversed,
        &scratch.b,
        &scratch.y,
        feature.charge,
        &mut scratch.contexts,
    );
    scratch.reversed = reversed;
    true
}

/// Count the signal and noise outcomes of the training PSMs
fn train(
    db: &IndexedDatabase,
    features: &[Feature],
    training: &[(usize, &IonEvidence)],
    tolerance: Tolerance,
) -> Model {
    let tables = || (vec![0u32; NCTX * NOUT], vec![0u32; NCTX * NOUT]);
    let (signal, noise) = training
        .par_iter()
        .fold(
            || (tables(), Scratch::default()),
            |((mut signal, mut noise), mut scratch), &(ix, evidence)| {
                let feature = &features[ix];
                psm_contexts(db, feature, evidence, tolerance, false, &mut scratch);
                for &(ctx, outcome) in &scratch.contexts {
                    signal[ctx as usize * NOUT + outcome as usize] += 1;
                }
                if psm_contexts(db, feature, evidence, tolerance, true, &mut scratch) {
                    for &(ctx, outcome) in &scratch.contexts {
                        noise[ctx as usize * NOUT + outcome as usize] += 1;
                    }
                }
                ((signal, noise), scratch)
            },
        )
        .map(|(tables, _)| tables)
        .reduce(tables, |(mut s1, mut n1), (s2, n2)| {
            s1.iter_mut().zip(&s2).for_each(|(a, b)| *a += b);
            n1.iter_mut().zip(&n2).for_each(|(a, b)| *a += b);
            (s1, n1)
        });
    Model::new(&signal, &noise)
}

/// Rank-1 targets at `TRAIN_FDR` of a target-decoy competition on the normalised HyperScore,
/// one competition per precursor charge class (<= 2, >= 3). `psms` are feature indices.
fn select_training(features: &[Feature], psms: &[usize]) -> Vec<usize> {
    let mut selected = Vec::new();
    for high in [false, true] {
        let mut class = psms
            .iter()
            .copied()
            .filter(|&ix| (features[ix].charge >= 3) == high)
            .collect::<Vec<_>>();
        // best first; on equal scores decoys first (conservative); stable otherwise
        class.sort_by(|&a, &b| {
            features[b]
                .normalized_hyperscore
                .total_cmp(&features[a].normalized_hyperscore)
                .then((features[a].label == 1).cmp(&(features[b].label == 1)))
        });
        let (mut targets, mut decoys) = (0usize, 0usize);
        let mut q = class
            .iter()
            .map(|&ix| {
                if features[ix].label == 1 {
                    targets += 1;
                } else {
                    decoys += 1;
                }
                (decoys + 1) as f64 / targets.max(1) as f64
            })
            .collect::<Vec<_>>();
        for i in (0..q.len().saturating_sub(1)).rev() {
            q[i] = q[i].min(q[i + 1]);
        }
        selected.extend(
            class
                .iter()
                .zip(&q)
                .filter(|(&ix, &q)| features[ix].label == 1 && q <= TRAIN_FDR)
                .map(|(&ix, _)| ix),
        );
    }
    selected
}

/// Outcome of the ion model for one file
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IonModelSummary {
    pub file_id: usize,
    /// Training PSMs of the two folds
    pub training: [usize; 2],
    /// Were the models trained (otherwise the file's features are 0)?
    pub trained: bool,
}

/// Compute `ion_llr` and `ion_explained` for all `features`, per file, from the
/// [`IonEvidence`] of their `spectra`. Features without evidence keep 0.
pub fn annotate(
    features: &mut [Feature],
    spectra: &[ProcessedSpectrum],
    db: &IndexedDatabase,
    tolerance: Tolerance,
    min_training_psms: usize,
) -> Vec<IonModelSummary> {
    // spectrum of each feature, found by its file and native id. An id that more than one
    // spectrum of a file carries (e.g. repeated MGF titles) cannot tell which of them a PSM
    // came from: such PSMs are left out (features 0, not used for training), whatever the
    // order of the spectra.
    let mut by_id: HashMap<(usize, &str), Option<usize>> = HashMap::with_capacity(spectra.len());
    let mut ambiguous = 0usize;
    for (ix, spectrum) in spectra.iter().enumerate() {
        match by_id.entry((spectrum.file_id, spectrum.id.as_str())) {
            Entry::Vacant(entry) => {
                entry.insert(spectrum.ion_evidence.is_some().then_some(ix));
            }
            Entry::Occupied(mut entry) => {
                entry.insert(None);
                ambiguous += 1;
            }
        }
    }
    if ambiguous > 0 {
        log::warn!(
            "ion model: {} spectra repeat the native id of another spectrum of their file; \
             the PSMs of those ids get no ion features",
            ambiguous
        );
    }
    let mut files: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    for (ix, feature) in features.iter().enumerate() {
        if let Some(&Some(spectrum)) = by_id.get(&(feature.file_id, feature.spec_id.as_str())) {
            files
                .entry(feature.file_id)
                .or_default()
                .push((ix, spectrum));
        }
    }
    drop(by_id);

    let mut values = vec![(0f32, 0f32); features.len()];
    let mut summaries = Vec::with_capacity(files.len());
    for (file_id, psms) in files {
        // fold of a spectrum: parity of its index among the file's spectra with PSMs
        let mut with_psms = psms.iter().map(|&(_, s)| s).collect::<Vec<_>>();
        with_psms.sort_unstable();
        with_psms.dedup();
        let fold = |spectrum: usize| {
            with_psms
                .binary_search(&spectrum)
                .expect("spectrum of a PSM")
                % 2
        };
        let evidence = |spectrum: usize| {
            spectra[spectrum]
                .ion_evidence
                .as_deref()
                .expect("spectrum with evidence")
        };

        let mut rank1: [Vec<usize>; 2] = Default::default();
        for &(ix, spectrum) in &psms {
            if features[ix].rank == 1 {
                rank1[fold(spectrum)].push(ix);
            }
        }
        let has_decoys = rank1.iter().flatten().any(|&ix| features[ix].label != 1);
        let training = [
            select_training(features, &rank1[0]),
            select_training(features, &rank1[1]),
        ];
        let sizes = [training[0].len(), training[1].len()];
        let trained = has_decoys && sizes.iter().all(|&n| n >= min_training_psms);
        summaries.push(IonModelSummary {
            file_id,
            training: sizes,
            trained,
        });
        if !trained {
            log::warn!(
                "ion model: file {}: {} + {} training PSMs{} (need {} per fold); ion features are 0",
                file_id,
                sizes[0],
                sizes[1],
                if has_decoys { "" } else { ", no decoys" },
                min_training_psms
            );
            continue;
        }
        log::info!(
            "ion model: file {}: {} + {} training PSMs",
            file_id,
            sizes[0],
            sizes[1]
        );

        let spectrum_of = psms.iter().copied().collect::<HashMap<_, _>>();
        let with_evidence = |set: &[usize]| {
            set.iter()
                .map(|&ix| (ix, evidence(spectrum_of[&ix])))
                .collect::<Vec<_>>()
        };
        let (fold0, fold1) = (with_evidence(&training[0]), with_evidence(&training[1]));
        let (model0, model1) = rayon::join(
            || train(db, features, &fold0, tolerance),
            || train(db, features, &fold1, tolerance),
        );
        let models = [model0, model1];

        let scored = psms
            .par_iter()
            .map_init(Scratch::default, |scratch, &(ix, spectrum)| {
                // scored by the model of the other fold
                let model = &models[1 - fold(spectrum)];
                psm_contexts(
                    db,
                    &features[ix],
                    evidence(spectrum),
                    tolerance,
                    false,
                    scratch,
                );
                let (llr, explained) = model.score(&scratch.contexts);
                (ix, llr as f32, explained as f32)
            })
            .collect::<Vec<_>>();
        for (ix, llr, explained) in scored {
            values[ix] = (llr, explained);
        }
    }
    for (feature, (llr, explained)) in features.iter_mut().zip(values) {
        feature.ion_llr = llr;
        feature.ion_explained = explained;
    }
    summaries
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::{Builder, EnzymeBuilder, PeptideIx};
    use crate::ion_series::{IonSeries, Kind};
    use crate::mass::PROTON;
    use crate::spectrum::{Precursor, RawSpectrum, Representation, SpectrumProcessor};

    #[test]
    fn evidence_sorts_by_mass_and_ranks_by_intensity() {
        let masses = [300.0, 100.0, 200.0, 100.0];
        let intensities = [5.0, 1.0, 5.0, 7.0];
        let ev = IonEvidence::new(&masses, &intensities);
        assert_eq!(ev.masses, vec![100.0, 100.0, 200.0, 300.0]);
        // ranks: 7 (1st), then the two 5s by mass (200 2nd, 300 3rd), then 1 (4th)
        assert_eq!(ev.rank_bins, vec![1, 0, 0, 1]);

        // beyond rank 80 everything is bin 6; equal intensities rank by mass
        let masses = (0..200).map(|i| 1000.0 - i as f32).collect::<Vec<_>>();
        let intensities = vec![1.0; 200];
        let ev = IonEvidence::new(&masses, &intensities);
        let expected = (1..=200).map(rank_bin).collect::<Vec<_>>();
        assert_eq!(ev.rank_bins, expected);
        assert_eq!(rank_bin(80), 5);
        assert_eq!(rank_bin(81), 6);
    }

    #[test]
    fn nearest_prefers_lower_on_ties_and_respects_width() {
        let ev = IonEvidence {
            masses: vec![100.0, 102.0, 104.0],
            rank_bins: vec![0, 1, 2],
        };
        assert_eq!(ev.nearest(101.0, 1.0), Some(0));
        assert_eq!(ev.nearest(101.5, 1.0), Some(1));
        assert_eq!(ev.nearest(99.0, 0.5), None);
        assert_eq!(ev.nearest(105.0, 1.0), Some(2));
        assert_eq!(ev.nearest(110.0, 1.0), None);
        assert_eq!(IonEvidence::default().nearest(100.0, 1.0), None);

        let m = Matcher {
            evidence: &ev,
            tolerance: Tolerance::Da(-0.4, 0.4),
        };
        assert_eq!(m.outcome(102.05), EB); // rank bin 1, error < 1/4
        assert_eq!(m.outcome(102.15), EB + 1); // error < 1/2
        assert_eq!(m.outcome(102.3), EB + 2);
        assert_eq!(m.outcome(103.0), ABSENT);
    }

    fn peptide(seq: &str) -> crate::peptide::Peptide {
        crate::peptide::Peptide::try_from(crate::enzyme::Digest {
            sequence: seq.into(),
            ..Default::default()
        })
        .unwrap()
    }

    #[test]
    fn ladder_matches_ion_series() {
        let mut p = peptide("PEPTMIDEK");
        p.nterm = Some(229.16293);
        p.modifications[4] = 15.994915;
        p.monoisotopic += 229.16293 + 15.994915;
        let (mut b, mut y) = (vec![], vec![]);
        let residues = p
            .sequence
            .iter()
            .copied()
            .zip(p.modifications.iter().copied());
        ladder(residues, p.nterm.unwrap(), p.monoisotopic, &mut b, &mut y);
        let sb = IonSeries::new(&p, Kind::B)
            .map(|i| i.monoisotopic_mass)
            .collect::<Vec<_>>();
        let sy = IonSeries::new(&p, Kind::Y)
            .map(|i| i.monoisotopic_mass)
            .collect::<Vec<_>>();
        assert_eq!(b, sb);
        assert_eq!(y, sy);
        assert_eq!(b.len(), 8);
    }

    #[test]
    fn contexts_follow_the_specification() {
        // "ADPK" at z = 3: ions b1 y3 b2 y2 b3 y1, fragment charges 1 and 2
        let seq = b"ADPK";
        let p = peptide("ADPK");
        let (mut b, mut y) = (vec![], vec![]);
        ladder(
            seq.iter().map(|&r| (r, 0.0)),
            0.0,
            p.monoisotopic,
            &mut b,
            &mut y,
        );
        // only b2 (singly) and y2 (singly) are present
        let ev = IonEvidence::new(&[b[1], y[1]], &[10.0, 5.0]);
        let m = Matcher {
            evidence: &ev,
            tolerance: Tolerance::Ppm(-10.0, 10.0),
        };
        let mut out = vec![];
        m.contexts(seq, &b, &y, 3, &mut out);
        assert_eq!(out.len(), 12);
        let ctx = |s: usize, bucket: usize, c: usize, pos: usize, site: usize, comp: usize| {
            (((((s * PB + bucket) * FC + c) * POS + pos) * SITES + site) * COMP + comp) as u16
        };
        // b1: A|D, site other, complement y3 absent
        assert_eq!(out[0], (ctx(0, 1, 0, 2, 0, 0), ABSENT as u8));
        assert_eq!(out[1], (ctx(0, 1, 1, 2, 0, 0), ABSENT as u8));
        // y3: A|D
        assert_eq!(out[2], (ctx(1, 1, 0, 7, 0, 0), ABSENT as u8));
        // b2: D|P -> N-terminal to proline; complement y2 found; rank 1, exact
        assert_eq!(out[4], (ctx(0, 1, 0, 5, 1, 1), 0));
        assert_eq!(out[5].1, ABSENT as u8);
        // y2: D|P, complement b2 found, rank 2 (bin 0)
        assert_eq!(out[6], (ctx(1, 1, 0, 5, 1, 1), 0));
        // b3: P|K, site other (N-side P is not D/E)
        assert_eq!(out[8], (ctx(0, 1, 0, 7, 0, 0), ABSENT as u8));

        // z = 2: fragment charge 1 only
        m.contexts(seq, &b, &y, 2, &mut out);
        assert_eq!(out.len(), 6);
        assert_eq!(out[2], (ctx(0, 0, 0, 5, 1, 1), 0));
    }

    #[test]
    fn smoothing_gives_distributions() {
        let mut counts = vec![0u32; NCTX * NOUT];
        for (i, c) in counts.iter_mut().enumerate() {
            *c = ((i * 7919) % 13) as u32;
        }
        counts[5 * NOUT..6 * NOUT].iter_mut().for_each(|c| *c = 0);
        let logp = smooth(&counts);
        for ctx in 0..NCTX {
            let total = logp[ctx * NOUT..(ctx + 1) * NOUT]
                .iter()
                .map(|l| l.exp())
                .sum::<f64>();
            assert!((total - 1.0).abs() < 1e-9, "context {ctx}: {total}");
        }
        // an unseen context falls back to its parents, which are non-degenerate
        assert!(logp[5 * NOUT..6 * NOUT].iter().all(|l| l.is_finite()));
    }

    fn feature(ix: usize, score: f64, label: i32, charge: u8) -> Feature {
        Feature {
            peptide_idx: PeptideIx(ix as u32),
            normalized_hyperscore: score,
            label,
            charge,
            rank: 1,
            ..Default::default()
        }
    }

    #[test]
    fn training_selection_competes_per_charge_class() {
        // charge 2: 300 targets above 1 decoy, then mixed; charge 3: decoys on top
        let mut features = vec![];
        for i in 0..300 {
            features.push(feature(i, 100.0 - i as f64 * 0.1, 1, 2));
        }
        features.push(feature(300, 69.0, -1, 2));
        features.push(feature(301, 68.0, 1, 2));
        features.push(feature(302, 50.0, -1, 3));
        features.push(feature(303, 49.0, 1, 3));
        // a tie between a target and a decoy: the decoy counts first
        features.push(feature(304, 100.0, -1, 3));
        features.push(feature(305, 100.0, 1, 3));
        let all = (0..features.len()).collect::<Vec<_>>();
        let selected = select_training(&features, &all);
        // charge 2: q of the first 300 targets is 1/300 <= 1%, after the decoy 2/301
        assert_eq!(selected.len(), 301);
        assert!(selected.iter().all(|&ix| features[ix].charge == 2));
    }

    /// A small search: target peptides of the test FASTA, spectra with their b/y ions plus
    /// noise peaks, scored with the real `Scorer`
    fn small_search() -> (IndexedDatabase, Vec<ProcessedSpectrum>, Vec<Feature>) {
        let fasta = crate::fasta::Fasta::parse(
            include_str!("../../../tests/Q99536.fasta").into(),
            "rev_",
            true,
        );
        let mut builder = Builder {
            enzyme: Some(EnzymeBuilder {
                missed_cleavages: Some(2),
                min_len: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        };
        builder.update_fasta("unused".into());
        let db = builder.make_parameters().build(fasta);
        let processor = SpectrumProcessor::new(30, true, 0.0).with_ion_evidence(true);
        let mut spectra = vec![];
        let mut seed = 12345u64;
        let mut rand = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as f32 / (1u64 << 31) as f32
        };
        let peptides = |decoy: bool| {
            db.peptides
                .iter()
                .enumerate()
                .filter(|(_, p)| p.decoy == decoy && p.sequence.len() >= 7)
                .map(|(ix, _)| ix)
                .collect::<Vec<_>>()
        };
        // 6 copies of a spectrum of every target peptide per file (70% of the b/y ions), and
        // one spectrum of every decoy peptide per file (15% of the ions: decoys win these
        // spectra, but with lower scores)
        let (targets, decoys) = (peptides(false), peptides(true));
        let jobs = (0..12)
            .flat_map(|copy| targets.iter().map(move |&ix| (copy, ix, 0.7)))
            .chain((0..2).flat_map(|copy| decoys.iter().map(move |&ix| (copy, ix, 0.15))));
        for (copy, ix, coverage) in jobs {
            {
                let p = &db.peptides[ix];
                let mut mz = vec![];
                let mut intensity = vec![];
                for kind in [Kind::B, Kind::Y] {
                    for ion in IonSeries::new(p, kind) {
                        if rand() < coverage {
                            mz.push(ion.monoisotopic_mass + PROTON);
                            intensity.push(100.0 + 1000.0 * rand());
                        }
                    }
                }
                for _ in 0..60 {
                    mz.push(150.0 + 1500.0 * rand());
                    intensity.push(300.0 * rand());
                }
                let mut order = (0..mz.len()).collect::<Vec<_>>();
                order.sort_by(|&a, &b| mz[a].total_cmp(&mz[b]));
                let raw = RawSpectrum {
                    ms_level: 2,
                    id: format!("scan={}", spectra.len() + 1),
                    file_id: copy % 2,
                    representation: Representation::Centroid,
                    precursors: vec![Precursor {
                        mz: p.monoisotopic / 2.0 + PROTON,
                        charge: Some(2),
                        ..Default::default()
                    }],
                    mz: order.iter().map(|&i| mz[i]).collect(),
                    intensity: order.iter().map(|&i| intensity[i]).collect(),
                    ..Default::default()
                };
                spectra.push(processor.process(raw));
            }
        }
        let scorer = crate::scoring::Scorer {
            db: &db,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 2,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: 2,
            max_precursor_charge: 2,
            override_precursor_charge: false,
            max_fragment_charge: None,
            chimera: false,
            report_psms: 3,
            wide_window: false,
            annotate_matches: false,
            score_type: crate::scoring::ScoreType::SageHyperScore,
        };
        let features = spectra
            .iter()
            .flat_map(|s| scorer.score(s))
            .collect::<Vec<_>>();
        (db, spectra, features)
    }

    #[test]
    fn evidence_keeps_all_deisotoped_peaks() {
        let (_, spectra, _) = small_search();
        for s in &spectra {
            let ev = s.ion_evidence.as_deref().unwrap();
            assert!(ev.len() > s.masses.len(), "top 30 cut vs all peaks");
            assert!(s.masses.iter().all(|m| ev.masses.contains(m)));
        }
        let off = SpectrumProcessor::new(30, true, 0.0).process(RawSpectrum {
            ms_level: 2,
            representation: Representation::Centroid,
            mz: vec![200.0, 300.0],
            intensity: vec![1.0, 2.0],
            ..Default::default()
        });
        assert!(off.ion_evidence.is_none());
    }

    #[test]
    fn annotate_is_cross_fitted_deterministic_and_thread_independent() {
        let (db, spectra, features) = small_search();
        assert!(features.iter().any(|f| f.label == -1));
        let tol = Tolerance::Ppm(-10.0, 10.0);
        let run = |threads: usize, min: usize| {
            let mut f = features.clone();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let summaries = pool.install(|| annotate(&mut f, &spectra, &db, tol, min));
            (f, summaries)
        };
        let (one, summaries) = run(1, 10);
        assert_eq!(summaries.len(), 2);
        assert!(summaries.iter().all(|s| s.trained), "{summaries:?}");
        let (four, summaries4) = run(4, 10);
        assert_eq!(summaries, summaries4);
        for (a, b) in one.iter().zip(&four) {
            assert_eq!(a.ion_llr.to_bits(), b.ion_llr.to_bits());
            assert_eq!(a.ion_explained.to_bits(), b.ion_explained.to_bits());
        }
        // the true peptides score higher than the other candidates
        let mean = |rank: u32| {
            let v = one
                .iter()
                .filter(|f| f.rank == rank)
                .map(|f| f.ion_llr as f64)
                .collect::<Vec<_>>();
            v.iter().sum::<f64>() / v.len() as f64
        };
        assert!(mean(1) > mean(2) + 1.0, "{} vs {}", mean(1), mean(2));
        assert!(one.iter().all(|f| (0.0..=1.0).contains(&f.ion_explained)));

        // too few training PSMs: neutral
        let (none, summaries) = run(2, 100_000);
        assert!(summaries.iter().all(|s| !s.trained));
        assert!(none
            .iter()
            .all(|f| f.ion_llr == 0.0 && f.ion_explained == 0.0));
    }

    #[test]
    fn no_psm_is_scored_by_a_model_it_trained() {
        // Removing training PSMs of fold 0 (lowering their label-free score) must leave the
        // features of every fold-0 PSM unchanged (they are scored by fold 1's model) and
        // change those of fold-1 PSMs (scored by fold 0's model)
        let (db, spectra, features) = small_search();
        let tol = Tolerance::Ppm(-10.0, 10.0);
        let mut base = features.clone();
        annotate(&mut base, &spectra, &db, tol, 10);

        let spectrum_ix = |f: &Feature| {
            spectra
                .iter()
                .position(|s| s.file_id == f.file_id && s.id == f.spec_id)
                .unwrap()
        };
        // fold of a feature: parity among its file's spectra with PSMs, in order
        let fold = |f: &Feature| {
            let mut with_psms = features
                .iter()
                .filter(|g| g.file_id == f.file_id)
                .map(spectrum_ix)
                .collect::<Vec<_>>();
            with_psms.sort_unstable();
            with_psms.dedup();
            with_psms.binary_search(&spectrum_ix(f)).unwrap() % 2
        };
        let folds = features.iter().map(fold).collect::<Vec<_>>();

        let mut changed = features.clone();
        let mut removed = 0;
        for (f, &fo) in changed.iter_mut().zip(&folds) {
            if fo == 0 && f.rank == 1 && f.label == 1 && removed % 10 == 0 {
                f.normalized_hyperscore = -1.0;
            }
            if fo == 0 && f.rank == 1 && f.label == 1 {
                removed += 1;
            }
        }
        let summaries = annotate(&mut changed, &spectra, &db, tol, 10);
        assert!(summaries.iter().all(|s| s.trained), "{summaries:?}");

        let mut fold1_changed = 0;
        for ((a, b), &fo) in base.iter().zip(&changed).zip(&folds) {
            if fo == 0 {
                assert_eq!(a.ion_llr.to_bits(), b.ion_llr.to_bits());
                assert_eq!(a.ion_explained.to_bits(), b.ion_explained.to_bits());
            } else {
                fold1_changed += (a.ion_llr != b.ion_llr) as usize;
            }
        }
        assert!(fold1_changed > 0);
    }

    #[test]
    fn psms_of_repeated_native_ids_get_no_ion_features() {
        let (db, mut spectra, mut features) = small_search();
        let tol = Tolerance::Ppm(-10.0, 10.0);
        // two spectra of file 0 with PSMs get the same native id (e.g. repeated MGF titles)
        let ids = spectra
            .iter()
            .filter(|s| {
                s.file_id == 0 && features.iter().any(|f| f.file_id == 0 && f.spec_id == s.id)
            })
            .map(|s| s.id.clone())
            .take(2)
            .collect::<Vec<_>>();
        let (a, b) = (ids[0].clone(), ids[1].clone());
        for s in spectra.iter_mut().filter(|s| s.file_id == 0 && s.id == b) {
            s.id = a.clone();
        }
        for f in features
            .iter_mut()
            .filter(|f| f.file_id == 0 && f.spec_id == b)
        {
            f.spec_id = a.clone();
        }
        let repeated = |f: &Feature| f.file_id == 0 && f.spec_id == a;
        let mut annotated = features.clone();
        let summaries = annotate(&mut annotated, &spectra, &db, tol, 10);
        assert!(summaries.iter().all(|s| s.trained), "{summaries:?}");
        assert!(annotated.iter().filter(|f| repeated(f)).count() >= 2);
        assert!(annotated
            .iter()
            .filter(|f| repeated(f))
            .all(|f| f.ion_llr == 0.0 && f.ion_explained == 0.0));
        assert!(annotated.iter().any(|f| !repeated(f) && f.ion_llr != 0.0));

        // the same whichever of the two spectra comes first
        let at = spectra
            .iter()
            .enumerate()
            .filter(|(_, s)| s.file_id == 0 && s.id == a)
            .map(|(ix, _)| ix)
            .collect::<Vec<_>>();
        assert_eq!(at.len(), 2);
        spectra.swap(at[0], at[1]);
        let mut swapped = features.clone();
        annotate(&mut swapped, &spectra, &db, tol, 10);
        for (x, y) in annotated.iter().zip(&swapped) {
            assert_eq!(x.ion_llr.to_bits(), y.ion_llr.to_bits());
            assert_eq!(x.ion_explained.to_bits(), y.ion_explained.to_bits());
        }
    }
}
