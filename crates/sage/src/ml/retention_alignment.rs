//! Perform global retention time alignment using a modified algorithm based on
//! Chen AT, et al. DART-ID https://pubmed.ncbi.nlm.nih.gov/31260443/ (Fig 1A/B)
//!
//! 0) Transform all RTs into unit-less percentages (0.0 - 1.0)
//! 1) Assume that the expected RT for a peptide can be estimated from the average
//!    RT across all runs
//! 2) For each run, calculate a linear regression between the observed peptide RTs
//!    and the global average. Transform all PSM retention times by the regression
//!    parameters
//!
//! If LFQ is enabled, MS1 apex times will be used instead of PSM retention times

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::sync::atomic::AtomicU32;

use super::matrix::Matrix;
use crate::database::PeptideIx;
use crate::scoring::Feature;
use dashmap::DashMap;
use fnv::FnvHasher;
use rayon::prelude::*;

type FnvDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FnvHasher>>;

fn max_rt_by_file(features: &[Feature], n_files: usize) -> Vec<f64> {
    let max_rt = (0..n_files)
        .map(|_| AtomicU32::new(0))
        .map(|_| AtomicU32::new(0))
        .collect::<Vec<_>>();

    features.par_iter().for_each(|feat| {
        max_rt[feat.file_id].fetch_max(feat.rt.ceil() as u32, std::sync::atomic::Ordering::SeqCst);
    });

    max_rt
        .into_iter()
        .map(|v| v.load(std::sync::atomic::Ordering::Acquire) as f64)
        .collect()
}

/// Return a map from PeptideIx to a map from File ID to average RT of the parent
/// PeptideIX
fn mean_rt_by_file(features: &[Feature]) -> FnvDashMap<PeptideIx, HashMap<usize, f64>> {
    let rts: FnvDashMap<PeptideIx, HashMap<usize, f64>> = DashMap::default();
    features
        .par_iter()
        .filter(|feat| feat.label == 1 && feat.spectrum_q <= 0.01)
        .for_each(|feat| {
            rts.entry(feat.peptide_idx)
                .or_default()
                .entry(feat.file_id)
                .and_modify(|f| *f = f.min(feat.rt as f64))
                .or_insert(feat.rt as f64);
        });
    rts
}

fn rt_matrix(features: &[Feature], max_rt: &[f64]) -> (HashMap<PeptideIx, f64>, Matrix) {
    let mean_rt = mean_rt_by_file(features);

    let (means, mat): (HashMap<PeptideIx, f64>, Vec<_>) = mean_rt
        .par_iter()
        .map(|entry| {
            let mut v = vec![f64::NAN; max_rt.len()];

            // While we're here, calculate the mean RT across all runs
            let mut sum = 0.0;
            let mut len = 0.0;
            for (&file_id, &rt) in entry.value() {
                // a file without retention times (max 0) is missing, not 0/0 = NaN, which
                // dropped the whole peptide from the alignment of the other files
                if max_rt[file_id] <= 0.0 {
                    continue;
                }
                let rt = rt / max_rt[file_id];
                v[file_id] = rt;
                sum += rt;
                len += 1.0;
            }

            ((*entry.key(), sum / len), v)
        })
        .filter(|((_, mean), _)| mean.is_normal())
        .unzip();
    let n = mat.len();
    let mat: Vec<f64> = mat.into_par_iter().flatten().collect();

    (means, Matrix::new(mat, n, max_rt.len()))
}

#[derive(Copy, Clone)]
pub struct Alignment {
    pub file_id: usize,
    pub max_rt: f32,
    pub slope: f32,
    pub intercept: f32,
}

pub fn global_alignment(features: &mut [Feature], n_files: usize) -> Vec<Alignment> {
    let max_rt = max_rt_by_file(features, n_files);
    let (_, rt) = rt_matrix(features, &max_rt);

    let mean_rts: Vec<f64> = (0..rt.rows)
        .into_par_iter()
        .map(|row| {
            // Don't include NaN values
            let (len, sum) = rt
                .row(row)
                .filter(|rt| rt.is_finite())
                .fold((0, 0.0f64), |(len, sum), x| (len + 1, sum + x));
            sum / len as f64
        })
        .collect();

    // for file_idx in 0..n_files {
    let alignments = (0..n_files)
        .into_par_iter()
        .map(|file_id| {
            // calculate dot product across all ID'ed peptides
            let (len, dot, sum_x, sum_y) = rt
                .col(file_id)
                .zip(mean_rts.iter())
                .filter(|(x, _)| x.is_finite())
                .fold(
                    (0, 0.0f64, 0.0f64, 0.0f64),
                    |(len, dot, sum_x, sum_y), (x, y)| (len + 1, dot + x * y, sum_x + x, sum_y + y),
                );

            let x_mean = sum_x / len as f64;
            let y_mean = sum_y / len as f64;
            let ssxy = dot - len as f64 * x_mean * y_mean;

            let sx2 = rt
                .col(file_id)
                .filter(|rt| rt.is_finite())
                .fold(1E-8f64, |sum, x| sum + (x - x_mean).powi(2));

            let mut slope = ssxy / sx2;
            let mut intercept = y_mean - slope * x_mean;

            if !slope.is_finite() {
                slope = 1.0;
            }

            if !intercept.is_finite() {
                intercept = 0.0;
            }

            log::info!(
                "aligning file #{file}: y = {m:.4}x + {b:.4}",
                file = file_id,
                m = slope,
                b = intercept
            );

            // Without retention times (e.g. MGF without RTINSECONDS) the maximum is 0 and
            // rt / max_rt = NaN, which made the LDA fail and fall back to the heuristic
            // score for *every* file of the run. Use a constant 0 feature instead.
            let file_max_rt = if max_rt[file_id] > 0.0 {
                max_rt[file_id]
            } else {
                log::warn!(
                    "file #{} has no retention times; RT-based features are not informative for it",
                    file_id
                );
                1.0
            };

            Alignment {
                file_id,
                max_rt: file_max_rt as f32,
                slope: slope as f32,
                intercept: intercept as f32,
            }
        })
        .collect::<Vec<Alignment>>();

    log::info!("aligned retention times across {} files", n_files);

    features.par_iter_mut().for_each(|feature| {
        let a = alignments[feature.file_id];

        // Calculate aligned RT
        // - Divide by maximum RT of this run
        // - Multiply by regression parameters
        feature.aligned_rt = (feature.rt / a.max_rt) * a.slope + a.intercept;
    });

    alignments
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn files_without_retention_times_do_not_produce_nan() {
        let mut features = (0..40)
            .map(|i| Feature {
                peptide_idx: PeptideIx(i % 10),
                file_id: (i % 2) as usize,
                // file 1 has no retention times
                rt: if i % 2 == 0 { 10.0 + i as f32 } else { 0.0 },
                label: 1,
                spectrum_q: 0.0,
                ..Default::default()
            })
            .collect::<Vec<_>>();
        global_alignment(&mut features, 2);
        assert!(features.iter().all(|f| f.aligned_rt.is_finite()));
    }

    #[test]
    fn rt_less_file_does_not_change_other_alignments() {
        // file 0 and file 1 have RTs (file 1 shifted), file 2 has none
        let make = |files: usize| {
            (0..60)
                .filter_map(|i| {
                    let file = i % 3;
                    (file < files).then(|| Feature {
                        peptide_idx: PeptideIx((i / 3) as u32),
                        file_id: file,
                        rt: match file {
                            0 => 10.0 + (i / 3) as f32,
                            1 => 12.0 + 1.1 * (i / 3) as f32,
                            _ => 0.0,
                        },
                        label: 1,
                        spectrum_q: 0.0,
                        ..Default::default()
                    })
                })
                .collect::<Vec<_>>()
        };
        let (mut two, mut three) = (make(2), make(3));
        let a2 = global_alignment(&mut two, 2);
        let a3 = global_alignment(&mut three, 3);
        for f in 0..2 {
            assert!((a2[f].slope - a3[f].slope).abs() < 1e-6, "file {f}");
            assert!((a2[f].intercept - a3[f].intercept).abs() < 1e-6, "file {f}");
        }
    }
}
