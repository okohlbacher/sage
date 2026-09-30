//! mzPeak input. Reading is done by the HUPO-PSI reference reader (`mzpeak_prototyping`)
//! into `mzdata`'s spectrum model; this module only maps those spectra to [`RawSpectrum`].

use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerSpectrum, RefPeakDataLevel, SignalContinuity};
use mzpeak_prototyping::MzPeakReader;
use sage_core::mass::Tolerance;
use sage_core::spectrum::{Precursor, RawSpectrum, Representation};

/// All spectra of an mzPeak reader; MS1 only if `keep_ms1` (LFQ)
pub fn read(mut reader: MzPeakReader, file_id: usize, keep_ms1: bool) -> Vec<RawSpectrum> {
    reader
        .iter()
        .filter(|s| keep_ms1 || s.ms_level() != 1)
        .map(|s| to_raw(&s, file_id))
        .collect()
}

fn to_raw(s: &MultiLayerSpectrum, file_id: usize) -> RawSpectrum {
    let mut raw = RawSpectrum::default_with_file_id(file_id);
    raw.ms_level = s.ms_level();
    raw.id = s.id().to_string();
    // minutes, as in mzML
    raw.scan_start_time = s.start_time() as f32;
    raw.ion_injection_time = s
        .acquisition()
        .first_scan()
        .map(|scan| scan.injection_time)
        .unwrap_or(0.0);
    raw.representation = match s.signal_continuity() {
        SignalContinuity::Centroid => Representation::Centroid,
        _ => Representation::Profile,
    };
    for p in s.precursor_iter() {
        let Some(ion) = p.ion() else { continue };
        let w = &p.isolation_window;
        // mzdata stores absolute window bounds, Sage offsets from the target
        let isolation_window = (w.lower_bound > 0.0 && w.upper_bound > 0.0)
            .then(|| Tolerance::Da(-(w.target - w.lower_bound), w.upper_bound - w.target));
        raw.precursors.push(Precursor {
            mz: ion.mz as f32,
            intensity: Some(ion.intensity),
            charge: ion.charge.map(|c| c as u8),
            spectrum_ref: p.precursor_id.clone(),
            isolation_window,
            inverse_ion_mobility: ion.ion_mobility().map(|v| v as f32),
        });
    }
    match s.peaks() {
        RefPeakDataLevel::Centroid(peaks) => {
            raw.mz = peaks.iter().map(|p| p.mz as f32).collect();
            raw.intensity = peaks.iter().map(|p| p.intensity).collect();
        }
        RefPeakDataLevel::RawData(arrays) => {
            raw.mz = arrays
                .mzs()
                .map(|mz| mz.iter().map(|&x| x as f32).collect())
                .unwrap_or_default();
            raw.intensity = arrays
                .intensities()
                .map(|int| int.to_vec())
                .unwrap_or_default();
        }
        // deconvoluted (neutral mass) peaks or no data: nothing Sage can search
        _ => {}
    }
    raw.total_ion_current = raw.intensity.iter().sum();
    raw
}

#[cfg(test)]
mod test {
    use crate::util::{read_mzml, read_spectra};

    /// The same spectra as the mzML the fixtures were converted from (mzPeak reference
    /// converter; point and chunked layouts)
    #[test]
    fn mzpeak_matches_mzml() {
        let mzml = crate::to_url("../../tests/LQSRPAAPPAPGPGQLTLR.mzML").unwrap();
        let expected = read_mzml(&mzml, 0, None).unwrap();
        for fixture in [
            "../../tests/LQSRPAAPPAPGPGQLTLR.mzpeak",
            "../../tests/LQSRPAAPPAPGPGQLTLR.chunked.mzpeak",
        ] {
            let url = crate::to_url(fixture).unwrap();
            let actual = read_spectra(&url, 0, None, Default::default(), true).unwrap();
            assert_eq!(actual.len(), expected.len(), "{fixture}");
            for (a, e) in actual.iter().zip(&expected) {
                assert_eq!(a.id, e.id);
                assert_eq!(a.ms_level, e.ms_level);
                assert!((a.scan_start_time - e.scan_start_time).abs() < 1e-4);
                assert_eq!(a.ion_injection_time, e.ion_injection_time);
                assert_eq!(a.representation, e.representation);
                assert_eq!(a.mz, e.mz);
                assert_eq!(a.intensity, e.intensity);
                assert_eq!(a.precursors.len(), e.precursors.len());
                for (p, q) in a.precursors.iter().zip(&e.precursors) {
                    assert_eq!(p.mz, q.mz);
                    assert_eq!(p.charge, q.charge);
                    assert_eq!(p.intensity, q.intensity);
                    assert_eq!(p.isolation_window, q.isolation_window);
                    assert_eq!(p.spectrum_ref, q.spectrum_ref);
                }
            }
        }
    }
}
