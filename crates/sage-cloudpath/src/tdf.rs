use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, path::Path};
use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};
use timsrust::readers::SpectrumReader;
use timsrust::readers::SpectrumReaderConfig as TimsrustSpectrumConfig;
pub struct TdfReader;

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerMS1CentoidingConfig {
    pub mz_ppm: f32,
    pub ims_pct: f32,
}

impl Default for BrukerMS1CentoidingConfig {
    fn default() -> Self {
        BrukerMS1CentoidingConfig {
            mz_ppm: 5.0,
            ims_pct: 3.0,
        }
    }
}

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerProcessingConfig {
    pub ms2: TimsrustSpectrumConfig,
    pub ms1: BrukerMS1CentoidingConfig,
}

impl TdfReader {
    /// Spectra of a .d file, each passed through `f` (see [`crate::util::read_processed`])
    pub fn parse<T: Send>(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerProcessingConfig,
        requires_ms1: bool,
        f: impl Fn(RawSpectrum) -> T + Sync,
    ) -> Result<Vec<T>, timsrust::TimsRustError> {
        let mut spectra = match dda_frames_once(path_name.as_ref(), file_id, config.ms2, &f) {
            Some(spectra) => spectra?,
            None => {
                let spectrum_reader = timsrust::readers::SpectrumReader::build()
                    .with_path(path_name.as_ref())
                    .with_config(config.ms2)
                    .finalize()?;
                self.read_msn_spectra(file_id, &spectrum_reader)?
                    .into_par_iter()
                    .map(&f)
                    .collect()
            }
        };
        if requires_ms1 {
            let ms1s = self.read_ms1_spectra(&path_name, file_id, config.ms1)?;
            spectra.extend(ms1s.into_par_iter().map(&f).collect::<Vec<_>>());
        }

        Ok(spectra)
    }

    fn read_ms1_spectra(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerMS1CentoidingConfig,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let start = std::time::Instant::now();
        let frame_reader = timsrust::readers::FrameReader::new(path_name.as_ref())?;
        let metadata = timsrust::readers::MetadataReader::new(path_name.as_ref())?;
        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;
        let tol_ppm = config.mz_ppm;
        let im_tol_pct = config.ims_pct;

        let ms1_spectra: Vec<RawSpectrum> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS1)
            .map_init(
                || PeakBuffer::with_capacity(2 * MAX_PEAKS),
                |buffer, frame| match frame {
                    Ok(frame) => {
                        buffer.clear();
                        buffer.with_frame(&frame, &ims_converter, &mz_converter);

                        // Squash the mobility dimension
                        let (mz, (intensity, mobility)): (Vec<f32>, (Vec<f32>, Vec<f32>)) =
                            buffer.fastcentroid_frame(tol_ppm, im_tol_pct);

                        let scan_start_time = frame.rt_in_seconds as f32 / 60.0;
                        let ion_injection_time = 100.0; // This is made up, in theory we can read
                                                        // if from the tdf file
                        let total_ion_current = intensity.iter().sum::<f32>();
                        let id = frame.index.to_string();

                        let spec = RawSpectrum {
                            file_id,
                            precursors: vec![],
                            representation: Representation::Centroid,
                            scan_start_time,
                            ion_injection_time,
                            mz,
                            ms_level: 1,
                            id,
                            intensity,
                            total_ion_current,
                            mobility: Some(mobility),
                        };
                        Ok(spec)
                    }
                    // a frame that cannot be read fails the file (it used to be logged and
                    // dropped, so a corrupt .d was searched partially without failing)
                    Err(x) => {
                        log::error!("cannot read an MS1 frame: {}", x);
                        Err(x)
                    }
                },
            )
            .collect::<Result<_, _>>()?;
        log::info!(
            "read {} ms1 spectra in {:#?}",
            ms1_spectra.len(),
            start.elapsed()
        );
        Ok(ms1_spectra)
    }

    fn read_msn_spectra(
        &self,
        file_id: usize,
        spectrum_reader: &SpectrumReader,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let spectra: Vec<RawSpectrum> = (0..spectrum_reader.len())
            .into_par_iter()
            .map(|index| match spectrum_reader.get(index) {
                Ok(dda_spectrum) => Ok(match dda_spectrum.precursor {
                    Some(dda_precursor) => {
                        let mut precursor = Self::parse_precursor(dda_precursor);
                        precursor.isolation_window = Option::from(Tolerance::Da(
                            -dda_spectrum.isolation_width as f32 / 2.0,
                            dda_spectrum.isolation_width as f32 / 2.0,
                        ));
                        let spectrum: RawSpectrum = RawSpectrum {
                            file_id,
                            precursors: vec![precursor],
                            representation: Representation::Centroid,
                            scan_start_time: dda_precursor.rt as f32 / 60.0,
                            // timsrust does not expose the accumulation time; this used to
                            // be filled with the retention time (seconds)
                            ion_injection_time: f32::NAN,
                            total_ion_current: 0.0,
                            mz: dda_spectrum.mz_values.iter().map(|&x| x as f32).collect(),
                            ms_level: 2,
                            id: dda_spectrum.index.to_string(),
                            intensity: dda_spectrum.intensities.iter().map(|&x| x as f32).collect(),
                            mobility: None,
                        };
                        Some(spectrum)
                    }
                    None => None,
                }),
                // a spectrum that cannot be read fails the file instead of being dropped
                Err(e) => {
                    log::error!("cannot read MS2 spectrum {}: {}", index, e);
                    Err(e)
                }
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect();
        Ok(spectra)
    }

    fn parse_precursor(dda_precursor: timsrust::Precursor) -> Precursor {
        let mut precursor: Precursor = Precursor::default();
        precursor.mz = dda_precursor.mz as f32;
        precursor.charge = dda_precursor.charge.map(|x| x as u8);
        precursor.intensity = dda_precursor.intensity.map(|x| x as f32);
        precursor.spectrum_ref = Option::from(dda_precursor.frame_index.to_string());
        precursor.inverse_ion_mobility = Option::from(dda_precursor.im as f32);
        precursor
    }
}

#[derive(Clone, Copy)]
struct ImsPeak {
    mz: f32,
    intensity: f32,
    im: f32,
}
const MAX_PEAKS: usize = 10_000;

/// Buffer that gets re-used on each thread to store the intermediates
/// of the centroiding for a single frame.
#[derive(Clone)]
struct PeakBuffer {
    peaks: Vec<ImsPeak>,
    order: Vec<usize>,
    agg_buff: Vec<ImsPeak>,
}

impl PeakBuffer {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            order: Vec::with_capacity(capacity),
            agg_buff: Vec::with_capacity(MAX_PEAKS),
        }
    }

    fn with_frame(
        &mut self,
        frame: &timsrust::Frame,
        ims_converter: &Scan2ImConverter,
        mz_converter: &Tof2MzConverter,
    ) {
        let expect_len = frame.tof_indices.len();
        self.expand_to_capacity(expect_len);

        let mz_iter = frame
            .tof_indices
            .iter()
            .map(|&x| mz_converter.convert(x as f64) as f32);
        let intensities_iter = frame.intensities.iter().map(|&x| x as f32);
        let imss_iter = Self::expand_mobility_iter(&frame.scan_offsets, ims_converter);

        let peak_iter = mz_iter
            .zip(intensities_iter)
            .zip(imss_iter)
            .map(|((mz, intensity), im)| ImsPeak { mz, intensity, im });
        self.peaks.extend(peak_iter);
        assert_eq!(self.peaks.len(), expect_len);

        // sort by mz ... bc binary searching on the mz space
        // for neighbors is the fastest way to find neighbors that I have tried.
        self.peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

        // The "order" is sorted by intensity
        // This will be used later during the centroiding (for details check that implementation)
        self.order.extend(0..self.len());
        self.order.sort_unstable_by(|&a, &b| {
            self.peaks[b]
                .intensity
                .partial_cmp(&self.peaks[a].intensity)
                .unwrap_or(Ordering::Equal)
        });
    }

    fn clear(&mut self) {
        self.peaks.clear();
        self.order.clear();
        self.agg_buff.clear();
    }

    fn expand_to_capacity(&mut self, capacity: usize) {
        if capacity <= self.len() {
            return;
        }
        let diff = capacity - self.len();
        // Grow by whatever is the largest 20% of the current capacity
        // or the difference.
        let diff = diff.max(self.len() / 5);

        self.peaks.reserve(diff);
        self.order.reserve(diff);
        self.agg_buff.reserve(capacity);
    }

    fn len(&self) -> usize {
        self.peaks.len()
    }

    /// Expand the scan offset slice to mobilities.
    ///
    /// The scan offsets is in essence a run-length
    /// encoded vector of scan numbers that can be converter to the 1/k0
    /// values.
    ///
    /// Essentially ... the slice [0,4,5,5], would expand to
    /// [0,0,0,0,1]; 0 to 4 have index 0, 4 to 5 have index 1, 5 to 5 would
    /// have index 2 but its empty!
    ///
    /// Then this index can be converted using the Scan2ImConverter.convert
    ///
    /// ... This should problably be implemented and exposed in timsrust.
    fn expand_mobility_iter<'a>(
        scan_offsets: &'a [usize],
        ims_converter: &'a Scan2ImConverter,
    ) -> impl Iterator<Item = f32> + 'a {
        let ims_iter = scan_offsets
            .windows(2)
            .enumerate()
            .filter_map(|(i, w)| {
                let num = w[1] - w[0];
                if num == 0 {
                    return None;
                }
                let lo = w[0];
                let hi = w[1];

                let im = ims_converter.convert(i as f64) as f32;
                Some((im, lo, hi))
            })
            .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im));
        ims_iter
    }

    /// Centroiding of the IM-containing spectra
    ///
    /// This is a very rudimentary centroiding algorithm but... it seems to work well.
    /// It iterativelty goes over the peaks in decreasing intensity order and
    /// accumulates the intensity of the peaks surrounding the peak. (sort of
    /// like the first pass in dbscan).
    ///
    /// The preserved mobility and mz are the ones from the apex peak.
    /// A more complex version where the weighted mean is preserved is possible
    /// but I have seen only marginal gains and a lot more complexity + time.
    ///
    /// This dramatically reduces the number of peaks in the spectra
    /// which saves a ton of memory and time when doing LFQ, since we
    /// iterate over each peak.
    fn fastcentroid_frame(
        &mut self,
        mz_tol_ppm: f32,
        im_tol_pct: f32,
    ) -> (Vec<f32>, (Vec<f32>, Vec<f32>)) {
        // Make sure the array is mz sorted ... I should delete
        // this assertions once I am confident of the implementation.
        // but tbh, its not that slow and its simple.
        assert!(
            self.peaks.windows(2).all(|x| x[0].mz <= x[1].mz),
            "mz_array is not sorted"
        );
        assert!(self.agg_buff.is_empty(), "agg_buff is not empty");

        let mut global_num_included = 0;

        let utol = mz_tol_ppm / 1e6;
        let im_tol = im_tol_pct / 100.0;

        for &idx in &self.order {
            if self.peaks[idx].intensity <= 0.0 {
                continue;
            }
            if self.agg_buff.len() > MAX_PEAKS {
                let curr_loc_int = self.peaks[idx].intensity;
                if curr_loc_int > 200.0 {
                    log::debug!(
                        "Reached limit of the agg buffer at index {}/{} curr int={}",
                        idx,
                        self.len(),
                        curr_loc_int
                    );
                }
                break;
            }

            let mz = self.peaks[idx].mz;
            let im = self.peaks[idx].im;
            let da_tol = mz * utol;
            let left_e = mz - da_tol;
            let right_e = mz + da_tol;

            let ss_start = self.peaks.partition_point(|&x| x.mz < left_e);
            let ss_end = self.peaks.partition_point(|&x| x.mz <= right_e);

            let abs_im_tol = im * im_tol;
            let left_im = im - abs_im_tol;
            let right_im = im + abs_im_tol;

            let mut curr_intensity = 0.0;

            let mut num_includable = 0;
            for i in ss_start..ss_end {
                let im_i = self.peaks[i].im;
                if (self.peaks[i].intensity > 0.0) && im_i >= left_im && im_i <= right_im {
                    curr_intensity += self.peaks[i].intensity;
                    self.peaks[i].intensity = -1.0;
                    num_includable += 1;
                }
            }

            assert!(num_includable > 0, "At least 'itself' should be included");

            self.agg_buff.push(ImsPeak {
                mz,
                intensity: curr_intensity,
                im,
            });
            global_num_included += num_includable;

            if global_num_included == self.len() {
                log::debug!("All peaks were included in the centroiding");
                break;
            }
        }

        self.agg_buff
            .sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        // println!("Centroiding: Start len: {}; end len: {};", arr_len, result.len());
        // Ultra data is usually start: 40k end 10k,
        // HT2 data is usually start 400k end 40k, limiting to 10k
        // rarely leaves peaks with intensity > 200 ... ive never seen
        // it happen. -JSP 2025-Jan

        self.agg_buff
            .drain(..)
            .map(|x| (x.mz, (x.intensity, x.im)))
            .unzip()
    }
}

/// One row of the PASEF table: a precursor's scan range in one MS2 frame
struct PasefEntry {
    frame: usize,
    scan_start: usize,
    scan_end: usize,
    isolation_width: f64,
    precursor: usize,
}

/// DDA MS2 spectra as timsrust's `SpectrumReader` returns them (same spectra, same
/// order, same peaks), but each frame decompressed once per block of spectra.
///
/// timsrust assembles a spectrum by decompressing every frame its precursor spans, and a
/// PASEF MS2 frame holds ~9 precursors, so every MS2 frame was decompressed ~9 times.
/// Returns `None` where this does not apply (not DDA-PASEF, calibration requested, an
/// unexpected table): the caller then uses timsrust's reader.
fn dda_frames_once<T: Send>(
    path: &Path,
    file_id: usize,
    config: TimsrustSpectrumConfig,
    f: &(impl Fn(RawSpectrum) -> T + Sync),
) -> Option<Result<Vec<T>, timsrust::TimsRustError>> {
    use timsrust::readers::{FrameReader, MetadataReader, PrecursorReader, TimsTofPath};
    if config.spectrum_processing_params.calibrate {
        return None;
    }
    // only TDF: probing a MiniTDF directory's SQL sidecar could abort in timsrust
    let tims_path = TimsTofPath::new(path).ok()?;
    if !matches!(
        tims_path.file_type(),
        timsrust::readers::TimsTofFileType::TDF
    ) {
        return None;
    }
    let frames = FrameReader::new(path).ok()?;
    if frames.get_acquisition() != timsrust::AcquisitionType::DDAPASEF {
        return None;
    }
    let mz_converter = MetadataReader::new(path).ok()?.mz_converter;
    let precursors = PrecursorReader::build()
        .with_path(path)
        .with_config(config.frame_splitting_params)
        .finalize()
        .ok()?;

    let tdf = tims_path.tdf().ok()?;
    let db = rusqlite::Connection::open_with_flags(tdf, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .ok()?;
    // timsrust's exact query: SQLite may return rows in a different order for another
    // projection (a covering index), and the last row of a precursor sets its width
    let mut stmt = db
        .prepare("SELECT Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy, Precursor FROM PasefFrameMsMsInfo")
        .ok()?;
    // missing values read as 0, as in timsrust
    let mut pasef = stmt
        .query_map([], |row| {
            Ok(PasefEntry {
                frame: row.get(0).unwrap_or_default(),
                scan_start: row.get(1).unwrap_or_default(),
                scan_end: row.get(2).unwrap_or_default(),
                isolation_width: row.get(4).unwrap_or_default(),
                precursor: row.get(6).unwrap_or_default(),
            })
        })
        .ok()?
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    // spectrum i = the i-th precursor id; its frames in table order (stable sort)
    pasef.sort_by_key(|p| p.precursor);
    let mut groups = Vec::new();
    let mut start = 0;
    for i in 1..=pasef.len() {
        if i == pasef.len() || pasef[i].precursor != pasef[start].precursor {
            groups.push(start..i);
            start = i;
        }
    }
    if groups.is_empty() || groups.len() != precursors.len() {
        return None;
    }

    let smoothing = config.spectrum_processing_params.smoothing_window;
    let centroiding = config.spectrum_processing_params.centroiding_window;
    // Blocks of up to 4096 spectra referencing at most 1024 distinct frames: bounds the
    // decoded frames alive at a time (~1000 for a typical ddaPASEF block).
    // ponytail: bounded by frame count, not bytes; frames of a DDA run are small
    const MAX_SPECTRA: usize = 4096;
    const MAX_FRAMES: usize = 1024;
    let mut blocks = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut block_start = 0;
    for (ix, group) in groups.iter().enumerate() {
        let new = pasef[group.clone()]
            .iter()
            .filter(|p| !seen.contains(&p.frame))
            .count();
        if ix > block_start && (ix - block_start == MAX_SPECTRA || seen.len() + new > MAX_FRAMES) {
            blocks.push(block_start..ix);
            block_start = ix;
            seen.clear();
        }
        seen.extend(pasef[group.clone()].iter().map(|p| p.frame));
    }
    blocks.push(block_start..groups.len());

    let mut spectra = Vec::with_capacity(groups.len());
    for block_range in blocks {
        let block = &groups[block_range.clone()];
        let mut needed = block
            .iter()
            .flat_map(|g| pasef[g.clone()].iter().map(|p| p.frame - 1))
            .collect::<Vec<_>>();
        needed.sort_unstable();
        needed.dedup();
        let decoded = match needed
            .par_iter()
            .map(|&f| frames.get(f))
            .collect::<Result<Vec<_>, _>>()
        {
            Ok(decoded) => decoded,
            Err(e) => {
                log::error!("cannot read an MS2 frame: {}", e);
                return Some(Err(e.into()));
            }
        };
        let block_spectra = block
            .par_iter()
            .enumerate()
            .map(|(k, group)| {
                let index = block_range.start + k;
                let mut isolation_width = 0.0;
                let mut tof = Vec::new();
                let mut intensity = Vec::new();
                for p in &pasef[group.clone()] {
                    isolation_width = p.isolation_width;
                    let frame = &decoded[needed.binary_search(&(p.frame - 1)).unwrap()];
                    if frame.intensities.is_empty() {
                        continue;
                    }
                    let lo = frame.scan_offsets[p.scan_start];
                    let hi = frame.scan_offsets[p.scan_end];
                    tof.extend_from_slice(&frame.tof_indices[lo..hi]);
                    intensity.extend(frame.intensities[lo..hi].iter().map(|&x| x as u64));
                }
                let (tof, intensity) = group_and_sum(tof, intensity);
                let intensity = smooth(&tof, intensity, smoothing);
                let keep = local_maxima(&tof, &intensity, centroiding);
                // exact-size peak arrays: spectra are kept for the whole run, and a
                // collected `filter` would carry up to 2x spare capacity
                let kept = keep.iter().filter(|&&k| k).count();
                let mut mz = Vec::with_capacity(kept);
                let mut int = Vec::with_capacity(kept);
                for i in (0..tof.len()).filter(|&i| keep[i]) {
                    mz.push(mz_converter.convert(tof[i]) as f32);
                    int.push(intensity[i] as f64 as f32);
                }
                let tims_precursor = precursors.get(index).unwrap();
                let rt = tims_precursor.rt;
                let mut precursor = TdfReader::parse_precursor(tims_precursor);
                precursor.isolation_window = Some(Tolerance::Da(
                    -isolation_width as f32 / 2.0,
                    isolation_width as f32 / 2.0,
                ));
                f(RawSpectrum {
                    file_id,
                    precursors: vec![precursor],
                    representation: Representation::Centroid,
                    scan_start_time: rt as f32 / 60.0,
                    ion_injection_time: f32::NAN,
                    total_ion_current: 0.0,
                    mz,
                    ms_level: 2,
                    id: index.to_string(),
                    intensity: int,
                    mobility: None,
                })
            })
            .collect::<Vec<_>>();
        spectra.extend(block_spectra);
    }
    Some(Ok(spectra))
}

// The three steps below are timsrust 0.4.2's (`vec_utils::group_and_sum`,
// `RawSpectrum::smooth`, `find_sparse_local_maxima_mask`), which are not public.

/// Sum intensities of equal TOF indices; result sorted by TOF index
fn group_and_sum(tof: Vec<u32>, intensity: Vec<u64>) -> (Vec<u32>, Vec<u64>) {
    let mut order = (0..tof.len()).collect::<Vec<_>>();
    order.sort_by_key(|&i| tof[i]);
    let mut out_tof: Vec<u32> = Vec::with_capacity(order.len());
    let mut out_int: Vec<u64> = Vec::with_capacity(order.len());
    for i in order {
        if out_tof.last() == Some(&tof[i]) {
            *out_int.last_mut().unwrap() += intensity[i];
        } else {
            out_tof.push(tof[i]);
            out_int.push(intensity[i]);
        }
    }
    (out_tof, out_int)
}

/// Add the intensity of every neighbour within `window` TOF indices
fn smooth(tof: &[u32], intensity: Vec<u64>, window: u32) -> Vec<u64> {
    let mut smoothed = intensity.clone();
    for i in 0..tof.len() {
        for j in i + 1..tof.len() {
            if tof[j] - tof[i] > window {
                break;
            }
            smoothed[i] += intensity[j];
            smoothed[j] += intensity[i];
        }
    }
    smoothed
}

/// Peaks that are not lower than a neighbour within `window` TOF indices
fn local_maxima(tof: &[u32], intensity: &[u64], window: u32) -> Vec<bool> {
    let mut keep = vec![true; tof.len()];
    for i in 0..tof.len() {
        for j in i + 1..tof.len() {
            if tof[j] - tof[i] > window {
                break;
            }
            if intensity[i] < intensity[j] {
                keep[i] = false;
            } else {
                keep[j] = false;
            }
        }
    }
    keep
}
