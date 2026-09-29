use async_compression::tokio::bufread::ZlibDecoder;
use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;
use sage_core::spectrum::{Precursor, Representation};
use sage_core::{mass::Tolerance, spectrum::RawSpectrum};
use std::collections::HashMap;
use tokio::io::{AsyncBufRead, AsyncReadExt};

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Which tag are we inside?
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
    Noise,
}

#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

// MUST supply only one of the following
const ZLIB_COMPRESSION: &[u8] = b"MS:1000574";
const NO_COMPRESSION: &[u8] = b"MS:1000576";

// MUST supply only one of the following
const NUMPRESS_LINEAR: &[u8] = b"MS:1002312";
const NUMPRESS_PIC: &[u8] = b"MS:1002313";
const NUMPRESS_SLOF: &[u8] = b"MS:1002314";
const NUMPRESS_LINEAR_ZLIB: &[u8] = b"MS:1002746";
const NUMPRESS_PIC_ZLIB: &[u8] = b"MS:1002747";
const NUMPRESS_SLOF_ZLIB: &[u8] = b"MS:1002748";

const INTENSITY_ARRAY: &[u8] = b"MS:1000515";
const MZ_ARRAY: &[u8] = b"MS:1000514";
const NOISE_ARRAY: &[u8] = b"MS:1002744";

// MUST supply only one of the following
const FLOAT_64: &[u8] = b"MS:1000523";
const FLOAT_32: &[u8] = b"MS:1000521";

const MS_LEVEL: &[u8] = b"MS:1000511";
const PROFILE: &[u8] = b"MS:1000128";
const CENTROID: &[u8] = b"MS:1000127";
const TOTAL_ION_CURRENT: &[u8] = b"MS:1000285";

const SCAN_START_TIME: &[u8] = b"MS:1000016";
const UNIT_SECONDS: &[u8] = b"UO:0000010";
const UNIT_MINUTES: &[u8] = b"UO:0000031";
const ION_INJECTION_TIME: &[u8] = b"MS:1000927";

const SELECTED_ION_MZ: &[u8] = b"MS:1000744";
const SELECTED_ION_INT: &[u8] = b"MS:1000042";
const SELECTED_ION_CHARGE: &[u8] = b"MS:1000041";

const ISO_WINDOW_TARGET: &[u8] = b"MS:1000827";
const ISO_WINDOW_LOWER: &[u8] = b"MS:1000828";
const ISO_WINDOW_UPPER: &[u8] = b"MS:1000829";

const INVERSE_ION_MOBILITY: &[u8] = b"MS:1002815";

pub struct MzMLReader {
    ms_level: Option<u8>,
    // If set to Some(level) and noise intensities are present in the MzML file,
    // divide intensities at this MS-level by noise to calculate S/N
    signal_to_noise: Option<u8>,

    file_id: usize,
}

impl MzMLReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_file_id_and_level_filter(file_id: usize, ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
            file_id,
            signal_to_noise: None,
        }
    }

    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            ms_level: None,
            signal_to_noise: None,
            file_id,
        }
    }

    pub fn set_file_id(&mut self, file_id: usize) -> &mut Self {
        self.file_id = file_id;
        self
    }

    pub fn set_signal_to_noise(&mut self, sn: Option<u8>) -> &mut Self {
        self.signal_to_noise = sn;
        self
    }

    /// Here be dragons -
    /// Seriously, this kinda sucks because it's a giant imperative, stateful loop.
    /// But I also don't want to spend any more time working on an mzML parser...
    pub async fn parse<B: AsyncBufRead + Unpin>(
        &self,
        b: B,
    ) -> Result<Vec<RawSpectrum>, MzMLError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut output_buffer = Vec::with_capacity(4096);
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = None;

        let mut spectrum = RawSpectrum::default_with_file_id(self.file_id);
        let mut precursor = Precursor::default();
        let mut iso_window_lo: Option<f32> = None;
        let mut iso_window_hi: Option<f32> = None;
        let mut spectra = Vec::new();

        let mut noise_array = Vec::new();

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or(MzMLError::Malformed)?
                    .value
            };
        }

        macro_rules! extract_value {
            ($ev:expr) => {{
                let s = $ev
                    .try_get_attribute(b"value")?
                    .ok_or(MzMLError::Malformed)?
                    .value;
                std::str::from_utf8(&s)?.parse()?
            }};
        }

        let mut param_groups: HashMap<Vec<u8>, Vec<BytesStart<'static>>> = HashMap::new();
        let mut current_group: Option<Vec<u8>> = None;

        // cvParam handling, shared by inline cvParams and by cvParams replayed from a
        // referenceableParamGroup (upstream #232: some converters put e.g. the centroid
        // flag or the MS level into a group)
        macro_rules! handle_cv {
            ($ev:expr) => {{
                let ev = $ev;
                match state {
                    Some(State::BinaryDataArray) => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                            MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                            NOISE_ARRAY => binary_array = Some(BinaryKind::Noise),
                            NUMPRESS_LINEAR | NUMPRESS_PIC | NUMPRESS_SLOF
                            | NUMPRESS_LINEAR_ZLIB | NUMPRESS_PIC_ZLIB | NUMPRESS_SLOF_ZLIB => {
                                // Decoding these bytes as floats would give garbage
                                return Err(MzMLError::UnsupportedCV(format!(
                                    "{} (MS-Numpress compression); convert without numpress",
                                    String::from_utf8_lossy(&accession)
                                )));
                            }
                            // Other cvParams (units, array names...) don't change the array
                            // type; an array whose type is never recognised stays `None` and
                            // is skipped
                            _ => {}
                        }
                    }
                    Some(State::Spectrum) => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            MS_LEVEL => {
                                let level = extract_value!(ev);
                                if let Some(filter) = self.ms_level {
                                    if level != filter {
                                        spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                        state = None;
                                    }
                                }
                                spectrum.ms_level = level;
                            }
                            PROFILE => spectrum.representation = Representation::Profile,
                            CENTROID => spectrum.representation = Representation::Centroid,
                            TOTAL_ION_CURRENT => {
                                let value = extract_value!(ev);
                                if value == 0.0 {
                                    // No ion current, break out of current state
                                    spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                    state = None;
                                } else {
                                    spectrum.total_ion_current = value;
                                }
                            }
                            _ => {}
                        }
                    }
                    Some(State::Precursor) => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ISO_WINDOW_TARGET => {
                                // use isolation window target for precursor m/z, e.g. to handle
                                // DIA setups where the mzML conversion software doesn't write
                                // a selection ion tag
                                if precursor.mz == 0.0 {
                                    precursor.mz = extract_value!(ev)
                                }
                            }
                            ISO_WINDOW_LOWER => iso_window_lo = Some(extract_value!(ev)),
                            ISO_WINDOW_UPPER => iso_window_hi = Some(extract_value!(ev)),
                            _ => {}
                        }
                    }
                    Some(State::SelectedIon) => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SELECTED_ION_CHARGE => {
                                precursor.charge = Some(extract_value!(ev));
                            }
                            SELECTED_ION_MZ => {
                                let val = extract_value!(ev);
                                if val != 0.0 {
                                    precursor.mz = val;
                                }
                            }
                            SELECTED_ION_INT => {
                                precursor.intensity = Some(extract_value!(ev));
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }
                    Some(State::Scan) => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SCAN_START_TIME => {
                                let scan_start_time = extract_value!(ev);
                                let unit = extract!(ev, b"unitAccession");

                                spectrum.scan_start_time = match unit.as_ref() {
                                    UNIT_SECONDS => scan_start_time / 60.0,
                                    UNIT_MINUTES => scan_start_time,
                                    _ => return Err(MzMLError::Malformed),
                                };
                            }
                            ION_INJECTION_TIME => {
                                spectrum.ion_injection_time = extract_value!(ev);
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }

                    _ => {}
                }
            }};
        }

        // A document whose `<mzML>` root is never closed is truncated (e.g. an incomplete
        // download or a cut-off .mzML.gz); reading it as if complete would silently drop
        // spectra. (Bare `<spectrum>` fragments without a root are accepted.)
        let (mut root_opened, mut root_closed) = (false, false);
        loop {
            match reader.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) => {
                    root_opened |= ev.name().into_inner() == b"mzML";
                    match ev.name().into_inner() {
                        b"referenceableParamGroup" => {
                            current_group = Some(extract!(ev, b"id").to_vec())
                        }
                        // `<cvParam ...></cvParam>` and `<referenceableParamGroupRef ...>
                        // </referenceableParamGroupRef>` are legal too
                        b"cvParam" => match &current_group {
                            Some(id) => param_groups
                                .entry(id.clone())
                                .or_default()
                                .push(ev.clone().into_owned()),
                            None => handle_cv!(ev),
                        },
                        b"referenceableParamGroupRef" => {
                            let id = extract!(ev, b"ref").to_vec();
                            let params = param_groups.get(&id).cloned().unwrap_or_default();
                            for param in &params {
                                handle_cv!(param);
                            }
                        }
                        _ => {}
                    }
                    // State transition into child tag
                    state = match (ev.name().into_inner(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
                        (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
                        _ => state,
                    };
                    match ev.name().into_inner() {
                        b"spectrum" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?;
                            spectrum.id = id.to_string();
                            // nothing may leak from the previous spectrum
                            precursor = Precursor::default();
                            iso_window_lo = None;
                            iso_window_hi = None;
                            noise_array.clear();
                        }
                        b"binaryDataArray" => {
                            compression = false;
                            binary_dtype = Dtype::F64;
                            binary_array = None;
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                                let scan = std::str::from_utf8(&scan.value)?;
                                precursor.spectrum_ref = Some(scan.to_string())
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match ev.name().into_inner() {
                    b"cvParam" => match &current_group {
                        Some(id) => param_groups
                            .entry(id.clone())
                            .or_default()
                            .push(ev.clone().into_owned()),
                        None => handle_cv!(ev),
                    },
                    b"referenceableParamGroupRef" => {
                        let id = extract!(ev, b"ref").to_vec();
                        let params = param_groups.get(&id).cloned().unwrap_or_default();
                        for param in &params {
                            handle_cv!(param);
                        }
                    }
                    _ => {}
                },
                Ok(Event::Text(text)) => {
                    if let Some(State::Binary) = state {
                        if let Some(filter) = self.ms_level {
                            if spectrum.ms_level != filter {
                                continue;
                            }
                        }
                        let raw = text.unescape()?;
                        // There are occasionally empty binary data arrays, or unknown CVs
                        if raw.is_empty() || binary_array.is_none() {
                            continue;
                        }
                        let decoded = base64::decode(raw.as_bytes())?;
                        let bytes = match compression {
                            false => &decoded,
                            true => {
                                let mut r = ZlibDecoder::new(decoded.as_slice());
                                let n = r.read_to_end(&mut output_buffer).await?;
                                &output_buffer[..n]
                            }
                        };

                        let width = match binary_dtype {
                            Dtype::F32 => 4,
                            Dtype::F64 => 8,
                        };
                        if bytes.len() % width != 0 {
                            // a truncated/corrupt array would otherwise be silently shortened
                            return Err(MzMLError::CorruptArray(spectrum.id.clone()));
                        }
                        let array = match binary_dtype {
                            Dtype::F32 => {
                                let mut buf: [u8; 4] = [0; 4];
                                bytes
                                    .chunks_exact(4)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f32::from_le_bytes(buf)
                                    })
                                    .collect::<Vec<f32>>()
                            }
                            Dtype::F64 => {
                                let mut buf: [u8; 8] = [0; 8];
                                // `chunks(8)` panicked on a truncated/corrupt array
                                bytes
                                    .chunks_exact(8)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f64::from_le_bytes(buf) as f32
                                    })
                                    .collect::<Vec<f32>>()
                            }
                        };
                        output_buffer.clear();

                        match binary_array {
                            Some(BinaryKind::Intensity) => {
                                spectrum.intensity = array;
                            }
                            Some(BinaryKind::Mz) => {
                                spectrum.mz = array;
                            }
                            Some(BinaryKind::Noise) => {
                                noise_array = array;
                            }
                            None => {}
                        }

                        binary_array = None;
                    }
                }
                Ok(Event::End(ev)) => {
                    root_closed |= ev.name().into_inner() == b"mzML";
                    if ev.name().into_inner() == b"referenceableParamGroup" {
                        current_group = None;
                    }
                    state = match (state, ev.name().into_inner()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
                        (Some(State::Precursor), b"precursor") => {
                            if precursor.mz != 0.0 {
                                precursor.isolation_window = match (iso_window_lo, iso_window_hi) {
                                    (Some(lo), Some(hi)) => Some(Tolerance::Da(-lo, hi)),
                                    _ => None,
                                };
                                spectrum.precursors.push(precursor);
                                precursor = Precursor::default();
                            }
                            Some(State::Spectrum)
                        }
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"spectrum") => {
                            let allow = self
                                .ms_level
                                .as_ref()
                                .map(|&level| level == spectrum.ms_level)
                                .unwrap_or(true);

                            match (allow, self.signal_to_noise) {
                                (true, Some(level))
                                    if level == spectrum.ms_level && !noise_array.is_empty() =>
                                {
                                    spectrum
                                        .intensity
                                        .iter_mut()
                                        .zip(noise_array.iter())
                                        .for_each(|(int, noise)| *int /= noise);
                                    noise_array.clear();
                                    spectra.push(spectrum);
                                }
                                (true, _) => {
                                    spectra.push(spectrum);
                                }
                                (false, _) => {}
                            }
                            spectrum = RawSpectrum::default_with_file_id(self.file_id);
                            None
                        }
                        _ => state,
                    };
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(err) => return Err(err.into()),
            }
            buf.clear();
        }
        if root_opened && !root_closed {
            return Err(MzMLError::Truncated(spectra.len()));
        }
        Ok(spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MzMLError {
    #[error("malformed MzML")]
    Malformed,
    #[error("truncated MzML: file ends before </mzML> (after {0} spectra)")]
    Truncated(usize),
    #[error(
        "corrupt binary array in spectrum {0}: byte count is not a multiple of the value width"
    )]
    CorruptArray(String),
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
    #[error("XML parsing error: {0}")]
    XMLError(#[from] quick_xml::Error),
    #[error("io error: {0}")]
    IOError(#[from] std::io::Error),
    #[error("utf8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("error parsing float: {0}")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("error parsing int: {0}")]
    IntError(#[from] std::num::ParseIntError),
    #[error("error decoding base64: {0}")]
    Base64Error(#[from] base64::DecodeError),
}

#[cfg(test)]
mod test {
    use sage_core::{mass::Tolerance, spectrum::Representation};

    use super::{MzMLError, MzMLReader};

    #[tokio::test]
    async fn referenceable_param_groups_are_resolved() {
        let full = include_str!("../../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
        let ms_level = r#"<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>"#;
        let centroid =
            r#"<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>"#;
        assert!(full.contains(ms_level) && full.contains(centroid));
        // move MS level and centroid flag of the spectrum into a param group
        let grouped = full
            .replacen(ms_level, r#"<referenceableParamGroupRef ref="ms2"/>"#, 1)
            .replacen(centroid, "", 1)
            .replacen(
                r#"<referenceableParamGroupList count="1">"#,
                &format!(
                    r#"<referenceableParamGroupList count="2">
                    <referenceableParamGroup id="ms2">{ms_level}{centroid}</referenceableParamGroup>"#
                ),
                1,
            );
        let expected = MzMLReader::with_file_id(0)
            .parse(full.as_bytes())
            .await
            .unwrap();
        let parsed = MzMLReader::with_file_id(0)
            .parse(grouped.as_bytes())
            .await
            .unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].ms_level, 2);
        assert_eq!(parsed[0].representation, Representation::Centroid);
        assert_eq!(parsed[0].mz, expected[0].mz);
        assert_eq!(parsed[0].precursors.len(), expected[0].precursors.len());
    }

    const FIXTURE: &str = include_str!("../../../tests/LQSRPAAPPAPGPGQLTLR.mzML");

    #[tokio::test]
    async fn numpress_is_an_error_not_garbage() {
        let numpress = FIXTURE.replacen("MS:1000574", "MS:1002312", 1);
        let result = MzMLReader::with_file_id(0).parse(numpress.as_bytes()).await;
        assert!(
            matches!(result, Err(MzMLError::UnsupportedCV(_))),
            "{result:?}"
        );
    }

    #[tokio::test]
    async fn unknown_cv_params_do_not_drop_arrays() {
        let expected = MzMLReader::with_file_id(0)
            .parse(FIXTURE.as_bytes())
            .await
            .unwrap();
        // an unrelated cvParam after the array-type cvParam used to reset the array type,
        // silently dropping the m/z array
        let extra = r#"<cvParam cvRef="MS" accession="MS:1000786" name="non-standard data array" value=""/>"#;
        let with_extra = FIXTURE.replace(
            r#"unitName="m/z"/>"#,
            &format!(r#"unitName="m/z"/>{extra}"#),
        );
        assert_ne!(with_extra, FIXTURE);
        let parsed = MzMLReader::with_file_id(0)
            .parse(with_extra.as_bytes())
            .await
            .unwrap();
        assert!(!parsed[0].mz.is_empty());
        assert_eq!(parsed[0].mz, expected[0].mz);
        assert_eq!(parsed[0].intensity, expected[0].intensity);
    }

    #[tokio::test]
    async fn isolation_window_does_not_leak_into_next_spectrum() {
        let start = FIXTURE.find("<spectrum ").unwrap();
        let end = FIXTURE.find("</spectrum>").unwrap() + "</spectrum>".len();
        let without_window = FIXTURE[start..end]
            .lines()
            .filter(|l| !l.contains("MS:1000828") && !l.contains("MS:1000829"))
            .collect::<Vec<_>>()
            .join("\n")
            .replace("scan=30069", "scan=30070");
        let two = format!("{}{}{}", &FIXTURE[..end], without_window, &FIXTURE[end..]);
        let parsed = MzMLReader::with_file_id(0)
            .parse(two.as_bytes())
            .await
            .unwrap();
        assert_eq!(parsed.len(), 2);
        assert!(parsed[0].precursors[0].isolation_window.is_some());
        assert!(parsed[1].precursors[0].isolation_window.is_none());
    }

    #[tokio::test]
    async fn param_groups_in_all_element_forms() {
        let ms_level = r#"<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>"#;
        let centroid =
            r#"<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>"#;
        let expand = |s: &str| s.replacen("/>", "></cvParam>", 1);
        for (group_form, ref_form) in [(false, false), (true, false), (false, true), (true, true)] {
            let (ms, ce) = match group_form {
                true => (expand(ms_level), expand(centroid)),
                false => (ms_level.to_string(), centroid.to_string()),
            };
            let reference = match ref_form {
                true => r#"<referenceableParamGroupRef ref="ms2"></referenceableParamGroupRef>"#,
                false => r#"<referenceableParamGroupRef ref="ms2"/>"#,
            };
            let grouped = FIXTURE
                .replacen(ms_level, reference, 1)
                .replacen(centroid, "", 1)
                .replacen(
                    r#"<referenceableParamGroupList count="1">"#,
                    &format!(
                        r#"<referenceableParamGroupList count="2">
                        <referenceableParamGroup id="ms2">{ms}{ce}</referenceableParamGroup>"#
                    ),
                    1,
                );
            let parsed = MzMLReader::with_file_id(0)
                .parse(grouped.as_bytes())
                .await
                .unwrap();
            assert_eq!(parsed[0].ms_level, 2, "group {group_form} ref {ref_form}");
            assert_eq!(parsed[0].representation, Representation::Centroid);
        }
    }

    #[tokio::test]
    async fn truncated_binary_array_is_an_error() {
        // drop the last base64 quantum of the first array: 3 bytes -> not a multiple of 4/8
        let start = FIXTURE.find("<binary>").unwrap() + "<binary>".len();
        let end = FIXTURE[start..].find("</binary>").unwrap() + start;
        let b64 = &FIXTURE[start..end];
        let mut bytes = base64::decode(b64).unwrap();
        // the fixture's arrays are zlib-compressed: decompress, truncate, recompress
        use std::io::Read;
        let mut raw = Vec::new();
        flate2::read::ZlibDecoder::new(bytes.as_slice())
            .read_to_end(&mut raw)
            .unwrap();
        raw.truncate(raw.len() - 3);
        let mut enc = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut enc, &raw).unwrap();
        bytes = enc.finish().unwrap();
        let corrupt = format!(
            "{}{}{}",
            &FIXTURE[..start],
            base64::encode(bytes),
            &FIXTURE[end..]
        );
        let result = MzMLReader::with_file_id(0).parse(corrupt.as_bytes()).await;
        assert!(
            matches!(result, Err(MzMLError::CorruptArray(_))),
            "{result:?}"
        );
    }

    #[tokio::test]
    async fn truncated_file_is_an_error() {
        let full = include_str!("../../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
        let complete = MzMLReader::with_file_id(0)
            .parse(full.as_bytes())
            .await
            .unwrap();
        assert_eq!(complete.len(), 1);

        // cut in the header, inside the spectrum, and after it but before </mzML>
        for cut in [
            full.find("<run").unwrap(),
            full.len() / 2,
            full.find("</spectrumList>").unwrap(),
        ] {
            let result = MzMLReader::with_file_id(0)
                .parse(full[..cut].as_bytes())
                .await;
            assert!(
                matches!(
                    result,
                    Err(MzMLError::Truncated(_)) | Err(MzMLError::XMLError(_))
                ),
                "cut at {cut}: {result:?}"
            );
        }
    }

    #[tokio::test]
    async fn parse_spectrum_issue_78() -> Result<(), MzMLError> {
        let s = r#"
        <spectrum id="spectrum=2442" index="286" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="147.290603637695"/>
            <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="769.255798339843977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <userParam name="preset scan configuration" type="xsd:string" value="2"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                            <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();

        assert_eq!(s.id, "spectrum=2442");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_117() -> Result<(), MzMLError> {
        // The issue was that some converters write the ion mobility as part of the selected ion (as in the last test)
        // and some write it as part of the scan, as in this test. This test checks that it can be read
        // fom the scan section.
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });

        // The rest of these assertions just make sure the integrity of the spectrum is maintained
        assert_eq!(s.id, "spectrum=8678309");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_210() -> Result<(), MzMLError> {
        // Handle cases where both isolation window target m/z is set and different than selected ion m/z
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );

        // Check different ordering of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );

        // Check fallback keeping iso window m/z of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.75) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        Ok(())
    }
}
