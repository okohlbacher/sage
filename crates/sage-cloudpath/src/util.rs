use crate::{read_and_execute, tdf::BrukerProcessingConfig, Error};
use sage_core::spectrum::RawSpectrum;
use serde::Serialize;
use tokio::io::AsyncReadExt;
use url::Url;

#[derive(Debug, PartialEq, Eq)]
pub enum FileFormat {
    MzML,
    MGF,
    TDF,
    Unidentified,
}

impl FileFormat {
    /// Does this file format support parallel reading?
    /// By this I mean that there is 'within' file parallelism that
    /// would make it faster to read than reading mutiple files in
    /// parallel. (is giving 4 cores to read 1 file 4 times, faster
    /// than giving 1 cores to read 1 file and read 4 at the same time)
    pub fn within_file_parallel(&self) -> bool {
        match self {
            FileFormat::MzML => false,
            FileFormat::MGF => false,
            FileFormat::TDF => true,
            FileFormat::Unidentified => false,
        }
    }
}

impl From<&str> for FileFormat {
    fn from(s: &str) -> Self {
        // ignore the query of URLs such as presigned S3 links (`...file.mzML?X-Amz-...`)
        let path = s.split('?').next().unwrap_or(s);
        let path_lower = path.to_lowercase();
        if path_lower.ends_with(".mgf.gz") || path_lower.ends_with(".mgf") {
            FileFormat::MGF
        } else if is_bruker(&path_lower) {
            FileFormat::TDF
        } else if path_lower.ends_with(".mzml.gz") || path_lower.ends_with(".mzml") {
            FileFormat::MzML
        } else {
            FileFormat::Unidentified
        }
    }
}

const BRUKER_EXTENSIONS: [&str; 5] = [".d", ".tdf", ".tdf_bin", "ms2", "raw"];

fn is_bruker(path: &str) -> bool {
    BRUKER_EXTENSIONS.iter().any(|ext| {
        if path.ends_with(std::path::MAIN_SEPARATOR) {
            path.strip_suffix(std::path::MAIN_SEPARATOR)
                .unwrap()
                .ends_with(ext)
        } else {
            path.ends_with(ext)
        }
    })
}

pub fn read_spectra(
    url: &Url,
    file_id: usize,
    sn: Option<u8>,
    bruker_processor: BrukerProcessingConfig,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    match FileFormat::from(url.as_ref()) {
        FileFormat::MzML => read_mzml_levels(url, file_id, sn, !requires_ms1),
        FileFormat::MGF => read_mgf(url, file_id),
        FileFormat::TDF => read_tdf(url, file_id, bruker_processor, requires_ms1),
        FileFormat::Unidentified => Err(Error::UnsupportedFormat(url.to_string())),
    }
}

pub fn read_mzml(
    url: &Url,
    file_id: usize,
    signal_to_noise: Option<u8>,
) -> Result<Vec<RawSpectrum>, Error> {
    read_mzml_levels(url, file_id, signal_to_noise, false)
}

/// Like [`read_mzml`], optionally without decoding/keeping MS1 spectra (which only LFQ
/// uses; they are ~90% of the peaks of typical DDA files)
pub fn read_mzml_levels(
    url: &Url,
    file_id: usize,
    signal_to_noise: Option<u8>,
    skip_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(url, |bf| async move {
        Ok(crate::mzml::MzMLReader::with_file_id(file_id)
            .set_signal_to_noise(signal_to_noise)
            .set_skip_ms1(skip_ms1)
            .parse(bf)
            .await?)
    })
}

/// [`read_spectra`], then `f` applied to every spectrum. For Bruker .d files `f` runs
/// while the file is read, so only a block of unprocessed spectra is alive at a time
/// instead of the whole file's (~1.5 GB for a 60 min ddaPASEF run).
pub fn read_processed<T: Send>(
    url: &Url,
    file_id: usize,
    sn: Option<u8>,
    bruker_processor: BrukerProcessingConfig,
    requires_ms1: bool,
    f: impl Fn(RawSpectrum) -> T + Sync,
) -> Result<Vec<T>, Error> {
    use rayon::prelude::*;
    match FileFormat::from(url.as_ref()) {
        FileFormat::TDF => read_tdf_with(url, file_id, bruker_processor, requires_ms1, f),
        _ => Ok(
            read_spectra(url, file_id, sn, bruker_processor, requires_ms1)?
                .into_par_iter()
                .map(&f)
                .collect(),
        ),
    }
}

pub fn read_tdf(
    url: &Url,
    file_id: usize,
    bruker_spectrum_processor: BrukerProcessingConfig,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    read_tdf_with(url, file_id, bruker_spectrum_processor, requires_ms1, |s| s)
}

fn read_tdf_with<T: Send>(
    url: &Url,
    file_id: usize,
    bruker_spectrum_processor: BrukerProcessingConfig,
    requires_ms1: bool,
    f: impl Fn(RawSpectrum) -> T + Sync,
) -> Result<Vec<T>, Error> {
    if url.scheme() != "file" {
        log::error!("Bruker files must be local: {}", url);
        return Err(Error::InvalidUri);
    }

    let path = url.to_file_path().map_err(|_| Error::InvalidUri)?;
    let res =
        crate::tdf::TdfReader.parse(&path, file_id, bruker_spectrum_processor, requires_ms1, f);
    match res {
        Ok(t) => Ok(t),
        Err(e) => Err(Error::TDF(e)),
    }
}

pub fn read_mgf(url: &Url, file_id: usize) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(url, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        let res = crate::mgf::MgfReader::with_file_id(file_id).parse(contents);
        match res {
            Ok(m) => Ok(m),
            Err(e) => Err(Error::MGF(e)),
        }
    })
}

pub fn read_fasta<S>(
    url: &Url,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    read_and_execute(url, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        Ok(sage_core::fasta::Fasta::parse(
            contents,
            decoy_tag.as_ref(),
            generate_decoys,
        ))
    })
}

pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents).await?;
        Ok(serde_json::from_str(&contents)?)
    })
}

/// Send telemetry data
pub fn send_data<T>(url: &str, data: &T) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
where
    T: Serialize,
{
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let client = reqwest::ClientBuilder::default().https_only(true).build()?;
        let res = client.post(url).json(data).send().await?;
        res.error_for_status()?;
        Ok(())
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn multi_member_gzip_is_read_completely() {
        use std::io::Write;
        let fixture = include_str!("../../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
        let (head, tail) = fixture.split_at(fixture.len() / 2);
        let mut file = Vec::new();
        for part in [head, tail] {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            enc.write_all(part.as_bytes()).unwrap();
            file.extend(enc.finish().unwrap());
        }
        let path =
            std::env::temp_dir().join(format!("sage-multimember-{}.mzML.gz", std::process::id()));
        std::fs::write(&path, file).unwrap();
        let url = crate::to_url(path.to_str().unwrap()).unwrap();
        let spectra = read_mzml(&url, 0, None);
        std::fs::remove_file(&path).ok();
        assert_eq!(spectra.unwrap().len(), 1);
    }

    #[test]
    fn test_identify_format() {
        assert_eq!(FileFormat::from("foo.mzml"), FileFormat::MzML);
        assert_eq!(FileFormat::from("foo.mzML"), FileFormat::MzML);
        assert_eq!(FileFormat::from("foo.mgf"), FileFormat::MGF);
        assert_eq!(FileFormat::from("foo.mgf.gz"), FileFormat::MGF);
        assert_eq!(FileFormat::from("foo.tdf"), FileFormat::TDF);
        assert_eq!(FileFormat::from("./tomato/foo.d"), FileFormat::TDF);
        assert_eq!(FileFormat::from("./tomato/foo.d/"), FileFormat::TDF);
        assert_eq!(
            FileFormat::from("s3://bucket/foo.mzML?X-Amz-Signature=abc"),
            FileFormat::MzML
        );
        assert_eq!(FileFormat::from("foo.mzXML"), FileFormat::Unidentified);
    }
}
