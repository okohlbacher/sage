use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::Url;
use sage_core::scoring::ScoreType;
use sage_core::{
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Clone)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub version: String,
    pub database: Parameters,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub predict_rt: bool,
    pub mzml_paths: Vec<Url>,
    pub output_paths: Vec<Url>,
    pub bruker_config: BrukerProcessingConfig,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f32,
    /// Self-trained, cross-fitted fragment-ion model: `ion_llr` and `ion_explained` columns
    pub ion_model: bool,

    #[serde(skip_serializing)]
    pub output_directory: Url,

    #[serde(skip_serializing)]
    pub write_pin: bool,

    #[serde(skip_serializing)]
    pub write_report: bool,

    #[serde(skip_serializing)]
    pub annotate_matches: bool,

    pub score_type: ScoreType,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Builder,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub report_psms: Option<usize>,
    pub chimera: Option<bool>,
    pub wide_window: Option<bool>,
    pub min_peaks: Option<usize>,
    pub max_peaks: Option<MaxPeaks>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: Option<u16>,
    pub precursor_charge: Option<(u8, u8)>,
    pub override_precursor_charge: Option<bool>,
    pub isotope_errors: Option<(i8, i8)>,
    pub deisotope: Option<bool>,
    pub quant: Option<QuantOptions>,
    pub predict_rt: Option<bool>,
    pub output_directory: Option<String>,
    pub mzml_paths: Option<Vec<String>>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,
    pub ion_model: Option<bool>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub score_type: Option<ScoreType>,
}

/// `max_peaks` as written in the configuration: a number of peaks, or `"auto"`.
/// Not set means [`MaxPeaks::DEFAULT`] peaks (see [`MaxPeaks::default`]).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MaxPeaks {
    /// Depends on the fragment tolerance and the search mode (see [`MaxPeaks::resolve`])
    Auto,
    Count(usize),
}

impl Default for MaxPeaks {
    /// `max_peaks` not set: a fixed [`MaxPeaks::DEFAULT`] peaks
    fn default() -> Self {
        MaxPeaks::Count(Self::DEFAULT)
    }
}

impl MaxPeaks {
    /// Peaks kept per MS2 spectrum when `max_peaks` is not set
    pub const DEFAULT: usize = 150;
    /// Peaks kept per MS2 spectrum by `"auto"` for high-resolution fragment spectra
    pub const AUTO_HIGH_RES: usize = 400;
    /// Peaks kept per MS2 spectrum by `"auto"` for low-resolution fragment spectra
    /// (e.g. ion trap CID searched at 0.5 Da)
    pub const AUTO_LOW_RES: usize = 80;
    /// Peaks kept per MS2 spectrum by `"auto"` for low-resolution fragment spectra when TMT
    /// reporter ions are quantified from the MS2 spectra (`quant.tmt_settings.level` 2): the
    /// fixed default, because 80 peaks drop reporter ions
    pub const AUTO_LOW_RES_TMT_MS2: usize = 150;
    /// Peaks kept per MS2 spectrum by `"auto"` in a wide-window search (`wide_window: true`),
    /// whatever the resolution: the fixed default. With 400 peaks, wrong candidates from the
    /// wide isolation window win more often (3-46% fewer PSMs on high-resolution data)
    pub const AUTO_WIDE_WINDOW: usize = 150;
    /// A Da fragment tolerance reaching this value (on at least one side) means low resolution
    pub const LOW_RES_DA: f32 = 0.1;

    /// Does `"auto"` treat this fragment tolerance as low resolution? Only Da tolerances
    /// of at least [`MaxPeaks::LOW_RES_DA`]; ppm and pct tolerances are high resolution.
    pub fn low_resolution(fragment_tol: &Tolerance) -> bool {
        match fragment_tol {
            Tolerance::Da(lo, hi) => lo.abs().max(hi.abs()) >= Self::LOW_RES_DA,
            Tolerance::Ppm(..) | Tolerance::Pct(..) => false,
        }
    }

    /// The number of peaks to keep per MS2 spectrum: an explicit count as given; `"auto"`
    /// gives [`MaxPeaks::AUTO_WIDE_WINDOW`] if `wide_window` (a wide-window search), else
    /// [`MaxPeaks::AUTO_HIGH_RES`] for a ppm or pct fragment tolerance or a narrow Da
    /// tolerance (such as ±0.02 Da); for a Da tolerance reaching [`MaxPeaks::LOW_RES_DA`] it
    /// gives [`MaxPeaks::AUTO_LOW_RES`], or [`MaxPeaks::AUTO_LOW_RES_TMT_MS2`] if `tmt_ms2`
    /// (TMT reporter ions are quantified from the MS2 spectra).
    /// [`Input::build`] then raises a resolved `"auto"` to `min_peaks` if it is lower.
    pub fn resolve(self, fragment_tol: &Tolerance, tmt_ms2: bool, wide_window: bool) -> usize {
        match self {
            MaxPeaks::Count(n) => n,
            MaxPeaks::Auto if wide_window => Self::AUTO_WIDE_WINDOW,
            MaxPeaks::Auto if !Self::low_resolution(fragment_tol) => Self::AUTO_HIGH_RES,
            MaxPeaks::Auto if tmt_ms2 => Self::AUTO_LOW_RES_TMT_MS2,
            MaxPeaks::Auto => Self::AUTO_LOW_RES,
        }
    }
}

impl<'de> Deserialize<'de> for MaxPeaks {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Visitor;
        impl serde::de::Visitor<'_> for Visitor {
            type Value = MaxPeaks;

            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a non-negative integer or \"auto\"")
            }

            fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<MaxPeaks, E> {
                usize::try_from(v)
                    .map(MaxPeaks::Count)
                    .map_err(|_| E::invalid_value(serde::de::Unexpected::Unsigned(v), &self))
            }

            fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<MaxPeaks, E> {
                u64::try_from(v)
                    .map_err(|_| E::invalid_value(serde::de::Unexpected::Signed(v), &self))
                    .and_then(|v| self.visit_u64(v))
            }

            fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<MaxPeaks, E> {
                match v {
                    "auto" => Ok(MaxPeaks::Auto),
                    _ => Err(E::invalid_value(serde::de::Unexpected::Str(v), &self)),
                }
            }
        }
        deserializer.deserialize_any(Visitor)
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    pub peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    pub integration: Option<sage_core::lfq::IntegrationStrategy>,
    pub spectral_angle: Option<f64>,
    pub ppm_tolerance: Option<f32>,
    pub mobility_pct_tolerance: Option<f32>,
    pub combine_charge_states: Option<bool>,
    pub peptide_q_value: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
            peptide_q_value: value.peptide_q_value.unwrap_or(default.peptide_q_value),
            mobility_pct_tolerance: value
                .mobility_pct_tolerance
                .unwrap_or(default.mobility_pct_tolerance),
            combine_charge_states: value
                .combine_charge_states
                .unwrap_or(default.combine_charge_states),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance > 4.0 {
            log::warn!("lfq_settings.mobility_pct_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.mobility_pct_tolerance is smaller than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }
        if settings.peptide_q_value > 0.01 {
            log::info!("lfq_settings.peptide_q_value is higher than expected, expect increased runtime and memory usage");
        }
        if settings.peptide_q_value < 0.01 {
            log::warn!("lfq_settings.peptide_q_value is lower than expected, not all identified peptides will have MS1 intensities extracted");
        }

        settings
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    pub level: Option<u8>,
    pub sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize, Debug)]
pub struct TmtSettings {
    pub level: u8,
    pub sn: bool,
}

impl From<TmtOptions> for TmtSettings {
    fn from(value: TmtOptions) -> Self {
        let default = Self::default();
        Self {
            level: value.level.unwrap_or(default.level),
            sn: value.sn.unwrap_or(default.sn),
        }
    }
}

impl Default for TmtSettings {
    fn default() -> Self {
        Self {
            level: 3,
            sn: false,
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct QuantOptions {
    pub tmt: Option<Isobaric>,
    #[serde(rename = "tmt_settings")]
    pub tmt_options: Option<TmtOptions>,

    pub lfq: Option<bool>,
    #[serde(rename = "lfq_settings")]
    pub lfq_options: Option<LfqOptions>,
}

#[derive(Serialize, Default, Clone)]
pub struct QuantSettings {
    pub tmt: Option<Isobaric>,
    pub tmt_settings: TmtSettings,
    pub lfq: bool,
    pub lfq_settings: LfqSettings,
}

impl From<QuantOptions> for QuantSettings {
    fn from(value: QuantOptions) -> Self {
        Self {
            tmt: value.tmt,
            tmt_settings: value.tmt_options.map(Into::into).unwrap_or_default(),

            lfq: value.lfq.unwrap_or(false),
            lfq_settings: value.lfq_options.map(Into::into).unwrap_or_default(),
        }
    }
}

impl Input {
    pub fn from_arguments(matches: ArgMatches) -> anyhow::Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(output_directory) = matches.get_one::<String>("output_directory") {
            input.output_directory = Some(output_directory.into());
        }
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.fasta = Some(fasta.into());
        }
        if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
            input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
        }

        if let Some(write_pin) = matches.get_one::<bool>("write-pin").copied() {
            input.write_pin = Some(write_pin);
        }

        if let Some(write_report) = matches.get_one::<bool>("write-report").copied() {
            input.write_report = Some(write_report);
        }

        if let Some(annotate_matches) = matches.get_one::<bool>("annotate-matches").copied() {
            input.annotate_matches = Some(annotate_matches);
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some(),
            "`database.fasta` must be set. For more information try '--help'"
        );
        ensure!(
            input
                .mzml_paths
                .as_ref()
                .map(|p| p.len())
                .unwrap_or_default()
                > 0,
            "`mzml_paths` must be set. For more information try '--help'"
        );

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        sage_cloudpath::util::read_json(path).map_err(anyhow::Error::from)
    }

    fn check_mass_tolerances(tolerance: &Tolerance) {
        let (lo, hi) = match tolerance {
            Tolerance::Ppm(lo, hi) => (*lo, *hi),
            Tolerance::Pct(lo, hi) => {
                log::warn!(
                    "Pct tolerances are very rarely used for mass tolerances, did you mean ppm?"
                );
                (*lo, *hi)
            }
            Tolerance::Da(lo, hi) => (*lo, *hi),
        };
        if hi.abs() > lo.abs() {
            log::warn!(
                "Tolerances are applied to experimental masses, not theoretical: [{}, {}]",
                lo,
                hi
            );
        }
        if lo > 0.0 {
            log::warn!(
                "The `left` tolerance should probably be negative, for example: [{}, {}]",
                -lo,
                hi.abs()
            )
        }
        if hi < 0.0 {
            log::warn!(
                "The `right` tolerance should probably be positive, for example: [{}, {}]",
                -lo.abs(),
                hi
            )
        }
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        if let Some(isotope_errors) = self.isotope_errors {
            if isotope_errors.0 > isotope_errors.1 {
                log::error!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
                std::process::exit(1);
            }
        }
        if let Some(charges) = self.precursor_charge {
            if charges.0 > charges.1 {
                log::error!(
                    "Precursor charges should be specified [low, high], user provided: [{}, {}]",
                    charges.0,
                    charges.1
                );
                std::process::exit(1);
            }
        }

        if !self.predict_rt.unwrap_or(true)
            && self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false)
        {
            log::warn!(
                "`predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`"
            );
            self.predict_rt = Some(true);
        }

        let mzml_paths = self
            .mzml_paths
            .expect("'mzml_paths' must be provided!")
            .iter()
            .map(|s| sage_cloudpath::to_url(s))
            .collect::<Result<Vec<_>, _>>()?;

        // Fail before the (possibly long) database build, not halfway through the run
        let unsupported = mzml_paths
            .iter()
            .filter(|url| !sage_cloudpath::FileFormat::from(url.as_str()).supported())
            .map(|url| url.to_string())
            .collect::<Vec<_>>();
        if !unsupported.is_empty() {
            anyhow::bail!(
                "unsupported input file format (expected .mzML[.gz], .mgf[.gz], .mzpeak (feature `mzpeak`) or a Bruker .d): {}",
                unsupported.join(", ")
            );
        }

        let output_directory = match self.output_directory {
            Some(path) => {
                match sage_cloudpath::try_parse_url(&path) {
                    Some(mut url) => {
                        // Valid URL, might still be a local directory that doesn't exist
                        if url.scheme() == "file" {
                            let path = url.to_file_path().expect("url scheme is file");
                            std::fs::create_dir_all(path)?;
                        }

                        if !url.path().ends_with("/") {
                            url.set_path(&format!("{}/", url.path()));
                        }
                        url
                    }
                    None => {
                        // Treat as a local path (covers Windows `C:\...` which
                        // otherwise parses as a URL with scheme `c`).
                        let path = std::path::Path::new(&path);
                        std::fs::create_dir_all(path)?;
                        Url::from_directory_path(path.canonicalize()?).expect("valid path")
                    }
                }
            }
            None => {
                let dir = std::env::current_dir()?;
                Url::from_directory_path(dir).expect("valid path")
            }
        };

        let score_type = self.score_type.unwrap_or(ScoreType::SageHyperScore);

        let quant: QuantSettings = self.quant.map(Into::into).unwrap_or_default();
        // TMT reporter ions quantified from the MS2 spectra must survive the peak cap
        let tmt_ms2 = quant.tmt.is_some() && quant.tmt_settings.level == 2;
        let wide_window = self.wide_window.unwrap_or(false);
        let min_peaks = self.min_peaks.unwrap_or(15);
        let max_peaks = self.max_peaks.unwrap_or_default();
        let mut max_peaks_resolved = max_peaks.resolve(&self.fragment_tol, tmt_ms2, wide_window);
        if max_peaks == MaxPeaks::Auto {
            // A spectrum is searched only if at least `min_peaks` peaks survive the top-N cut
            // (`Runner::searchable`): "auto" never keeps fewer than that
            let raised = max_peaks_resolved < min_peaks;
            max_peaks_resolved = max_peaks_resolved.max(min_peaks);
            let low_resolution = MaxPeaks::low_resolution(&self.fragment_tol);
            let reason = if wide_window {
                "wide_window: the fixed default, whatever the resolution".to_string()
            } else {
                format!(
                    "{:?} fragment tolerance: {} resolution{}",
                    self.fragment_tol,
                    if low_resolution { "low" } else { "high" },
                    if low_resolution && tmt_ms2 {
                        ", TMT reporter ions quantified from MS2"
                    } else {
                        ""
                    },
                )
            };
            log::info!(
                "max_peaks: auto -> {} peaks per MS2 spectrum ({}{})",
                max_peaks_resolved,
                reason,
                if raised {
                    format!("; raised to min_peaks {}", min_peaks)
                } else {
                    String::new()
                },
            );
        } else if max_peaks_resolved < min_peaks {
            log::warn!(
                "max_peaks {} is below min_peaks {}: no MS2 spectrum can be searched",
                max_peaks_resolved,
                min_peaks
            );
        }

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            quant,
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: max_peaks_resolved,
            min_peaks,
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            annotate_matches: self.annotate_matches.unwrap_or(false),
            precursor_charge: self.precursor_charge.unwrap_or((2, 4)),
            override_precursor_charge: self.override_precursor_charge.unwrap_or(false),
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            wide_window,
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
            write_pin: self.write_pin.unwrap_or(false),
            bruker_config: self.bruker_config.unwrap_or_default(),
            write_report: self.write_report.unwrap_or(false),
            protein_grouping: self.protein_grouping.unwrap_or(true),
            protein_grouping_peptide_fdr: self.protein_grouping_peptide_fdr.unwrap_or(0.01),
            ion_model: self.ion_model.unwrap_or(false),
            score_type,
        })
    }
}

#[cfg(test)]
mod test {

    use super::{Input, MaxPeaks, Search};
    use sage_core::{database::EnzymeBuilder, enzyme::EnzymeParameters, mass::Tolerance};

    /// The resolved search parameters for a minimal configuration with `extra` keys added
    fn resolved(
        fragment_tol: serde_json::Value,
        extra: serde_json::Value,
    ) -> anyhow::Result<Search> {
        let mut config = serde_json::json!({
            "database": { "fasta": "proteins.fasta" },
            "precursor_tol": { "ppm": [-10.0, 10.0] },
            "fragment_tol": fragment_tol,
            // a URL is not resolved on disk (a local path must exist)
            "mzml_paths": ["s3://bucket/spectra.mzML"],
        });
        for (key, value) in extra.as_object().expect("extra keys are an object") {
            config[key] = value.clone();
        }
        let input: Input = serde_json::from_value(config)?;
        input.build()
    }

    /// `max_peaks` of the resolved search parameters for a minimal configuration
    fn resolved_max_peaks(
        fragment_tol: serde_json::Value,
        max_peaks: Option<serde_json::Value>,
    ) -> anyhow::Result<usize> {
        resolved_max_peaks_quant(fragment_tol, max_peaks, None)
    }

    /// As [`resolved_max_peaks`], with a `quant` section
    fn resolved_max_peaks_quant(
        fragment_tol: serde_json::Value,
        max_peaks: Option<serde_json::Value>,
        quant: Option<serde_json::Value>,
    ) -> anyhow::Result<usize> {
        let mut extra = serde_json::json!({});
        if let Some(max_peaks) = max_peaks {
            extra["max_peaks"] = max_peaks;
        }
        if let Some(quant) = quant {
            extra["quant"] = quant;
        }
        Ok(resolved(fragment_tol, extra)?.max_peaks)
    }

    #[test]
    fn max_peaks_default_is_150() -> anyhow::Result<()> {
        use serde_json::json;
        let tmt_ms2 = json!({ "tmt": "Tmt16", "tmt_settings": { "level": 2 } });
        assert_eq!(MaxPeaks::default(), MaxPeaks::Count(150));
        assert_eq!(MaxPeaks::DEFAULT, 150);
        // not set and null are the fixed default, whatever the tolerance, TMT or search mode
        for tol in [
            json!({ "ppm": [-20.0, 20.0] }),
            json!({ "da": [-0.5, 0.5] }),
            json!({ "da": [-0.02, 0.02] }),
            json!({ "pct": [-0.002, 0.002] }),
        ] {
            for max_peaks in [None, Some(json!(null))] {
                assert_eq!(resolved_max_peaks(tol.clone(), max_peaks.clone())?, 150);
                assert_eq!(
                    resolved_max_peaks_quant(
                        tol.clone(),
                        max_peaks.clone(),
                        Some(tmt_ms2.clone())
                    )?,
                    150
                );
                let mut extra = json!({ "wide_window": true });
                if let Some(max_peaks) = max_peaks {
                    extra["max_peaks"] = max_peaks;
                }
                assert_eq!(resolved(tol.clone(), extra)?.max_peaks, 150);
            }
        }
        Ok(())
    }

    #[test]
    fn max_peaks_auto_resolves_by_fragment_tolerance() -> anyhow::Result<()> {
        use serde_json::json;
        let ppm = json!({ "ppm": [-20.0, 20.0] });
        let da = json!({ "da": [-0.5, 0.5] });
        let narrow_da = json!({ "da": [-0.02, 0.02] });
        let pct = json!({ "pct": [-0.002, 0.002] });

        // "auto": 400 for ppm, pct and narrow Da, 80 for wide Da
        let auto = Some(json!("auto"));
        assert_eq!(resolved_max_peaks(ppm.clone(), auto.clone())?, 400);
        assert_eq!(resolved_max_peaks(pct.clone(), auto.clone())?, 400);
        assert_eq!(resolved_max_peaks(narrow_da.clone(), auto.clone())?, 400);
        assert_eq!(resolved_max_peaks(da.clone(), auto)?, 80);
        // an explicit number is used as given, whatever the tolerance
        for tol in [&ppm, &da, &narrow_da, &pct] {
            assert_eq!(resolved_max_peaks(tol.clone(), Some(json!(150)))?, 150);
            assert_eq!(resolved_max_peaks(tol.clone(), Some(json!(0)))?, 0);
        }

        let auto = |tol: Tolerance| MaxPeaks::Auto.resolve(&tol, false, false);
        assert_eq!(auto(Tolerance::Ppm(-5.0, 5.0)), MaxPeaks::AUTO_HIGH_RES);
        assert_eq!(auto(Tolerance::Ppm(-500.0, 500.0)), MaxPeaks::AUTO_HIGH_RES);
        assert_eq!(auto(Tolerance::Da(-0.05, 0.05)), MaxPeaks::AUTO_HIGH_RES);
        assert_eq!(auto(Tolerance::Da(-0.099, 0.099)), MaxPeaks::AUTO_HIGH_RES);
        // low resolution from +-0.1 Da on, on either side
        assert_eq!(auto(Tolerance::Da(-0.1, 0.1)), MaxPeaks::AUTO_LOW_RES);
        assert_eq!(auto(Tolerance::Da(-1.0005, 1.0005)), MaxPeaks::AUTO_LOW_RES);
        assert_eq!(auto(Tolerance::Da(-0.5, 0.0)), MaxPeaks::AUTO_LOW_RES);
        assert_eq!(auto(Tolerance::Da(-0.01, 0.4)), MaxPeaks::AUTO_LOW_RES);
        assert_eq!(
            MaxPeaks::Count(42).resolve(&Tolerance::Da(-0.5, 0.5), false, false),
            42
        );
        Ok(())
    }

    #[test]
    fn max_peaks_auto_keeps_tmt_reporters_of_ms2() -> anyhow::Result<()> {
        use serde_json::json;
        let ppm = json!({ "ppm": [-20.0, 20.0] });
        let da = json!({ "da": [-0.5, 0.5] });
        let narrow_da = json!({ "da": [-0.02, 0.02] });
        let tmt_ms2 = json!({ "tmt": "Tmt16", "tmt_settings": { "level": 2 } });
        let tmt_ms3 = json!({ "tmt": "Tmt16", "tmt_settings": { "level": 3 } });
        let tmt_default = json!({ "tmt": "Tmt16" }); // level 3
        let lfq = json!({ "lfq": true });
        let auto = |tol: &serde_json::Value, quant: &serde_json::Value| {
            resolved_max_peaks_quant(tol.clone(), Some(json!("auto")), Some(quant.clone()))
        };

        // low resolution with reporter ions read from MS2: the fixed default
        assert_eq!(auto(&da, &tmt_ms2)?, MaxPeaks::AUTO_LOW_RES_TMT_MS2);
        // reporter ions from MS3 (not capped), or no TMT: unchanged
        assert_eq!(auto(&da, &tmt_ms3)?, MaxPeaks::AUTO_LOW_RES);
        assert_eq!(auto(&da, &tmt_default)?, MaxPeaks::AUTO_LOW_RES);
        assert_eq!(auto(&da, &lfq)?, MaxPeaks::AUTO_LOW_RES);
        // high resolution keeps more peaks than the fixed default anyway
        assert_eq!(auto(&ppm, &tmt_ms2)?, MaxPeaks::AUTO_HIGH_RES);
        assert_eq!(auto(&narrow_da, &tmt_ms2)?, MaxPeaks::AUTO_HIGH_RES);
        // an explicit number is used as given
        assert_eq!(
            resolved_max_peaks_quant(da.clone(), Some(json!(60)), Some(tmt_ms2.clone()))?,
            60
        );
        assert_eq!(
            MaxPeaks::Auto.resolve(&Tolerance::Da(-0.5, 0.5), true, false),
            MaxPeaks::AUTO_LOW_RES_TMT_MS2
        );
        Ok(())
    }

    #[test]
    fn max_peaks_auto_keeps_150_in_wide_window_searches() -> anyhow::Result<()> {
        use serde_json::json;
        let tmt_ms2 = json!({ "tmt": "Tmt16", "tmt_settings": { "level": 2 } });
        for tol in [
            json!({ "ppm": [-20.0, 20.0] }),
            json!({ "da": [-0.5, 0.5] }),
            json!({ "da": [-0.02, 0.02] }),
            json!({ "pct": [-0.002, 0.002] }),
        ] {
            let ww = |extra: serde_json::Value| -> anyhow::Result<usize> {
                let mut extra = extra;
                extra["wide_window"] = json!(true);
                Ok(resolved(tol.clone(), extra)?.max_peaks)
            };
            assert_eq!(
                ww(json!({ "max_peaks": "auto" }))?,
                MaxPeaks::AUTO_WIDE_WINDOW
            );
            assert_eq!(
                ww(json!({ "max_peaks": "auto", "quant": tmt_ms2 }))?,
                MaxPeaks::AUTO_WIDE_WINDOW
            );
            // an explicit number is used as given
            assert_eq!(ww(json!({ "max_peaks": 400 }))?, 400);
            assert_eq!(ww(json!({ "max_peaks": 80 }))?, 80);
            // wide_window false: the resolution rule
            let narrow = resolved(
                tol.clone(),
                json!({ "max_peaks": "auto", "wide_window": false }),
            )?;
            assert_ne!(narrow.max_peaks, MaxPeaks::AUTO_WIDE_WINDOW);
            assert!(narrow.max_peaks == 400 || narrow.max_peaks == 80);
        }
        assert_eq!(MaxPeaks::AUTO_WIDE_WINDOW, 150);
        for tol in [Tolerance::Ppm(-10.0, 10.0), Tolerance::Da(-0.5, 0.5)] {
            for tmt_ms2 in [false, true] {
                assert_eq!(
                    MaxPeaks::Auto.resolve(&tol, tmt_ms2, true),
                    MaxPeaks::AUTO_WIDE_WINDOW
                );
                assert_eq!(MaxPeaks::Count(7).resolve(&tol, tmt_ms2, true), 7);
            }
        }
        Ok(())
    }

    #[test]
    fn max_peaks_auto_keeps_at_least_min_peaks() -> anyhow::Result<()> {
        use sage_core::spectrum::{Precursor, RawSpectrum, Representation, SpectrumProcessor};
        use serde_json::json;
        let da = json!({ "da": [-0.5, 0.5] });
        let ppm = json!({ "ppm": [-20.0, 20.0] });

        // "auto" -> 80 for 0.5 Da is raised to min_peaks 100, so spectra stay searchable
        let search = resolved(da.clone(), json!({ "max_peaks": "auto", "min_peaks": 100 }))?;
        assert_eq!((search.max_peaks, search.min_peaks), (100, 100));
        // at or below the resolved value: unchanged
        let search = resolved(da.clone(), json!({ "max_peaks": "auto", "min_peaks": 80 }))?;
        assert_eq!((search.max_peaks, search.min_peaks), (80, 80));
        let search = resolved(
            ppm.clone(),
            json!({ "max_peaks": "auto", "min_peaks": 100 }),
        )?;
        assert_eq!((search.max_peaks, search.min_peaks), (400, 100));
        // wide window: 150, raised to a higher min_peaks
        let search = resolved(
            ppm.clone(),
            json!({ "max_peaks": "auto", "wide_window": true, "min_peaks": 200 }),
        )?;
        assert_eq!((search.max_peaks, search.min_peaks), (200, 200));
        // an explicit number or the default is used as given (a warning is logged)
        let search = resolved(da.clone(), json!({ "max_peaks": 50, "min_peaks": 100 }))?;
        assert_eq!((search.max_peaks, search.min_peaks), (50, 100));
        let search = resolved(da.clone(), json!({ "min_peaks": 200 }))?;
        assert_eq!((search.max_peaks, search.min_peaks), (150, 200));

        // end to end: a 300-peak MS2 spectrum keeps min_peaks peaks with "auto"
        let search = resolved(da, json!({ "max_peaks": "auto", "min_peaks": 100 }))?;
        let mz: Vec<f32> = (0..300).map(|i| 200.0 + 4.0 * i as f32).collect();
        let intensity: Vec<f32> = (0..300).map(|i| 1.0 + i as f32).collect();
        let spectrum = RawSpectrum {
            ms_level: 2,
            id: "scan=1".into(),
            precursors: vec![Precursor {
                mz: 800.0,
                charge: Some(2),
                ..Default::default()
            }],
            representation: Representation::Centroid,
            mz,
            intensity,
            ..Default::default()
        };
        let processed = SpectrumProcessor::new(search.max_peaks, false, 0.0).process(spectrum);
        assert_eq!(processed.masses.len(), 100);
        assert!(processed.masses.len() >= search.min_peaks);
        Ok(())
    }

    #[test]
    fn max_peaks_rejects_other_values() {
        use serde_json::json;
        for bad in [
            json!(-1),
            json!(1.5),
            json!("Auto"),
            json!("150"),
            json!(true),
            json!([150]),
        ] {
            let err = serde_json::from_value::<MaxPeaks>(bad.clone())
                .expect_err(&format!("{bad} must be rejected"));
            assert!(
                err.to_string()
                    .contains("a non-negative integer or \"auto\""),
                "{bad}: {err}"
            );
        }
        assert_eq!(
            serde_json::from_value::<MaxPeaks>(json!(400)).unwrap(),
            MaxPeaks::Count(400)
        );
        assert_eq!(
            serde_json::from_value::<MaxPeaks>(json!("auto")).unwrap(),
            MaxPeaks::Auto
        );
    }

    #[test]
    fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
        let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
        }))?;
        let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "P",
        }))?;
        let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "",
        }))?;

        let a: EnzymeParameters = a.into();
        let b: EnzymeParameters = b.into();
        let c: EnzymeParameters = c.into();

        assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
        {
            let mut expected = [false; 26];
            expected[(b'P' - b'A') as usize] = true;
            assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
        }
        assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

        // default trypsin (no `cleave_at`) keeps its proline rule in a partial block
        let d: EnzymeParameters = serde_json::from_value::<EnzymeBuilder>(serde_json::json!({
            "missed_cleavages": 2,
        }))?
        .into();
        let mut trypsin = [false; 26];
        trypsin[(b'P' - b'A') as usize] = true;
        assert_eq!(d.enzyme.map(|e| e.skip_suffix), Some(trypsin));

        Ok(())
    }
}
