//! Send a minimal telemetry report

use sage_core::tmt::Isobaric;
use serde::Serialize;
use sysinfo::{System, SystemExt};

#[derive(Debug, Serialize)]
pub struct Telemetry {
    // Which version of Sage?
    version: String,
    // How many peptides are in the fragment index?
    peptides: usize,
    // How many fragments are in the index? When every input file is read in the first
    // batch, the index holds only the fragments of the peptides that a precursor window
    // of the searched spectra can reach (see `Runner::new`), so this can be well below
    // the number of fragments of all peptides.
    fragments: usize,
    // How many files are being processed?
    files: usize,
    // How long did analysis take?
    runtime_secs: u64,

    // Is LFQ used?
    lfq: bool,
    // Which kind of TMT tags are used, if any?
    tmt: Option<Isobaric>,
    // Are results written in parquet format?
    parquet: bool,

    // Details about the operating system and computer:
    // - Which OS?
    // - Total memory available
    // - Number of CPU cores
    os_name: String,
    total_memory: u64,
    cpus: usize,
}

impl Telemetry {
    pub fn new(
        settings: crate::input::Search,
        peptides: usize,
        fragments: usize,
        parquet: bool,
        runtime_secs: u64,
    ) -> Telemetry {
        // The system details are filled in by `send`: querying them here (it was a full
        // `sysinfo` refresh of every process, disk, network and sensor) took 0.1-0.2 s at
        // the end of every run, also with telemetry disabled.
        Telemetry {
            version: settings.version,
            peptides,
            fragments,
            files: settings.mzml_paths.len(),
            runtime_secs,
            lfq: settings.quant.lfq,
            tmt: settings.quant.tmt,
            parquet,
            os_name: String::new(),
            total_memory: 0,
            cpus: num_cpus::get(),
        }
    }

    pub fn send(mut self) {
        // only memory needs a refresh; the OS version is read directly
        let mut system = System::new();
        system.refresh_memory();
        self.os_name = system.long_os_version().unwrap_or_default();
        self.total_memory = system.total_memory();

        log::trace!("sending telemetry...");
        // doesn't matter if it fails
        match sage_cloudpath::util::send_data(
            "https://pax3h44gubc6o5ci23knddnw2i0qnuaz.lambda-url.us-west-2.on.aws/",
            &self,
        ) {
            Ok(_) => {
                log::trace!("telemetry data sent successfully!")
            }
            Err(e) => {
                log::trace!("error while sending telemetry: {}", e)
            }
        }
    }
}
