# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [v0.15.0-fork.3] - 2026-10-05

Faster, with the same results: with the same settings, every result and PIN file in the benchmark runs below is byte-identical to fork.2's. The defaults are unchanged. `max_peaks` stays 150, and the two new quality options are opt-in.

Data behind the numbers in this section:
- **The 20 files**: the public benchmark of OpenMS issue #10364, 8,000 MS2 spectra from each of 20 public runs (PXD001819 Velos UPS1, ProteoBench HF-X and Astral, PXD007683 Lumos LFQ and TMT, PXD054559 Eclipse TMTpro, PXD053748 timsTOF), plus larger inputs built from these files.
- **Held-out**: another 8,000 spectra from each of the same 20 runs, disjoint from the 20 files (ranks 8,001-16,000 of the benchmark's own selection rule).
- **Identifications**: Percolator 3.09 on Sage's PIN file, 20 seeds, target PSMs at 1% FDR, mean per file.
- **Error control**: the doubled-database entrapment FDP, pooled over three shuffles, and the natural entrapment FDP of the three Velos UPS1 runs.
- **Timings**: AMD EPYC 7763 (64 cores, 128 threads), idle, warm page cache. fork.2 and fork.3 alternate in paired ABBA blocks; the numbers are medians.

### Added
- `max_peaks: "auto"` (opt-in; the default stays 150). It resolves to:
  - 400 peaks per MS2 spectrum for a ppm or pct `fragment_tol`, or for a Da tolerance narrower than 0.1 Da;
  - 80 for a Da tolerance that reaches 0.1 Da on one side (low-resolution fragments), or 150 in that case if TMT reporter ions are quantified at MS2;
  - 150 with `wide_window`;
  - never fewer than `min_peaks`: it is raised to it, and the log says so.

  The resolved number is logged and written to `results.json`. Other strings, negative and fractional numbers are rejected with an error that names the accepted values.
  - Identifications: +2.26% PSMs on the 20 files (Astral +10.8%, Lumos LFQ -0.1%) and +2.38% on held-out (Astral +7.7%, Lumos LFQ +0.4%).
  - Error control: on the 20 files, the doubled-database FDP fell (1.08% -> 1.03%) and the Velos UPS1 FDP rose (1.20% -> 1.25%). On held-out, the Velos UPS1 FDP fell (1.06% -> 0.95%), but the doubled-database FDP rose by 0.070 points (1.00% -> 1.07%). That is above the pre-registered limit of +0.05 points, so `"auto"` did not become the default. The rise is about 1.5 standard errors above zero (0.4-0.6 standard errors above the limit) and is not significant.
  - Wide-window and DIA-style searches: with 400 peaks, wrong candidates from the wide isolation window win more often. Percolator PSMs fell by 3-46% on HF-X, timsTOF and Astral files (Astral -46%). On a public diaPASEF run (PXD017703), they fell by 12-16% at 1.9x the run time. `"auto"` therefore resolves to 150 with `wide_window`.
  - Time: on the 20 files at 4 threads, the total stays about the same (+0.2%). High-resolution files take up to 7.7% longer and 0.5 Da files up to 7.0% less, because building the index takes most of those runs. The search phase itself takes 70-80% longer on Astral spectra (0.46-0.48 -> 0.82-0.84 s per file), so larger high-resolution inputs take longer. On 24,000 Astral spectra, runs took 7-12% longer at 4-64 threads; on 24,000 Velos spectra (80 peaks), they were 9-11% shorter.
  - Memory: the extra peaks are held for every MS2 spectrum of a batch. That is about 0.30 GiB per full-size Astral run (about 120,000 MS2 spectra) in a batch, or 0.50 GiB with `report_psms` 10. By default a batch holds CPUs/2 files (`--batch-size`), usually all files of a run. A smaller `--batch-size` or `max_peaks: 150` bounds it.
  - `"auto"` decides by the unit of `fragment_tol` only. For other cases (e.g. low-resolution spectra searched with a ppm tolerance), set a number.
- `ion_model` (default false; experimental): a self-trained, cross-fitted fragment-ion likelihood model (after ProSE's `FragmentIonLikelihoodModel`). It adds the features `ion_llr` and `ion_explained` to results.sage.tsv, results.sage.parquet and the PIN file.
  - For each file, the model learns how the fragment ions of confidently identified peptides show up in the spectra: all deisotoped peaks, before the `max_peaks` cut, compared against a reversed-sequence noise model.
  - The spectra are split into two folds, and every PSM is scored by the model of the other fold.
  - The features are meant for Percolator or mokapot; Sage's own LDA does not use them. All other output is unchanged.
  - Identifications (`max_peaks` 150): +5.2% PSMs on the 20 files (Astral +24%) and +5.5% on held-out (Astral +25%).
  - Error control was not confirmed on held-out data, so the option is not recommended. The Velos UPS1 FDP rose from 1.06% to 1.30%: all three Velos files and all 20 seeds went up, and the pre-registered limit was 1.16%. The doubled-database FDP rose by only 0.03 points (1.00% -> 1.03%). On the 20 files the two measures went from 1.20% to 1.27% and from 1.08% to 1.09%.
  - Together with `max_peaks: "auto"`: +6.4% PSMs on the 20 files and +6.0% on held-out, but the held-out doubled-database FDP rose by 0.077 points.
  - Cost (measured with `max_peaks: "auto"` in both arms): +2.4% wall time on the 20 files at 4 threads. On 24,000 Astral spectra: +2.5% at 16 threads and +4.5% at 64, with +0.13 GiB peak memory.
  - Memory: the peak lists of all MS2 spectra of a batch stay in memory until the batch is searched, about 0.35 GiB per 72,000 Astral spectra. Twelve such files in one batch peaked at 28.5 instead of 24.4 GiB. A smaller `--batch-size` bounds it.
  - DOCS.md describes the design limits: the match window is centred on the theoretical mass, the folds are assigned per spectrum rather than per peptide, and the training gate applies to the whole file.
- Parallel mzML parse: local, uncompressed mzML files are parsed in parallel chunks of whole spectra, and the spectra are the same as with the serial parse.
  - Gzipped and remote files, and files with an unusual layout (comments, CDATA or processing instructions between the spectra), are parsed serially.
  - Setting the environment variable `SAGE_MZML_SERIAL` to any value except `0` (e.g. `SAGE_MZML_SERIAL=1`) forces the serial parse. When a local, uncompressed mzML file falls back to the serial parse, the reason is logged at debug level (`SAGE_LOG=sage_cloudpath=debug`).
  - Local mzML, MGF and FASTA files are read in 2 MiB reads instead of object_store's 8 KiB stream. A local file that cannot be opened is named in the error.
  - Speed, measured on this change alone against fork.2 on 72,000 Astral spectra: 1.42x faster at 16 threads (13.9 -> 9.8 s), 1.63x at 64 (11.7 -> 7.2 s) and 1.61x at 128.
  - In this release, the parallel parse is neutral on the 8,000-spectrum files at 64-128 threads (20-file sums 1.007x and 1.001x against the same binary with `SAGE_MZML_SERIAL=1`). It is 1-6% faster when one call reads three such files. The small costs measured on this change alone in those cases (0.5-2.8%) do not occur in the release.
  - Memory is bounded by the chunks being parsed, not by the whole file's raw spectra. This also holds when a long tail follows the spectra (chromatograms, index).
  - Each file is read twice: a scan for spectrum boundaries, then the chunks. The second read comes from the page cache, unless the files of a batch do not fit in free memory; then the device reads them twice.
  - On storage slower than about 300 MB/s (HDD, 1 GbE network shares), that can make the parallel parse slower than the serial one. Use `SAGE_MZML_SERIAL=1` or a smaller `--batch-size` there. On local NVMe and on Ceph, the parallel parse stayed 3.3-4.5x faster than the serial one even when the file was read twice (cold page cache).

### Changed
- Performance, with identical results (`max_peaks` 150 in both builds):
  - the 20 files at 4 threads: 220.8 -> 131.0 s (1.69x; every file 1.59-1.78x). CPU time 720 -> 470 s (-35%); peak memory 0.7-1.0 GiB lower per file (3.3-4.5 -> 2.6-3.7 GiB);
  - by thread count, on 4 of the files (Velos, HF-X, Astral, timsTOF):
    - 1.61x at 1 thread (129.6 -> 80.7 s);
    - 2.10x at 16 (21.9 -> 10.4 s);
    - 2.45x at 64 (18.5 -> 7.5 s);
    - 2.29x at 128 (18.6 -> 8.1 s). On files this small, 128 threads are 8% slower than 64.
  - larger inputs:
    - 24,000 spectra (three Astral or three Velos files concatenated): 1.44x, 1.79x and 2.01x at 4, 16 and 64 threads;
    - 72,000 Astral spectra (the three Astral files, three times over): 1.29x at 4 threads, 2.0-2.1x at 16 and 2.8-3.2x at 64-128;
    - three files in one call (default batch): 1.83x at 16 threads and 2.02x at 64;
    - a gzipped 24,000-spectrum Astral file: 1.38x, 1.32x and 1.19x at 8, 16 and 64 threads.
  - Where the gain comes from:
    - The fragment index is built by m/z bins with counting sorts, instead of one global sort. Its pages are radix-sorted by peptide.
    - Only the fragments of peptides that a precursor window can reach are indexed (see the next item).
    - Digest groups are built in parallel. The fork.2 notes already listed this, but that change had been reverted before fork.2 was released.
    - The peptide sort goes through 16-byte keys, and duplicates are merged in parallel.
    - Freed memory is returned to the OS before the peptide sort and before the fragment index is allocated.
    - Local mzML files are parsed in parallel (see Added), and the end of a run is faster (see below).
- Fragment index pruning applies when every input file is read in the first batch (at most `--batch-size` files, by default CPUs/2) and `database.prefilter` is off.
  - The index then holds only the fragments of peptides that a precursor window of those spectra can reach. On the 20 files, that drops 35-64% of the fragments. The peptide list stays complete, and the results are identical.
  - The build does not wait for the read. It prunes at the first of three points (before counting the fragments, after counting, before bucketing) at which all spectra have been read.
  - If the read is still going after the last point, the full index is built. This is typical for gzipped, remote or otherwise serially parsed input at high thread counts. Multi-batch runs and the prefilter also build the full index, as before.
  - The telemetry field `fragments` now reports the indexed (pruned) count.
- Peak memory is 0.6-2.0 GiB below fork.2 in every default-mode benchmark run of plain mzML files (e.g. 72,000 Astral spectra: 5.6-6.0 -> 4.5-5.4 GiB). Exceptions:
  - Gzipped (or otherwise serially parsed) input whose read ends after the fragments have been counted: the fragment array is allocated at full size (often the full index), and peak memory is level with fork.2 or slightly above it. That is +0.03 to +0.11 GiB for a gzipped 24,000-spectrum Astral file at 8-64 threads, and up to +0.5 GiB for three full gzipped Lumos LFQ runs (with MS1 and LFQ) in one batch at 32 threads.
  - Several gzipped (or otherwise serially parsed) files in one batch at 64-128 threads: release candidates occasionally peaked 1.0-1.2 GiB above fork.2 (5 of 59 runs that took a late prune or the full index). Freed memory is now also returned right before the fragment index is allocated. After that change, no such peak occurred in 148 benchmark runs (13 of 212 without the change), but one of 22 profiled runs still reached about 1 GiB above fork.2. Such peaks are now rare, but they still occur.
  - `max_peaks: "auto"` and `ion_model` hold more memory per MS2 spectrum of a batch (see Added).
- Prefilter mode (`database.prefilter: true`) is 1.8x faster than fork.2 at 64 threads, with 0.24-0.42 GiB lower peak memory (velos_125_R1 and astral_B1 of the 20 files). Each FASTA chunk's index is now returned to the OS before the next chunk is digested; without that step, the faster build peaked 0.09-0.13 GiB above fork.2.
- The end of a run is faster (about 0.1-0.3 s per run, more with more PSMs; output unchanged):
  - the system query for telemetry runs only when the report is sent;
  - the `sage` binary leaves the search results to the OS at exit;
  - protein names for protein-level FDR are formatted once per protein.
- The crates are versioned 0.15.0-fork.3, so `sage --version`, `results.json` and the telemetry record name the release. fork.2's binaries reported 0.15.0-beta.2.
- The declared minimum Rust versions now match a `--locked` build: sage-core 1.80, sage-cloudpath 1.85 and sage-cli 1.88 (also needed for `--features mzpeak`). They were 1.62, which the code had outgrown.
- Library API (sage-core, sage-cloudpath, sage-cli):
  - `sage_cli::input::Input::max_peaks` is an `Option<MaxPeaks>` (`MaxPeaks::Auto` or `MaxPeaks::Count(n)`) instead of an `Option<usize>`. `MaxPeaks::default()` is `Count(150)` (`MaxPeaks::DEFAULT`). `MaxPeaks::resolve(fragment_tol, tmt_ms2, wide_window)` returns the number of peaks.
  - New public fields:
    - `ProcessedSpectrum::ion_evidence` (`Option<Box<IonEvidence>>`);
    - `SpectrumProcessor::ion_evidence` (set with `with_ion_evidence`);
    - `Feature::{normalized_hyperscore, ion_llr, ion_explained}` (the JSON of a `Feature` gains `ion_llr` and `ion_explained`);
    - `Search::ion_model` and `Input::ion_model`.
  - New module `sage_core::ion_model`. `sage_cloudpath::parquet::build_schema` and `serialize_features` take an `ion_model` flag.
  - `Runner::new` may return a `database` whose fragment index covers only the peptides reachable from the first batch's spectra, under the parameters given to `new`. Call `run` with the same `parallel` and unchanged `parameters`, and do not use `runner.database` to score other spectra or to search with other tolerances.
  - `Runner::run` frees the database and the results. `Runner::run_then_exit` leaves them to the OS, for a process that exits right after it (the `sage` binary).
  - `sage_core::database::set_release_freed_memory` registers a function that the database build calls before its large allocations; `sage_core::database::release_freed_memory` calls it. sage-cli registers mimalloc's `mi_collect` (`sage_cli::release_freed_memory`).
  - `Telemetry` fills in the OS name and the total memory only in `send()`.

### Fixed
- An explicit `max_peaks` below `min_peaks` now logs a warning. With such a setting no MS2 spectrum can be searched, and the run ended with empty results and exit code 0 without saying why.
- `Runner::run` (library use) no longer leaks the database. fork.2 kept the database of every call in memory.

## [v0.15.0-fork.2] - 2026-10-01

First release of the fork, based on upstream master 2c9922e (crate version 0.15.0-beta.2).

### Added
- mzPeak input (`.mzpeak`, optional cargo feature `mzpeak`): read through the HUPO-PSI reference reader (`mzpeak_prototyping`) and `mzdata`. The reader is pinned to the public fork okohlbacher/mzPeak at bfec73a (branch `optional-bruker`), which makes its `bruker` feature optional (not upstream yet). Results are identical to searching the same data as mzML. Grid-encoded timsTOF mzPeak files are not supported in this configuration.
- From upstream 2c9922e: per-modification occurrence limits using `{"mass": <mass>, "max_count": <limit>}` entries in `database.variable_mods`; existing bare-mass entries remain supported.
- From upstream 2c9922e: `database.max_combinations` to cap the number of peptide variants (including the unmodified form) generated from variable modifications, preferring variants with fewer modifications.

### Changed
- Performance: about 4.6-5.8x faster end to end (128/64/16 threads: 35.1/41.7/88.7 s -> 7.6/7.2/15.3 s) and 2.2-3.6 GiB lower peak memory on a 3-file timsTOF DDA benchmark (human, 7.7 M peptides), with identical identifications. Isotope windows share one fragment lookup, a per-page skip table narrows fragment lookups, lookups are pipelined with software prefetch, full rescoring walks the peaks with a cursor instead of a binary search per ion, KDE fits parallelise over bins, the fragment index is built in one exact-size allocation on transparent huge pages, spectra are processed per file, the database build overlaps with reading the first batch, the database is not freed at exit, output is serialized in parallel, and mimalloc is the global allocator. Bruker DDA: each MS2 frame is decoded once per block of spectra (timsrust decoded it once per precursor, ~9 times), frames are read without a memory map, spectra are processed while the file is read, and up to 4 `.d` files are read at a time.
- mimalloc trades memory for speed at high thread counts; `MIMALLOC_ARENA_EAGER_COMMIT=0` lowers peak memory (by ~1.5 GB at 128 threads in the benchmark) for ~10% more time.
- Output is deterministic: identical searches give byte-identical result and PIN files, also across thread counts. `psm_id` is the row number of the output (it was a counter shared by the parallel search).
- mzML: MS1 spectra are no longer decoded or kept unless LFQ is enabled.
- A file that cannot be read completely (I/O error, truncated or corrupt mzML/MGF, unreadable Bruker frames) now fails the run with an error naming the file, instead of being searched partially or skipped with exit code 0. Unsupported file formats are rejected before the database is built.
- `posterior_error` now reports the PEP of a target PSM (min(1, p/(1-p))); it used to report the probability of being a decoy, which is ~2x too optimistic for low-scoring PSMs. It is 0 (PEP = 1) when the LDA falls back to the heuristic score. q-values are unchanged.
- The PIN file no longer contains `posterior_error` (it is derived from Sage's own label-trained model and leaked target/decoy labels into Percolator/mokapot); `sqrt(delta_mobility)` now holds the square root, as its header says.
- An `enzyme` block without `cleave_at` keeps trypsin's `restrict: "P"` default (e.g. `{"missed_cleavages": 2}` used to cleave before P).
- MGF spectra listing several charges (`CHARGE=2+ and 3+`) are searched at every listed charge (only the first was searched).
- The logged "target peptide-spectrum matches" count no longer includes passing decoys.
- FASTA sequences are upper-cased and a terminal `*` is removed; empty headers get `unnamed_protein_<n>` (renamed if a record is explicitly called that).
- MGF: a nested `BEGIN IONS`, a stray `END IONS` or a corrupt intensity fails the file; lines between blocks are ignored (they used to be parsed into the next spectrum); file-level `CHARGE`/`TOL` apply to the first spectrum too; a leading byte-order mark is ignored; a file with content but no block warns.
- Deisotoping uses the highest precursor charge of a spectrum (unknown counting as 3), independent of the order of listed charges.
- The ion mobility model ignores PSMs without a measured mobility (e.g. mzML next to `.d` files) and gives them the median residual.
- Output file names are percent-decoded (`my file.mzML`, not `my%20file.mzML`).

### Fixed
- Crashes of the whole run: `*`/lower-case/non-letters after a cleavage site, FASTA or prefilter chunks without peptides, profile MS2 spectra (now skipped with a warning), m/z and intensity arrays of different lengths, MS2 spectra without precursor m/z, batch size 0 on single-CPU machines, NaN values in the HTML report.
- mzML: `referenceableParamGroupRef` is resolved (#232); MS-Numpress arrays are rejected instead of decoded as garbage; unknown cvParams no longer drop arrays; per-spectrum state no longer leaks into the next spectrum; multi-member gzip files are read completely; `.GZ` is recognised.
- MGF: multi-digit charges (`CHARGE=10+`), empty files, `PEPMASS` of 0.
- Low-memory prefilter kept the highest peptide indices instead of the best-scoring candidates; decoys colliding with targets across prefilter chunks corrupted protein lists.
- `isotope_errors: [n, n]` with n != 0 searched isotope 0.
- Files without retention times (MGF without RTINSECONDS) disabled LDA rescoring for the whole run and distorted RT alignment and the RT model of the other files.
- Bruker `ion_injection_time` reported the retention time.
- mzML: scan-level ion mobility applies to every precursor of a spectrum, and isolation windows no longer carry over from one precursor to the next.
- HTML report: no abort when a discriminant score is not finite.
- Non-finite peaks are dropped from spectra with an ion mobility array too.

## [v0.15.0]
### Added
- IDPicker-based protein grouping with picked group FDR control (`protein_grouping` setting, enabled by default). Proteins are grouped using a bipartite graph greedy set cover approach, and protein group-level q-values are reported via target-decoy competition. New output columns: `protein_groups`, `num_protein_groups`, `protein_group_q`.
- `protein_grouping_peptide_fdr` parameter to control the peptide FDR threshold used for confident peptides during protein grouping (default: 0.01)
- Initial support for LFQ on data with ion mobility.
- Speedup on the generation of databases when large number of peptides are redundant.
- Initial support for searching diaPASEF data
- `override_precursor_charge` setting that forces multiple charge states to be searched
- Cross-cloud storage support: Replaced AWS SDK with the `object_store` crate. Sage now natively supports reading/writing from Amazon S3, Google Cloud Storage, and Azure Blob Storage using `s3://`, `gs://`, and `az://` URL schemes.
- HTML QC report generation (`--write-report`)
- `lfq_settings.peptide_q_value` parameter for controlling which peptides are quantified
- Stack size configuration for Rayon threads
- Allow zero or multiple amino-acids as cleavage restrictions
- CITATION.cff file

### Fixed
- Handle negative mass errors for .pin files
- Extract precursor m/z from 'isolation window target m/z' if missing
- Selected ion m/z of 0.0 was overwriting precursor.mz
- C/N-term mixup in modification handling
- Bruker `.tdf` filename handling (use parent directory name)
- Performance optimizations on prefiltering
- Picked protein FDR now correctly uses only proteotypic peptides for competition

### Breaking Changes
- `precursor_ppm` field reports the non-absoluted average mass error, rather than the absoluted average mass error.
- Don't deisotope reporter ion regions if MS2-based TMT/iTRAQ is used
- Removed `fragment_min_mz` and `fragment_max_mz` parameters. These were decreasing the accuracy of preliminary scoring estimation when attempting to annotate multiply-charged, high-m/z ions.
- `sage-cloudpath` no longer exposes the `CloudPath` type. All paths are represented as URLs (`url::Url`). Local paths are converted to `file://` URLs internally.

## [v0.14.7]
### Added
- Added columns missing from parquet output: `semi_enzymatic` and `missed_cleavages`
### Changed
- Fixed ion mobility parsing from some mzMLs
- MGF paths were being lowercased prior to parsing

## [v0.14.6]
### Added
- Support for MGF files
- Support for writing ion mobility measurements to output files: `ion_mobility`, `predicted_mobility`, `delta_mobility` added to primary tsv and parquet reports. Ion mobility is predicted in a similar manner to RT, using a linear model trained on the data from the search.

## [v0.14.5]
### Added
- Support for semi-enzymatic digests (`database.enzyme.semi_enzymatic` parameter)
- Ability to directly export matched fragment ions (e.g. for spectral library or rescoring) with the `--annotate-matches` CLI option. This is compatible with the `--parquet` CLI option as well. Annotations will be written to `matched_fragments.sage.tsv` or `matched_fragments.sage.parquet`
- Sage sends basic telemetry data (version of Sage, run time, OS, # of CPU cores, # of peptides in database, whether LFQ is used) to a remote server. No information about your actual data is sent - e.g. identifications, quantities, organism, or modifications are NOT tracked or reported.  This data will be used to help focus efforts on improving Sage and figuring which features are most used. Please take a look at `crates/sage-cli/src/telemetry.rs` to see exactly what is sent! You can disable sending telemetry data  by using the `--disable-telemetry-i-dont-want-to-improve-sage` CLI flag.
### Changed
- Modified visibility on some crate internals to support the [sagepy project](https://github.com/theGreatHerrLebert/sagepy)
- Added `psm_id` field to various output files to match the new `--annotate-matches` option.
### Removed
- Removed the `ms1_intensity` field from CSV output, since it is essentially useless


## [v0.14.4]
### Added
- **Unstable feature**: Preliminary support for reading Bruker .d folders (ddaPASEF; no MS1/LFQ support yet)
### Changed
- Retention times are converted to minutes
### Fixed
- Fixed bug where charge state 1 would never be searched

## [v0.14.3]
### Fixed
- Hotfix for bug in parquet LFQ writer

## [v0.14.2]
### Added
- `quant.lfq_settings.combine_charge_state` boolean option. By default this is set to `true`, and LFQ is performed on the peptide-level, where all charge states are treated as the same precursor. Setting this to `false` performs LFQ on the peptide-charge-level, where each charge state will be treated separately.
### Changed
- Percolator output format now contains the integer-valued charge state encoded in the `z=other` column, if the charge state is outside the range 2-6 (e.g. a value of 7 will appear in the `z=other` column, rather than it being one-hot encoded)
- LFQ uses the the charge state range from the `precursor_charge` configuration option for tracing MS1 peaks

## [v0.14.1]
### Added
- Added additional output showing search progress if `SAGE_LOG=trace` environment variable is set
- Added additional warnings about precursor tolerances
- Added configuration option `precursor_charge` to make it explicit what charge states are being searched in the case where the mzML does not contain charge state information, or where `wide_window` is turned on.
### Changed
- Added a warning message if variable modifications are specified as single values (e.g. `15.9949`) instead of lists of values (e.g. `[15.9949]`). By v0.15 this will become a hard error and will not parse, to simply some of the internal logic.

## [v0.14.0]
### Added
- Support for parquet file format output. Search results and reporter ion quantification will be written to one file (`results.sage.parquet`) and label-free quant will be written to another (`lfq.parquet`). Parquet files tend to be significantly smaller than TSV files, faster to parse, and are compatible with a variety of distributed SQL engines.
### Changed
- Implement heapselect algorithm for faster sorting of candidate matches (#80). This is a backwards-incompatible change with respect to output - small changes in PSM ranks will be present between v0.13.4 and v0.14.0

## [v0.13.4]
### Fixed
- Bug in mzML parser, where some older specification-compliant mzMLs would not parse. If your mzMLs previously parsed, then there will be no change in behavior. Added a test case

## [v0.13.3]
### Fixed
- Bug in `database.enzyme.restrict` parameter, where `null` values were being overriden with "P" (causing Trypsin/P to behave like Trypsin)

## [v0.13.2]
### Changed
- Subtle change to TMT integration tolerance, and selection of which ion to quantify (most intense). As a result, TMT integration should be more in agreement (if not 100% so) with ProteomeDiscover/FragPipe/etc
- Remove `delta_mass` (precursor ppm) LDA feature - instead, build a delta mass (or ppm) profile using KDE/posterior error calculation code, and use the P(decoy) as a feature for LDA.

## [v0.13.1]
### Changed
- Internal performance and stability improvements for RT prediction & LDA

## [v0.13.0]
### Added
- Better error reporting thanks to @Elendol
- Added support for multiple variable mods for the same amino acid
- Added support for N/C-terminal modifications specific to an individual amino acid

New syntax:
```json
"variable_mods": {
    "M": [15.9949],
    "^Q": -17.026549,
    "^E": -18.010565,
    "[": 42.010565
}
```

Either a single floating point number (-18.0) or a list of floating point numbers ([-18.0, -15.2]) can be supplied as modifications. Support for single values may eventually be phased out to simplify the parser.

### Changed
- Changed "_fdr" columns to "_q" (e.g. "spectrum_q") in "results.sage.tsv" file
- Changed internal data representation of `Peptide` struct to allow for sharing of sequences (using `Arc`) among modified peptides
- Fragment index creation should now be faster

## [0.12.0]
### Added
- Add `wide_window` option to configuration file. This option turns off `precursor_tol`, instead using the isolation window written in the mzML file.
### Changed
- Changed internal calculation of precursor tolerances when searching with `isotope_errors`. The new version should be more accurate. This change also enables a significant boost to search speed for open searches.

## [0.11.2]
### Added
- Add rank & charge features to LDA
### Changed
- One-hot encode charge state information for percolator `.pin` files
- Change PSMId -> SpecId for Mokapot compatibility with `.pin` files

## [0.11.1]
### Added
- Support for additional fragment ion types, via the "database.ion_kinds" configuration option. Valid values are "a", "b", "c", "x", "y", "z"
### Changed
- Sort protein names alphanumerically for each peptide entry. This should enhance stability across runs, and fixes a bug with picked-protein group FDR
- Fix another bug where picked-FDR approaches assume internal decoy generation

### Changed
- Modify order of operations during deisotoping. Deisotoped peaks can contribute intensity to only 1 parent peak now, rather than potentially multiple parent peaks

## [0.11.0]
### Added
- Support for percolator output files (`--write-pin` CLI flag)
- Support for modifying file batch size (`--batch-size N` CLI flag)
- Add `delta_best` feature, which reports the delta hyperscore from the best match to current ranked PSM
- Add Sage version to `results.json` files

### Changed
- Breaking changes to `quant` section of the configuration file format
- Rename `delta_hyperscore` to `delta_next`
- Altered internal scoring algorithm. Rather than consider all MS2 peaks within a m/z tolerance window to be matches to a theoretical spectrum, consider only the closest peak. This should increase the accuracy of # of matched peaks, and subsequent scores
- Overhaul of chimeric scoring, `report_psms` can now be used to search for multiple chimeric spectra
- Completely overhauled the LFQ algorithm: added match-between runs, peak scoring using normalized spectral angle relative to theoretical isotopic envelope, target decoy scoring of MS1 integration
- Fixed bug in picked-peptide FDR that could lead to liberal FDR
- Fixed bug in picked-protein FDR that could lead to conservative FDR
- Fixed bug where using variable protein terminal (e.g. protein N-terminal acetylation) modifications could cause some determinism. This also improves the accuracy of peptide => protein assignment. Unfortunately this fix has performance implications, causing creation of the fragment index to take up to ~2x as long.

### Removed
- Remove `no-parallel` CLI flag, and `parallel` configuration file entry

## [0.10.0]
### Added
- Retention times are now globally aligned across files
- RT prediction is then performed on all files at once (on aligned RTs), rather than one file at a time - previously, there were many instances where some files in a search could not have RTs predicted, decreasing the effectiveness of delta_rt as a feature for LDA.

### Changed
- Peptide sequences within a protein are now deduplicated - previously, repeated peptides would be called multiple times for the same protein (e.g. num_proteins > 1 even if the peptide was unique)

## [0.9.4]
### Changed
- Fix issues with RT prediction (and occasionally LDA) that arise from 0's being present on the diagonals of the covariance matrix (small amount of regularization added)

## [0.9.3]
### Added
- Allow users to set minimum number of matched b+y ions for reporting PSMs (`min_matched_peaks`)

### Changed
- Internal code for calculating factorials

## [0.9.2]
### Added
- Added option for TMT signal/noise quantification, if noise values are present in mzML

## [0.9.1]
### Changed
- FASTA file path, JSON configuration file can now be specified as "s3://" paths, allowing Sage to run completely disk-free

## [0.9.0]
### Added
- Support for non-specific digests, N-terminal enzymatic digestion

## [0.8.1]
### Added
- `quant.tmt_level` configuration option to enable MS2 (or MSn) isobaric quantification

## [0.8.0]
### Added
- Support for protein N-terminal ('['), C-terminal (']') as well as peptide C-terminal ('$') modifications
- Support for k-combinations of variable modifications. This can be specified with the `database.max_variable_mods` parameter

## [0.7.1] - 2022-11-04
### Changed
- Fix bug with in silico digest: Logic around overwriting decoys with target sequences was incorrect peptides shared between targets/decoys were being annotated as decoy peptides but assigned to non-decoy proteins. We now make sure that they are assigned to non-decoy proteins and also annotated as target sequences.

## [0.7.0] - 2022-11-03
### Added
- Add support for user-specified enzymes to JSON file.  `database.enzyme.sites` and `database.enzyme.restrict` are limited to valid amino acids
- Sage can now search MS2 spectra without annotated precursor charge states. Default behavior is to search with z=2, z=3, z=4, and then merge the PSMs for scoring

### Changed
- Configuration file schema changed. `peptide_min_len`, `peptide_max_len`, `missed_cleavages` are now specified under `database.enzyme` in the JSON file
- Internal behavior of Sage was changed to enable deterministic searching
- Docker file changed from Alpine to Debian


## [0.6.0] - 2022-11-01
### Added
- Changelog
- `rank` column added to output file
- `database.generate_decoys` parameter, which turns off internal decoy generation. This enables the use of FASTA databases for SearchGUI/PeptideShaker

### Changed
- Base ProForma v2 notation is used for peptide modifications, i.e. "\[+304.2071\]-PEPTIDEM\[+15.9949\]AAC\[+57.0214\]H"
- `scannr` column now contains the full nativeID/spectrum title from the mzML file, i.e. "controllerType=0 controllerNumber=1 scan=30069"
- `discriminant_score` column renamed to `sage_discriminant_score` for PeptideShaker recognition
- `database.decoy_prefix` JSON option changed to `database.decoy_tag`. This allows decoy tagging to occur anywhere within the accession: "sp|P01234_REVERSED|HUMAN"
- Output file renamed:  `results.pin` to `results.sage.tsv`
- Output file renamed: `quant.csv` to `quant.tsv`
- Rename `pin_paths` to `output_paths` in results.json file


## [0.5.1] - 2022-10-31
### Added
- Support for selenocysteine and pyrrolysine amino acids

## [0.5.0] - 2022-10-28
### Added
- Ability to directly read/write files from AWS S3

### Changed
- Processing files in parallel processes them in batches of `num_cpus / 2` to avoid memory issues
- Fixed bug where `protein_fdr` was erroneously assigned to `peptide_fdr` output field
- Additional parallelization for assignment of PEP, FDR, writing output files

## [0.4.0] - 2022-10-18
### Added
- Label free quantification can be enabled by turning on `quant.lfq` JSON parameter 
- Commmand line arguments can be used to override configuration file

## [0.3.1] - 2022-10-06
### Added
- Workflow contributions from [@wfondrie](https://github.com/wfondrie).

### Changed
- Don't parse empty MS2 spectra

## [0.3.0] - 2015-09-15
### Added
- Retention time prediction
- Ability to filter low-number b/y-ions for faster preliminary scoring (`database.min_ion_index` option)
- Ability to toggle retention time prediction (`predict_rt`)
