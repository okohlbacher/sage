# Sage Documentation

The most up-to-date documentation now lives [here](https://sage-docs.vercel.app/docs).


## Features & Information

### Assign multiple peptides to complex spectra

<img src="figures/chimera_27525.png" width="800">

- When chimeric searching is enabled, multiple peptide identifications can be reported for each MS2 scan

### Sage trains machine learning models for FDR refinement and posterior error probability calculation

- Retention times are globally aligned across runs
- Boosts PSM identifications using prediction of retention times with a [linear regression](https://doi.org/10.1021/ac070262k) model
- Hand-rolled, 100% pure Rust implementations of Linear Discriminant Analysis and KDE-mixture models for refinement of false discovery rates
- Models demonstrate 1:1 results with scikit-learn, but have increased performance
- No need for a second post-search pipeline step

<img src="figures/SageLDA.png" width="600px">

## Installation

Sage is distributed as source code, and as a standalone executable file.

### Installing via conda

Sage can be installed from [bioconda](https://anaconda.org/bioconda/sage-proteomics):

```
$ conda install -c bioconda -c conda-forge sage-proteomics
$ sage --help
```

### Compiling the development version

1. Install the [Rust programming language compiler](https://rustup.rs/)
2. Download Sage source code via git: `git clone https://github.com/lazear/sage.git` or by [zip file](https://github.com/lazear/sage/archive/refs/heads/master.zip)
3. Compile: `cargo build --release`
4. Run: `./target/release/sage config.json`

Once you have Rust installed, you can copy and paste the following lines into your terminal to complete the above instructions, and run Sage on the example mzML provided in the repository (a single scan from PXD016766)

```sh
git clone https://github.com/lazear/sage.git
cd sage
cargo run --release tests/config.json 
```

### Downloading the latest release

1. Visit the [Releases](https://github.com/lazear/sage/releases/latest) website.
2. Download the correct pre-compiled binary for your operating system.
3. Run: `sage <path/to/config.json>`

### Interfacing with AWS S3

Sage is capable of natively reading & writing files to AWS S3:

- S3 paths should be specified as `s3://bucket/prefix/key.mzML.gz` or `s3://bucket/prefix` for output folder
- See [AWS docs](https://docs.aws.amazon.com/sdk-for-rust/latest/dg/credentials.html) for configuring your credentials
- Using S3 may incur data transfer charges as well as multi-part upload request charges.

## Usage 

```shell
Usage: sage [OPTIONS] <parameters> [mzml_paths]...

🔮 Sage 🧙 - Proteomics searching so fast it feels like magic!

Arguments:
  <parameters>     Path to configuration parameters (JSON file)
  [mzml_paths]...  Paths to mzML files to process. Overrides mzML files listed in the configuration file.

Options:
  -f, --fasta <fasta>
          Path to FASTA database. Overrides the FASTA file specified in the configuration file.
  -o, --output_directory <output_directory>
          Path where search and quant results will be written. Overrides the directory specified in the configuration file.
      --batch-size <batch-size>
          Number of files to search in parallel (default = number of CPUs/2)
      --parquet
          Write parquet files instead of tab-separated files
      --write-pin
          Write percolator-compatible `.pin` output files
  -h, --help
          Print help information
  -V, --version
          Print version information
```

Sage is called from the command line using and requires a path to a JSON-encoded parameter file as an argument (see below). 

Example usage: `sage config.json`

Some options in the parameters file can be over-written using the command line interface. These are:

1. The paths to the mzML data
2. The path to the database (fasta file)
3. The output directory

For example: 

```
# Specify fasta and output dir:
sage -f proteins.fasta -o output_directory config.json

# Specify mzML files:
sage -f proteins.fasta config.json *.mzML

# Specify mzML file located in an S3 bucket
sage config.json s3://my-bucket/YYYY-MM-DD_expt_A_fraction_1.mzML.gz
```

Running Sage will produce several output files (located in either the current directory, or `output_directory` if that option is specified):
- Record of search parameters (`results.json`) will be created that details input/output paths and all search parameters used for the search
- MS2 search results will be stored as a tab-separated file (`results.sage.tsv`) file - this is a tab-separated file, which can be opened in Excel/Pandas/etc
- MS2 and MS3 quantitation results will be stored as a tab-separated file (`tmt.tsv`, `lfq.tsv`) if `quant.tmt` or `quant.lfq` options are used in the parameter file

If `--parquet` is passed as a command line argument, `results.sage.parquet` (and optionally, `lfq.parquet`) will be written. These have a similar set of columns, but TMT values are stored as a nested array alongside PSM features

## Configuration file schema

### Notes

- The majority of parameters are optional - only "database.fasta", "precursor_tol", and "fragment_tol" are required. Sage will try and use reasonable defaults for any parameters not supplied
- Tolerances are specified on the *experimental* m/z values. To perform a -100 to +500 Da open search (mass window applied to *theoretical*), you would use `"da": [-500, 100]`

### Decoys

Using decoy sequences is critical to controlling the false discovery rate in proteomics experiments. Sage can use decoy sequences in the supplied FASTA file, or it can generate internal sequences. Sage reverses tryptic peptides (not proteins), so that the [picked-peptide](https://pubmed.ncbi.nlm.nih.gov/36166314/) approach to FDR can be used.

If `database.generate_decoys` is set to true (or unspecified), then decoy sequences in the FASTA database matching `database.decoy_tag` will be *ignored*, and Sage will internally generate decoys. It is __critical__ that you ensure you use the proper `decoy_tag` if you are using a FASTA database containing decoys and have internal decoy generation turned on - otherwise Sage will treat the supplied decoys as hits!

Internally generated decoys will have protein accessions matching "{decoy_tag}{accession}", e.g. if `decoy_tag` is "rev_" then a protein accession like "rev_sp|P01234|HUMAN" will be listed in the output file.

### FASTA digestion

Sage will process a protein into peptides via several routes listed below. Currently, one and only one is supported.

- Enzymatic: `database.enzyme.cleave_at = "KR"` - configuration option set to a sequence of amino acids (e.g. "KR" for trypsin, "FWYL" for chymotrypsin)
- Non-enzymatic: `database.enzyme.cleave_at = ""` - All potential peptides between `min_len` and `max_len` will be generated from the sequence
- No digestion: `database.enzyme.cleave_at = "$"` - FASTA entries will be used as-is, subject to `min_len` and `max_len` options


### Example configuration file

For additional information about configuration options and output file formats, please see [the new documentation](https://sage-docs.vercel.app/docs)

```jsonc
// Note that json does not allow comments, they are here just as explanation
// but need to be removed in a real config.json file
{
  "database": {
    "bucket_size": 32768,           // How many fragments are in each internal mass bucket
    "enzyme": {               // Optional. Default is trypsin, using the parameters below
      "missed_cleavages": 2,  // Optional[int], Number of missed cleavages for tryptic digest
      "min_len": 5,           // Optional[int] {default=5}, Minimum AA length of peptides to search
      "max_len": 50,          // Optional[int] {default=50}, Maximum AA length of peptides to search
      "cleave_at": "KR",      // Optional[str] {default='KR'}. Amino acids to cleave at
      "restrict": "P",        // Optional[str] {default='P'}. Do not cleave if one of these AAs follows the cleavage site
      "c_terminal": false,      // Optional[bool] {default=true}. Cleave at c terminus of matching amino acid
      "semi_enzymatic": false      // Optional[bool] {default=false}. Generate semi-enzymatic peptides
    },
    "peptide_min_mass": 500.0,      // Optional[float] {default=500.0}, Minimum monoisotopic mass of peptides to fragment
    "peptide_max_mass": 5000.0,     // Optional[float] {default=5000.0}, Maximum monoisotopic mass of peptides to fragment
    "ion_kinds": ["b", "y"],        // Optional[List[str]] {default=["b","y"]} Which fragment ions to generate and search?
    "min_ion_index": 2,     // Optional[int] {default=2}, Do not generate b1/b2/y1/y2 ions for preliminary searching. Does not affect full scoring of PSMs
    "static_mods": {        // Optional[Dict[char, float]] {default={}}, static modifications
      "^": 304.207,         // Apply static modification to N-terminus of peptide
      "K": 304.207,         // Apply static modification to lysine
      "C": 57.0215          // Apply static modification to cysteine
    },
    "variable_mods": {    // Variable modification masses or objects with mass and optional max_count
      "M": [15.9949],     // Variable mods are applied *before* static mod
      "K": [{"mass": 42.0106, "max_count": 1}, 14.0157],
      "^Q": [-17.026549],
      "^E": [-18.010565], // Applied to N-terminal glutamic acid
      "$": [49.2, 22.9],  // Applied to peptide C-terminus
      "[": [42.0],          // Applied to protein N-terminus
      "]": [111.0]          // Applied to protein C-terminus
    },
    "max_variable_mods": 2, // Optional[int] {default=2} Limit modifications on each peptide
    "max_combinations": 8,  // Optional[int] {default=null} Limit total variants per peptide
    "decoy_tag": "rev_",    // Optional[str] {default="rev_"}: See notes above
    "generate_decoys": false, // Optional[bool] {default="true"}: Ignore decoys in FASTA database matching `decoy_tag`
    "fasta": "dual.fasta"   // str: mandatory path to FASTA file
  },
  "quant": {                // Optional - specify only if TMT or LFQ
    "tmt": "Tmt16",         // Optional[str] {default=null}, one of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18"
    "tmt_settings": {
      "level": 3,           // Optional[int] {default=3}, MS-level to perform TMT quantification on
      "sn": false           // Optional[bool] {default=false}, use Signal/Noise instead of intensity for TMT quant. Requires noise values in mzML
    },
    "lfq": true,            // Optional[bool] {default=null}, perform label-free quantification
    "lfq_settings": {
      "peak_scoring": "Hybrid", // See DOCS.md for details - recommend that you do not change this setting
      "integration": "Sum",   // Optional["Sum" | "Apex"], use sum of MS1 traces in peak, or MS1 intensity at peak apex
      "spectral_angle": 0.7,  // Optional[float] {default = 0.7}, normalized spectral angle cutoff for calling an MS1 peak
      "ppm_tolerance": 5.0,    // Optional[float] {default = 5.0}, tolerance (in p.p.m.) for DICE window around calculated precursor mass
      // Optional[bool] {default = true}. Combine all charge states for quantification. Setting this to false
      // quantifies each peptide-charge precursor in `precursor_charge` range (see below) separately
      "combine_charge_states": true
    }
  },
  "precursor_tol": {        // Tolerance can be either "ppm" or "da"
    "da": [
      -500,                 // This value is substracted from the experimental precursor to match theoretical peptides
      100                   // This value is added to the experimental precursor to match theoretical peptides
    ]
  },
  "fragment_tol": {         // Tolerance can be either "ppm" or "da"
    "ppm": [
     -10,                   // This value is subtracted from the experimental fragment to match theoretical fragments 
     10                     // This value is added to the experimental fragment to match theoretical fragments 
    ]
  },
  // Optional[Tuple[int, int]] {default=[2, 4]}
  // If charge states are not annotated in the mzML, or if `wide_window` mode is turned on, then consider
  // all precursors at z=2, z=3, z=4
  "precursor_charge": [2, 4]
  "isotope_errors": [       // Optional[Tuple[int, int]] {default=[0,0]}: C13 isotopic envelope to consider for precursor
    -1,                     // Consider -1 C13 isotope
    3                       // Consider up to +3 C13 isotope (-1/0/1/2/3) 
  ],
  "deisotope": false,       // Optional[bool] {default=true}: perform deisotoping and charge state deconvolution
  "chimera": false,         // Optional[bool] {default=false}: search for chimeric/co-fragmenting PSMS
  "wide_window": false,     // Optional[bool] {default=false}: _ignore_ `precursor_tol` and search in wide-window/DIA mode
  "predict_rt": false,    // Optional[bool] {default=true}: use retention time prediction model as an feature for LDA
  "min_peaks": 15,          // Optional[int] {default=15}: only process MS2 spectra with at least N peaks
  "max_peaks": 150,         // Optional[int | "auto"] {default=150}: take the top N most intense MS2 peaks to search;
                            // "auto" = 400 for ppm, 80 for a Da `fragment_tol` reaching 0.1 Da (see below)
  "min_matched_peaks": 6,   // Optional[int] {default=4}: minimum # of matched b+y ions to use for reporting PSMs
  "max_fragment_charge": 1, // Optional[int] {default=null}: maximum fragment ion charge states to consider,
  "report_psms": 1,         // Optional[int] {default=1}: number of PSMs to report for each spectra. Higher values might disrupt PSM rescoring.
  "ion_model": false,       // Optional[bool] {default=false}: experimental self-trained, cross-fitted fragment-ion model;
                            // adds the `ion_llr` and `ion_explained` columns to the results and the PIN file (for Percolator/mokapot)
  "output_directory": "s3://bucket/prefix" // Optional[str] {default=`.`}: Place output files in a given directory or S3 bucket/prefix
  "mzml_paths": [           // List[str]: representing paths to mzML (or gzipped-mzML) files for search
    "local/path.mzML",
    "s3://bucket/PXD0000001/foo.mzML.gz"
  ]       
}
```

## Using the docker image

Sage can be used from a docker image!

```shell
$ docker pull ghcr.io/lazear/sage:master
$ docker run -it --rm -v ${PWD}:/data ghcr.io/lazear/sage:master sage -o /data /data/config.json
# The sage executable is located in /app/sage in the image
```

> `-v ${PWD}:/data` means it will mount your current directory as `/data`
> in the docker image. Make sure all the paths in your command and configuration
> use the location in the image and not your local directory

# Further Details

This documentation covers the parameters in the JSON configuration file for the proteomics search engine. The configuration file contains information about the search engine's settings, including database, enzyme, modifications, and other settings. For a complete example of a configuration file, please see the [online docs](https://sage-docs.vercel.app/docs)

## Database

- **bucket_size**: Integer. The number of fragments in each internal mass bucket (default: 8192). Tweaking this parameter can increase search performance for wide precursor or fragment searches.
- **prefilter**: Boolean. Build the database in chunks of the FASTA file, search all spectra against each chunk, keep only the peptides that rank among the best candidates of some spectrum, and build the final fragment index from those peptides (default: false). Meant for databases whose full fragment index does not fit into memory; see [Memory, speed and the fragment index](#memory-speed-and-the-fragment-index).
- **prefilter_chunk_size**: Integer. Number of FASTA proteins per prefilter chunk (default: 0 = chosen so that a chunk holds about 8.4 million peptides, estimated from the unmodified digest and the variable modifications).
- **prefilter_low_memory**: Boolean. Keep only the best `report_psms` + 1 fully scored candidates of each spectrum (true), or every candidate of the preliminary search (false: a larger final index) (default: true).

### Enzyme

The enzyme section contains parameters related to the enzyme used for digestion. The default enzyme is trypsin, with the parameters specified below.

- **missed_cleavages**: Integer. The number of missed cleavages for tryptic digest (default: 1).
- **min_len**: Integer. The minimum amino acid (AA) length of peptides to search (default: 5).
- **max_len**: Integer. The maximum AA length of peptides to search (default: 50).
- **cleave_at**: String. Amino acids to cleave at (default: 'KR').
- **restrict**: String. Do not cleave if one of these amino acids follows the cleavage site. Default: 'P' when `cleave_at` is omitted (trypsin); no restriction when a custom `cleave_at` is given without `restrict`. Use `""` or `null` for no restriction.
- **c_terminal**: Boolean. Cleave at the C-terminus of matching amino acids (default:true).

Example: 
```json
"database": {
  "enzyme": {
    "missed_cleavages": 1,
    "min_len": 5,
    "max_len": 50,
    "cleave_at": "KR",
    "restrict": "P",
    "c_terminal": true
  }
}
```

### Fragment Settings

- **peptide_min_mass**: Float. The minimum monoisotopic mass of peptides to fragment *in silico* (default: 500.0).
- **peptide_max_mass**: Float. The maximum monoisotopic mass of peptides to fragment *in silico* (default: 5000.0).
- **ion_kinds**: List of strings. Which fragment ions to produce? Allowed values: "a", "b", "c", "x", "y", "z". (default: ["b", "y"])
- **min_ion_index**: Integer. Do not generate b1/bN/y1/yN ions for preliminary searching if `min_ion_index = N`. Does not affect full scoring of PSMs (default: 2).

Example:
```json
"database": {
  "peptide_min_mass": 500.0,
  "peptide_max_mass": 5000.0,
  "ion_kinds": ["b", "y"],
  "min_ion_index": 2
}
```

### Modifications

#### Static Modifications

- **static_mods**: Dictionary with characters as keys and floats as values. Represents static modifications applied to amino acids or termini (default: {}). Static modifications are applied after variable modifications
  - Example: Apply a static modification of 304.207 to the N-terminus of the peptide and lysine, and 57.0215 to cysteine.
    ```json
    "database": {
      "static_mods": {
        "^": 304.207,
        "K": 304.207,
        "C": 57.0215
      }
    }
    ```

#### Variable Modifications

- **max_variable_mods**: Integer. Limit the total variable modifications on each peptide (default: 2).
- **max_combinations**: Integer. Optional hard cap on the total variants generated per input peptide, including the unmodified form. Variants with fewer modifications are retained first. Values below 1 are treated as 1 (default: unlimited).
- **variable_mods**: Dictionary with characters as keys and lists containing bare masses or objects with a `mass` field and optional `max_count`. A bare mass remains unrestricted except by `max_variable_mods`; `max_count` limits occurrences of that specific modification on one peptide.
  - Example: Apply a variable modification of 15.9949 to methionine, 49.2022 to the C-terminus of the peptide, 42.0 to the N-terminus of the protein, and 111.0 to the C-terminus of the protein.
    ```jsonc
    "database": {
      "variable_mods": {
        "M": [15.9949],
        "K": [{"mass": 42.0106, "max_count": 1}, 14.0157],
        "^Q": [-17.026549],
        "^E": [-18.010565], // Applied to N-terminal glutamic acid
        "$": [49.2022],     // Applied to peptide C-terminus
        "[": 42.0,          // Applied to protein N-terminus
        "]": 111.0          // Applied to protein C-terminus
      }
    }
    ```
  - Syntax:
    "^X": Modification to be applied to amino acid X if it appears at the N-terminus of a peptide
    "$X": Modification to be applied to amino acid X if it appears at the C-terminus of a peptide
    "[X": Modification to be applied to amino acid X if it appears at the N-terminus of a protein
    "]X": Modification to be applied to amino acid X if it appears at the C-terminus of a protein

### Decoys

- **decoy_tag**: String. The tag used to identify decoy entries in the FASTA database (default: "rev_").
- **generate_decoys**: Boolean. If true, ignore decoys in the FASTA database matching `decoy_tag`, and generate internally reversed peptides (default: true).

### FASTA

- **fasta**: String. The path to the FASTA file, either a local path or s3 object URI.

## Quantification

The quant section is optional and should be specified only if TMT or LFQ is used.


- **tmt**: String. One of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18" (default: null).
- **tmt_settings**: Object containing TMT-specific settings.
  - **level**: Integer. The MS-level to perform TMT quantification on (default: 3).
  - **sn**: Boolean. Use Signal/Noise instead of intensity for TMT quantification. Requires noise values in mzML (default: false).
- **lfq**: Boolean. Perform label-free quantification (default: null).
- **lfq_settings**: Object containing LFQ-specific settings.
  - **peak_scoring**: String. The method used for scoring peaks in LFQ, one of: "Hybrid", "RetentionTime", "SpectralAngle" (default: "Hybrid").
  - **integration**: String. The method used for integrating peak intensities, either "Sum" or "Max" (default: "Sum").
  - **spectral_angle**: Float. Threshold for the spectral angle similarity measure, ranging from 0 to 1 (default: 0.7).
  - **ppm_tolerance**: Float. Tolerance for matching MS1 ions in parts per million (default: 5.0).

Example: 
```json
 "quant": {
    "tmt": "Tmt16",
    "tmt_settings": {
      "level": 3,
      "sn": false
    },
    "lfq": true,
    "lfq_settings": {
      "peak_scoring": "Hybrid",
      "integration": "Sum",
      "spectral_angle": 0.7,
      "ppm_tolerance": 5.0
    }
  }
```


## Precursor Tolerance

- **precursor_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-500, 100] in daltons.
    ```json
    "precursor_tol": {
      "da": [-500, 100]
    }
    ```

## Fragment Tolerance

- **fragment_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-10, 10] in parts per million.
    ```json
    "fragment_tol": {
      "ppm": [-10, 10]
    }
    ```

## Isotope Errors

- **isotope_errors**: List of two integers. The C13 isotopic envelope to consider for precursor (default: [0, 0]).
  - Example: Consider -1 and up to +3 C13 isotopes (-1/0/1/2/3).
    ```json
    "isotope_errors": [-1, 3]
    ```

**NOTE**: Searching with isotope errors is slower than searching with a wider precursor tolerance that encompasses the isotope errors, e.g. `"da": [-3.5, 1.25]`. Using the wider precursor tolerance will generally increase the number of confidently identified PSMs as well.

## Other Settings

Note on the settings below:

`predict_rt` is incompatible with `quant.lfq = true`. Setting `quant.lfq = true` will automatically turn on global retention time alignment and prediction, which are crucial for accurate direct ion current extraction.

- **deisotope**: Boolean. Perform deisotoping and charge state deconvolution on MS2 spectra (default: true). Recommended for high-resolution MS2 scans. This setting may interfere with TMT-MS2 quantification, use at your own risk.
- **chimera**: Boolean. Search for chimeric/co-fragmenting PSMs (default: false).
- **wide_window**: Boolean. Ignore `precursor_tol` and search spectra in wide-window/dynamic precursor tolerance mode (default: false).
- **predict_rt**: Boolean. Use retention time prediction model as a feature for LDA (default: true).
- **min_peaks**: Integer. Only process MS2 spectra with at least N peaks (default: 15).
- **max_peaks**: Integer or `"auto"`. Take the top N most intense MS2 peaks (after deisotoping) to search (default: 150). With a number below `min_peaks`, no MS2 spectrum can be searched; Sage logs a warning.

  `"auto"` (opt-in) chooses the number from `fragment_tol`:
  - 400 peaks for a ppm or pct tolerance, or a Da tolerance narrower than 0.1 Da (e.g. 0.02 Da). High-resolution spectra carry real fragment ions well below the 150th most intense peak.
  - 80 peaks for a Da tolerance that reaches at least 0.1 Da on one side (low-resolution fragment spectra, e.g. ion trap CID searched at 0.5 Da). At a wide tolerance, low-intensity peaks mostly add chance matches. With TMT quantification at MS2 (`quant.tmt_settings.level` 2) it keeps 150 instead, because the reporter ions are read from the capped spectrum and 80 peaks drop some of them.
  - 150 with `wide_window: true`, whatever the tolerance. In wide-window and DIA-style searches, 400 peaks let wrong candidates from the wide isolation window win more often. On HF-X, timsTOF and Astral files of the benchmark below, Percolator found 3-46% fewer PSMs with 400 peaks than with 150 (Astral -46%). On a public diaPASEF run (PXD017703), it found 12-16% fewer, and the run took 1.9x as long.
  - Never fewer than `min_peaks`: a smaller number is raised to `min_peaks`.

  The resolved number is logged (e.g. `max_peaks: auto -> 400 peaks per MS2 spectrum (...)`) and written to `results.json`. `"auto"` goes by the unit of `fragment_tol` only. It assumes that a Da tolerance of 0.1 Da or more means low-resolution fragment spectra and that a ppm or pct tolerance means high-resolution ones. For other cases (e.g. low-resolution data searched with a ppm tolerance), set a number.

  Measured on the public benchmark of OpenMS issue #10364 (8,000 MS2 spectra from each of 20 public runs) and on a held-out set (another 8,000 spectra from each of the same runs), with Percolator 3.09, 20 seeds and target PSMs at 1% FDR, against `max_peaks: 150`:

  | | 20 benchmark files | held-out spectra |
  |---|---|---|
  | PSMs per file | +2.26% (Astral +10.8%, Lumos LFQ -0.1%) | +2.38% (Astral +7.7%, Lumos LFQ +0.4%) |
  | doubled-database entrapment FDP (three shuffles) | 1.08% -> 1.03% | 1.00% -> 1.07% |
  | Velos UPS1 entrapment FDP | 1.20% -> 1.25% | 1.06% -> 0.95% |

  On the held-out set, the doubled-database FDP rose by 0.070 points, above the pre-registered limit of +0.05 for a default change, so `"auto"` stayed opt-in. The excess is about 1.5 standard errors and is not significant, and the two sets disagree in sign. Sage's own q-values found 2.2% more PSMs on the benchmark files.

  Cost, measured against `max_peaks: 150` (AMD EPYC 7763, paired runs):
  - Time: the search phase scores more peaks. On Astral spectra it takes 70-80% longer (0.46-0.48 -> 0.82-0.84 s per 8,000 spectra at 4 threads); on 0.5 Da data at 80 peaks it is shorter.
    - The 20 benchmark files at 4 threads took about the same total time (+0.2%; high-resolution files -1.4% to +7.7%, 0.5 Da files -3.4% to -7.0%), because building the fragment index takes most of these runs.
    - Larger inputs show more of it: 24,000 Astral spectra took 11.7%, 10.2% and 7.4% longer at 4, 16 and 64 threads, and 24,000 Velos spectra (80 peaks) took 11.0%, 8.6% and 10.5% less.
  - Memory: the extra peaks (9 bytes each) are kept for every MS2 spectrum of a batch until the batch has been searched.
    - That is about +0.30 GiB per full-size Astral run (about 120,000 MS2 spectra) in the batch, and +0.50 GiB with `report_psms` 10. This was measured with three and six such runs in one batch; it grows linearly with the number of files.
    - By default a batch holds CPUs/2 files (`--batch-size`), so on a large machine it holds all files of a run. `--batch-size` or `max_peaks: 150` bounds the extra memory.
    - A single gzipped file shows little of it (+0.12 GiB), because there the read and the index build set the peak.
- **min_matched_peaks**: Integer. The minimum number of matched b+y ions to use for reporting PSMs (default: 4).
- **max_fragment_charge**: Integer. The maximum fragment ion charge states to consider (default: null - use precursor z-1).
- **report_psms**: Integer. The number of PSMs to report for each spectrum. Higher values might disrupt LDA (default: 1).
- **ion_model**: Boolean. Learn a fragment-ion likelihood model from each file's own confident PSMs and report two extra features per PSM, `ion_llr` and `ion_explained`, in results.sage.tsv (or results.sage.parquet) and in the PIN file (default: false). They are meant for a rescorer such as Percolator or mokapot; Sage's own LDA does not use them. Experimental: its error control was not confirmed on held-out data, so it is not recommended for routine use. See [Fragment-ion model](#fragment-ion-model-ion_model).
- **parallel**: Boolean. Parse and search files in parallel. For large numbers of files or low RAM, setting this to false can reduce memory usage at the cost of running slower (default: true).

## mzML Paths

- **mzml_paths**: List of strings. The paths to mzML (or gzipped-mzML) files for search. Paths are either local, or point to an S3 object. Files ended in ".gz" or ".gzip" are inferred to be compressed.
  - Example:
    ```json
    "mzml_paths": [
      "local/path.mzML",
      "s3://my-mass-spec-data/PXD0000001/foo.mzML.gz"
    ]
    ```
  - Local, uncompressed mzML files are parsed in parallel, in chunks of whole spectra (the spectra are the same as with the serial parse).
    - Gzipped and remote files, and files with comments, CDATA sections or processing instructions between the spectra or another unusual layout, are parsed serially. They are read completely (all raw spectra of a file) before their spectra are processed.
    - Setting the environment variable `SAGE_MZML_SERIAL` to any value except `0` (e.g. `SAGE_MZML_SERIAL=1`) forces the serial parse. Why a file is parsed serially is logged at debug level (`SAGE_LOG=sage_cloudpath=debug`).
    - Speed: on 72,000 Astral spectra (AMD EPYC 7763), the parallel parse alone made a run 1.42x faster at 16 threads and 1.6x faster at 64-128. On 8,000-spectrum files it is about neutral.
    - The parallel parse holds only the chunks being parsed in memory. It reads every file twice: first a scan for the spectrum boundaries, then the chunks. The second read normally comes from the page cache.
    - If the files of a batch (all files read together, by default CPUs/2 of them) do not fit into free memory, the second read goes to the storage again. On storage slower than about 300 MB/s (a hard disk, a 1 GbE network share), the parallel parse can then be slower than the serial one; set `SAGE_MZML_SERIAL=1` or a smaller `--batch-size` there. On local NVMe and on Ceph storage, it stayed 3.3-4.5x faster than the serial parse even when the file was read twice.
  
## Output directory:

- **output_directory**: Local directory, or S3 location where output files will be written. If the local directory does not already exist, it will be created. Write permissions are required for the directory or S3 path.
  - Possible output files are: "results.json", "results.sage.tsv", "lfq.tsv", and "tmt.tsv"
  - Example:
  ```json
  "output_directory": "s3://my-mass-spec-results/PXD003881/"
  ```

## Memory, speed and the fragment index

- **Batches**: `--batch-size` files (default: CPUs/2) are read and searched together, and the processed spectra of all files of a batch stay in memory until the batch has been searched. On a large machine that is usually every file of a run. A smaller `--batch-size` lowers the peak memory of runs with many large files. That matters most with `max_peaks: "auto"` and `ion_model`, which keep more data per MS2 spectrum.
- **Pruned fragment index**: when every input file is read in the first batch and `database.prefilter` is off, the fragment index holds only the fragments of peptides that some precursor window of those spectra can reach. On the benchmark of OpenMS issue #10364, that drops 35-64% of the fragments.
  - The peptide list stays complete, and the results are the same as with the full index.
  - The build does not wait for the spectra. It prunes at the first of three points at which all spectra have been read: before counting the fragments, after counting them, or before bucketing them. Each step logs what it did.
  - The full index is built if the read is still going after the last point. That happens often with gzipped, remote or otherwise serially parsed files at high thread counts. The full index is also built when the files need more than one batch, and in prefilter mode.
  - The telemetry record's `fragments` field reports the number of fragments in the index, i.e. the pruned count.
- **Prefilter mode** (`database.prefilter: true`) is meant for databases whose full fragment index does not fit into memory. It searches the spectra against the database chunk by chunk (see [Database](#database)) and builds the final index from the peptides it keeps.
  - The fragment index of each chunk is returned to the operating system before the next chunk is digested.
  - On two of the benchmark files (Velos and Astral) at 64 threads, prefilter mode was 1.8x faster than v0.15.0-fork.2 and peaked 0.24-0.42 GiB lower.
- **Allocator**: Sage uses mimalloc, which keeps freed memory for a while before returning it to the operating system. Sage asks it to return freed memory before the peptide sort and before the fragment index is allocated.
  - With several gzipped (or otherwise serially parsed) files in one batch at 64-128 threads, where the index is often not pruned before it is allocated, peak memory occasionally rose 1.0-1.2 GiB above v0.15.0-fork.2. Since the second release point was added, this did not happen in 100 benchmark runs, but it did in one of 22 profiled runs (about 1 GiB above v0.15.0-fork.2).
  - `MIMALLOC_ARENA_EAGER_COMMIT=0` makes mimalloc commit memory as it is used, which lowers peak memory at high thread counts at some cost in speed.

# Interpreting Sage Output

The "results.sage.tsv" file contains the following columns (headers):

- `peptide`: Peptide sequence, including modifications (e.g., NC\[+57.021\]HKGSFK).
- `proteins`: Proteins containing the peptide sequence.
- `num_proteins`: Number of proteins assigned to the peptide sequence.
- `filename`: File containing this PSM
- `scannr`: Spectrum identifier from mzML file.
- `rank`: Rank of the PSM. If `report_psms > 1`, then the best match will have rank = 1, the second best match will have rank = 2, etc. 
- `label`: Target/Decoy label (-1: decoy, 1: target).
- `expmass`: Experimental mass of the peptide.
- `calcmass`: Calculated mass of the peptide.
- `charge`: Reported precursor charge.
- `pepide_len`: Length of the peptide sequence.
- `missed_cleavages`: Number of missed cleavages.
- `isotope_error`: C13 isotope error.
- `precursor_ppm`: Difference between experimental mass and calculated mass, reported in parts-per-million.
- `fragment_ppm`: Average parts-per-million (delta mass) for matched fragment ions compared to theoretical ions.
- `hyperscore`: X!Tandem hyperscore for the PSM.
- `delta_next`: Difference between the hyperscore of this candidate and the next best candidate.
- `delta_bext`: Difference between the hyperscore of the best candidate (rank=1) and this candidate.
- `rt`: Retention time.
- `aligned_rt`: Globally aligned retention time.
- `predicted_rt`: Predicted retention time, if enabled.
- `delta_rt_model`: Difference between predicted and observed retention time.
- `matched_peaks`: Number of matched theoretical fragment ions.
- `longest_b`: Longest b-ion series.
- `longest_y`: Longest y-ion series.
- `longest_y_pct`: Longest y-ion series, divided by peptide length (as a percentage).
- `matched_intensity_pct`: Fraction of MS2 intensity explained by matched b- and y-ions (as a percentage of total MS2 intensity for this spectrum).
- `scored_candidates`: Number of scored candidates for this spectrum.
- `poisson`: Probability of matching exactly N peaks across all candidates (Pr(x=k)).
- `sage_discriminant_score`: Combined score from linear discriminant analysis, used for FDR (False Discovery Rate) calculation.
- `posterior_error`: log10 of the posterior error probability (local FDR) of this PSM, i.e. 0 means PEP = 1 and -2 means PEP = 0.01. Estimated as min(1, p / (1 - p)) from the KDE's probability p that a PSM with this discriminant score is a decoy. 0 when the LDA model could not be fitted.
- `spectrum_q`: Assigned spectrum-level q-value.
- `peptide_q`: Assigned peptide-level q-value.
- `protein_q`: Assigned protein-level q-value.
- `ms1_intensity`: Intensity of the selected MS1 precursor ion (not label-free quant)
- `ms2_intensity`: Total intensity of MS2 spectrum
- `ion_llr`, `ion_explained`: Fragment-ion model features, only with `ion_model: true` (see below). They are also written to the PIN file, before `Peptide`.

These columns provide comprehensive information about each candidate peptide spectrum match (PSM) identified by the Sage search engine.

## Fragment-ion model (`ion_model`)

With `"ion_model": true`, Sage learns, for every input file, how the fragment ions of confidently identified peptides show up in that file's spectra, and scores every PSM against it. In short: a PSM whose b/y ions are found (or missing) the way they are for correct peptides gets a high `ion_llr`, one whose ions look like those of reversed sequences a low one, and each PSM is scored by a model trained on the other half of the file's spectra. The option is experimental and off by default. Its two features are only useful with a rescorer (Percolator, mokapot), because Sage's LDA does not use them. Details:

- **Evidence**: all deisotoped peaks of the MS2 spectrum (not only the `max_peaks` most intense ones), with their intensity rank.
- **Contexts and outcomes**: each theoretical b/y ion (fragment charge 1, and also 2 for precursor charge >= 3) has a context (ion series, precursor charge 2 / 3 / >= 4, fragment charge, position along the peptide in tenths, cleavage N-terminal to proline or C-terminal to D/E, complementary ion found or not) and an outcome (not found, or found at one of 7 intensity rank bins (1-2, 3-5, 6-10, 11-20, 21-40, 41-80, > 80) x 3 mass error bins (< 1/4, < 1/2, <= 1 of the window's half-width); the nearest peak counts).
  - The window is centred on the theoretical m/z, with half the total width of `fragment_tol` on each side: ±(hi - lo)/2. For a symmetric tolerance this is Sage's own matching window.
  - For an asymmetric `fragment_tol` it differs from Sage's windows, which themselves differ between the preliminary search and the rescoring. That is deliberate. On the benchmark, following the rescoring window instead lost 1.7% of the PSMs with ppm(-10, 30) and 28% with Da(-0.5, 0) when the peaks were centred near 0. It gained nothing measurable where the asymmetry matched a real mass offset.
- **Training**: rank-1 target PSMs at 1% FDR of a label-free target-decoy competition on the HyperScore computed with intensities normalised to the spectrum's most intense peak, separately for precursor charge <= 2 and >= 3. A *signal* table counts the outcomes of the identified peptides, a *noise* table those of the same peptides reversed except for the C-terminal residue, matched against the same spectra. Both are smoothed by back-off over coarser contexts (pseudo-count 20).
- **Features**: `ion_llr` = sum over the ions of log P(outcome | signal) - log P(outcome | noise); `ion_explained` = share of the ions the signal model expects to be found that were found (weighted by that probability).
- **Cross-fitting**: the spectra of a file are split into two folds by the parity of their index (among the spectra with PSMs, in file order). Each fold trains a model, and every PSM is scored by the model of the other fold, so no PSM is scored by a model trained on it.
  - The folds are assigned per spectrum, not per peptide. A peptide identified in spectra of both folds is therefore scored partly by a model trained on its other observations. On the held-out data, 7% of the rank-1 target PSMs are in that situation. Leaving such peptides out of the other fold's training changed the yield by +0.06% and the doubled-database FDP by +0.008 points, so no measurable effect.
  - A file gets features only if both folds have at least 100 training PSMs and the file has a rank-1 decoy PSM. Otherwise every PSM of the file gets 0, and a warning is logged. The gate is per file, as in ProSE: one bit per file, the same for targets and decoys. All 40 benchmark and held-out files were trained, with at least a 4x margin.
- **Cost**: time and memory grow with the number of MS2 spectra and peaks.
  - Time: building the peak lists while the files are read, then training and applying the model, took about 15-25 µs per MS2 spectrum at 4 threads (`report_psms` 10). On the release build (with `max_peaks: "auto"` in both arms), that was +2.4% wall time on the 20 benchmark files at 4 threads (per file +0.1% to +4.4%), and +2.5% / +4.5% on 24,000 Astral spectra at 16 / 64 threads.
  - Memory: the peak lists take about 5 bytes per deisotoped peak (about 4.5 KB per Astral MS2 spectrum). They are kept for every MS2 spectrum of a batch until the batch has been searched, about 40 MB per 8,000 Astral spectra and 0.35 GiB per 72,000.
    - Peak memory rose by 0.13 GiB on 24,000 Astral spectra.
    - A batch of twelve 72,000-spectrum Astral files peaked at 28.5 instead of 24.4 GiB.
    - By default a batch holds CPUs/2 files; a smaller `--batch-size` bounds the extra memory.

Measured on the public benchmark of OpenMS issue #10364 (8,000 MS2 spectra from each of 20 public runs) and on a held-out set (another 8,000 spectra from each of the same runs), with Percolator 3.09, 20 seeds, target PSMs at 1% FDR and `max_peaks` 150:

| | 20 benchmark files | held-out spectra |
|---|---|---|
| PSMs per file | +5.2% (Astral +24%, Velos +5.1%, HF-X +4.6%, TMT +3.4-3.6%, Lumos LFQ and timsTOF +0.2-0.3%) | +5.5% (Astral +25%) |
| doubled-database entrapment FDP (three shuffles) | 1.08% -> 1.09% | 1.00% -> 1.03% |
| Velos UPS1 entrapment FDP | 1.20% -> 1.27% | 1.06% -> 1.30% |

- On the held-out set, the Velos UPS1 FDP rose in all three Velos files and in all 20 seeds, above the pre-registered limit of 1.16%.
- With `max_peaks: "auto"`, the gain was +6.4% (benchmark) and +6.0% (held-out). But the held-out doubled-database FDP rose by 0.077 points, above the limit of +0.05.
- The error control of the two features is therefore not confirmed. The option stays experimental and is not recommended for routine use.

All other columns are the same with and without the option, and the output does not depend on the number of threads.
