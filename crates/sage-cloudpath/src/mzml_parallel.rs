//! Parallel parsing of local, uncompressed mzML files.
//!
//! The serial reader ([`MzMLReader::parse`]) scans and decodes a whole file on one thread
//! (~120-150 MB/s). On large inputs that is the critical path from 16 threads on: 8 s of an
//! 11 s run for 72,000 Astral spectra at 64 threads. Here the file is cut into chunks of whole
//! `<spectrum>` elements and every chunk is parsed by the unchanged [`MzMLReader`] as
//! `header + chunk + footer`, in parallel, and processed right away:
//!
//! 1. scan: pieces of the file are read with positioned reads in parallel and searched for
//!    `<spectrum` start tags;
//! 2. header: the bytes before the first spectrum (with the referenceableParamGroups) and a
//!    footer that closes the elements the header opens;
//! 3. parse: chunks of whole spectra, each read by the task that parses it and passed through
//!    `f`; the results are concatenated in file order. The last chunk runs to the real end of
//!    the file, so a truncated file is detected as before.
//!
//! Memory is bounded by the chunks in flight (one per worker, at most [`MAX_CHUNK`] bytes and
//! their raw spectra) instead of the whole file's raw spectra.
//!
//! The spectra are those of the serial parse: the parser keeps no state from one spectrum to
//! the next except the header's param groups (every per-spectrum field is reset at
//! `<spectrum>`), so a chunk parsed after the header starts in the state the serial parse is in
//! between two spectra. Anything this module cannot cut safely goes to the serial reader
//! instead: no `<spectrum` start tag, a header that does not end with the `<spectrumList>` start
//! tag, a start tag that does not directly follow `</spectrum>`, comments, CDATA sections or
//! processing instructions after the header, read errors, and every parse error, so that
//! errors (a truncated file, a corrupt array, MS-Numpress) are reported exactly as before.
use crate::mzml::MzMLReader;
use rayon::prelude::*;
use sage_core::spectrum::RawSpectrum;
use std::fs::File;
use std::path::Path;

/// Bytes per piece of the boundary scan
const PIECE: usize = 2 << 20;
/// Bytes read before each piece: enough for `</spectrum>` and the whitespace before the next
/// start tag (a start tag with more whitespace in front goes to the serial reader)
const LEAD: usize = 64;
/// Bytes read after each piece, so that a start tag beginning in the piece is seen whole
const TRAIL: usize = 16;
/// Chunk sizes: about 8 chunks per worker, within these bounds
const MIN_CHUNK: usize = 256 << 10;
const MAX_CHUNK: usize = 8 << 20;
/// Every chunk re-parses the header; larger headers go to the serial reader
const MAX_HEADER: usize = 256 << 10;

const START: &[u8] = b"<spectrum";
const END: &[u8] = b"</spectrum>";

/// Sizes used by [`read_processed`]; tests use small ones to cut files into many pieces/chunks
#[derive(Copy, Clone)]
pub(crate) struct Sizes {
    pub piece: usize,
    /// bytes per chunk; `None`: from the file size and the thread count
    pub chunk: Option<usize>,
}

impl Default for Sizes {
    fn default() -> Self {
        Self {
            piece: PIECE,
            chunk: None,
        }
    }
}

fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r')
}

fn trim_end(mut b: &[u8]) -> &[u8] {
    while let [rest @ .., last] = b {
        if !is_space(*last) {
            break;
        }
        b = rest;
    }
    b
}

#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    std::os::unix::fs::FileExt::read_exact_at(file, buf, offset)
}

#[cfg(windows)]
fn read_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> std::io::Result<()> {
    use std::os::windows::fs::FileExt;
    while !buf.is_empty() {
        match file.seek_read(buf, offset) {
            Ok(0) => return Err(std::io::ErrorKind::UnexpectedEof.into()),
            Ok(n) => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

#[cfg(not(any(unix, windows)))]
fn read_at(_: &File, _: &mut [u8], _: u64) -> std::io::Result<()> {
    Err(std::io::ErrorKind::Unsupported.into())
}

/// What the scan found in one piece
#[derive(Default)]
struct Piece {
    /// file offsets of `<spectrum` start tags, each with: does it directly follow `</spectrum>`
    /// (whitespace aside)?
    starts: Vec<(u64, bool)>,
    /// file offset of the last `<!` or `<?` (comment, CDATA, DOCTYPE, processing instruction)
    unusual: Option<u64>,
}

fn scan_piece(file: &File, len: u64, piece: usize, k: usize, buf: &mut Vec<u8>) -> Option<Piece> {
    let lo = (k * piece) as u64;
    let hi = (lo + piece as u64).min(len);
    let from = lo.saturating_sub(LEAD as u64);
    let to = (hi + TRAIL as u64).min(len);
    buf.resize((to - from) as usize, 0);
    read_at(file, buf, from).ok()?;
    let (begin, end) = ((lo - from) as usize, (hi - from) as usize);

    let mut out = Piece::default();
    // every `<` is markup: base64 has none, and text and attribute values must escape it
    for i in memchr::memchr_iter(b'<', &buf[begin..end]).map(|i| i + begin) {
        match buf.get(i + 1) {
            Some(b'!' | b'?') => out.unusual = Some(from + i as u64),
            _ if buf[i..].starts_with(START)
                && buf.get(i + START.len()).is_some_and(|&c| is_space(c)) =>
            {
                out.starts
                    .push((from + i as u64, trim_end(&buf[..i]).ends_with(END)))
            }
            _ => {}
        }
    }
    Some(out)
}

/// The closing tags of the elements the header leaves open, innermost first, or `None` if the
/// header does not end with the `<spectrumList>` start tag
fn footer(header: &[u8]) -> Option<Vec<u8>> {
    use quick_xml::events::Event;
    let h = trim_end(header);
    let tag = &h[memchr::memrchr(b'<', h)?..];
    let name = b"<spectrumList";
    let named = tag.starts_with(name)
        && tag
            .get(name.len())
            .is_some_and(|&c| is_space(c) || c == b'>');
    if !named || !tag.ends_with(b">") || tag.ends_with(b"/>") {
        return None;
    }
    let mut reader = quick_xml::Reader::from_reader(header);
    let (mut open, mut buf) = (Vec::new(), Vec::new());
    loop {
        match reader.read_event_into(&mut buf).ok()? {
            Event::Start(e) => open.push(e.name().as_ref().to_vec()),
            Event::End(_) => {
                open.pop()?;
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    let mut out = Vec::new();
    for name in open.iter().rev() {
        out.extend_from_slice(b"</");
        out.extend_from_slice(name);
        out.push(b'>');
    }
    Some(out)
}

/// Spectra of a local, uncompressed mzML file, parsed in parallel and passed through `f`, in
/// file order: the same as `MzMLReader::parse` followed by `f`. `Err(reason)`: the caller reads
/// the file with the serial reader (see the module documentation).
pub(crate) fn read_processed<T: Send>(
    path: &Path,
    file_id: usize,
    sn: Option<u8>,
    skip_ms1: bool,
    sizes: Sizes,
    f: &(impl Fn(RawSpectrum) -> T + Sync),
) -> Result<Vec<T>, &'static str> {
    let file = File::open(path).map_err(|_| "open failed")?;
    let len = file.metadata().map_err(|_| "no metadata")?.len();

    // 1. every `<spectrum` start tag
    let pieces = (len as usize).div_ceil(sizes.piece);
    let found = (0..pieces)
        .into_par_iter()
        .map_init(Vec::new, |buf, k| {
            scan_piece(&file, len, sizes.piece, k, buf)
        })
        .collect::<Option<Vec<Piece>>>()
        .ok_or("read error")?;
    let first = found
        .iter()
        .find_map(|p| p.starts.first())
        .ok_or("no <spectrum> start tag")?
        .0;
    if found.iter().filter_map(|p| p.unusual).any(|u| u > first) {
        return Err("comment, CDATA or processing instruction after the header");
    }
    let starts = found
        .iter()
        .flat_map(|p| p.starts.iter())
        .enumerate()
        .map(|(i, &(s, follows_end))| (i == 0 || follows_end).then_some(s))
        .collect::<Option<Vec<u64>>>()
        .ok_or("a <spectrum> start tag does not follow </spectrum>")?;
    drop(found);
    if starts.len() < 2 {
        return Err("fewer than two spectra");
    }

    // 2. header and footer
    if first as usize > MAX_HEADER {
        return Err("header too large");
    }
    let mut header = vec![0; first as usize];
    read_at(&file, &mut header, 0).map_err(|_| "read error")?;
    let footer = footer(&header).ok_or("header does not end with <spectrumList>")?;

    // 3. chunks of whole spectra, parsed in parallel
    let threads = rayon::current_num_threads().max(1);
    let floor = MIN_CHUNK.max(32 * header.len());
    let target = sizes
        .chunk
        .unwrap_or_else(|| (len as usize / (threads * 8)).clamp(floor, MAX_CHUNK.max(floor)))
        as u64;
    let mut bounds = vec![first];
    for &s in &starts[1..] {
        if s - bounds[bounds.len() - 1] >= target {
            bounds.push(s);
        }
    }
    bounds.push(len);
    let chunks = bounds.len() - 1;
    let parsed = (0..chunks)
        .into_par_iter()
        .map_init(Vec::new, |input, k| {
            // header + chunk + footer, in a buffer the worker reuses
            let (a, b) = (bounds[k], bounds[k + 1]);
            input.clear();
            input.extend_from_slice(&header);
            input.resize(header.len() + (b - a) as usize, 0);
            read_at(&file, &mut input[header.len()..], a).ok()?;
            if k + 1 < chunks {
                input.extend_from_slice(&footer);
            }
            let mut reader = MzMLReader::with_file_id(file_id);
            reader.set_signal_to_noise(sn).set_skip_ms1(skip_ms1);
            let raw = futures::executor::block_on(reader.parse(&input[..])).ok()?;
            Some(raw.into_iter().map(f).collect::<Vec<T>>())
        })
        .collect::<Option<Vec<Vec<T>>>>()
        .ok_or("read or parse error")?;

    let mut out = Vec::with_capacity(parsed.iter().map(Vec::len).sum());
    for part in parsed {
        out.extend(part);
    }
    log::debug!(
        "{}: parsed {} spectra in {} parallel chunks",
        path.display(),
        out.len(),
        chunks
    );
    Ok(out)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::mzml::MzMLError;

    const FIXTURE: &str = include_str!("../../../tests/LQSRPAAPPAPGPGQLTLR.mzML");

    fn spectrum_element() -> &'static str {
        let start = FIXTURE.find("<spectrum ").unwrap();
        let end = FIXTURE.find("</spectrum>").unwrap() + "</spectrum>".len();
        &FIXTURE[start..end]
    }

    /// The fixture with `n` spectra: every 5th is MS1, the others alternate a charge state; `sep`
    /// goes between two spectra
    fn synthetic(n: usize, sep: &str) -> String {
        let one = spectrum_element();
        let head = &FIXTURE[..FIXTURE.find("<spectrum ").unwrap()];
        let tail = &FIXTURE[FIXTURE.find("</spectrum>").unwrap() + "</spectrum>".len()..];
        let mut out = head.replacen(
            r#"<spectrumList count="1""#,
            &format!(r#"<spectrumList count="{n}""#),
            1,
        );
        for i in 0..n {
            let mut s = one.replace("scan=30069", &format!("scan={}", 1000 + i));
            if i % 5 == 0 {
                s = s.replace(
                    r#"name="ms level" value="2""#,
                    r#"name="ms level" value="1""#,
                );
            }
            if i % 2 == 1 {
                s = s.replacen(
                    r#"name="charge state" value="3""#,
                    r#"name="charge state" value="2""#,
                    1,
                );
            }
            if i > 0 {
                out.push_str(sep);
            }
            out.push_str(&s);
        }
        out.push_str(tail);
        out
    }

    fn temp(name: &str, content: &str) -> std::path::PathBuf {
        let path =
            std::env::temp_dir().join(format!("sage-par-{}-{}.mzML", std::process::id(), name));
        std::fs::write(&path, content).unwrap();
        path
    }

    fn serial(
        content: &str,
        skip_ms1: bool,
        sn: Option<u8>,
    ) -> Result<Vec<RawSpectrum>, MzMLError> {
        futures::executor::block_on(
            MzMLReader::with_file_id(3)
                .set_skip_ms1(skip_ms1)
                .set_signal_to_noise(sn)
                .parse(content.as_bytes()),
        )
    }

    fn parallel(
        name: &str,
        content: &str,
        skip_ms1: bool,
        sn: Option<u8>,
        sizes: Sizes,
    ) -> Result<Vec<RawSpectrum>, &'static str> {
        let path = temp(name, content);
        let out = read_processed(&path, 3, sn, skip_ms1, sizes, &|s| s);
        std::fs::remove_file(&path).ok();
        out
    }

    fn same(a: &[RawSpectrum], b: &[RawSpectrum]) -> bool {
        format!("{a:?}") == format!("{b:?}")
    }

    #[test]
    fn equals_the_serial_parse_for_any_piece_and_chunk_size() {
        for sep in ["", "\n      ", "\r\n\t"] {
            let content = synthetic(60, sep);
            for skip_ms1 in [false, true] {
                let expected = serial(&content, skip_ms1, None).unwrap();
                assert_eq!(expected.len(), if skip_ms1 { 48 } else { 60 });
                for (piece, chunk) in [
                    (4096, Some(1)),
                    (777, Some(20_000)),
                    (1 << 20, None),
                    (100, Some(1 << 30)),
                ] {
                    let sizes = Sizes { piece, chunk };
                    let got = parallel("eq", &content, skip_ms1, None, sizes).unwrap();
                    assert!(
                        same(&expected, &got),
                        "sep {sep:?} piece {piece} chunk {chunk:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn start_tags_with_other_whitespace_are_cut_too() {
        let content = synthetic(10, "\n").replace("<spectrum index=", "<spectrum\n        index=");
        let expected = serial(&content, false, None).unwrap();
        let sizes = Sizes {
            piece: 1000,
            chunk: Some(1),
        };
        assert!(same(
            &expected,
            &parallel("ws", &content, false, None, sizes).unwrap()
        ));
    }

    #[test]
    fn param_groups_and_noise_arrays_reach_every_chunk() {
        // MS level and centroid flag of every spectrum through a param group of the header, and a
        // noise array (S/N) in every spectrum
        let ms_level = r#"<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>"#;
        let centroid =
            r#"<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>"#;
        let base = synthetic(12, "\n")
            .replace(r#"name="ms level" value="1""#, r#"name="ms level" value="2""#)
            .replace(ms_level, r#"<referenceableParamGroupRef ref="ms2"/>"#)
            .replace(centroid, "")
            .replacen(
                r#"<referenceableParamGroupList count="1">"#,
                &format!(
                    r#"<referenceableParamGroupList count="2">
                    <referenceableParamGroup id="ms2">{ms_level}{centroid}</referenceableParamGroup>"#
                ),
                1,
            );
        // the intensity array a second time as a noise array
        let i = base.find(r#"accession="MS:1000515""#).unwrap();
        let start = base[..i].rfind("<binaryDataArray ").unwrap();
        let end = base[i..].find("</binaryDataArray>").unwrap() + i + "</binaryDataArray>".len();
        let noise = base[start..end].replace("MS:1000515", "MS:1002744");
        let content = base.replace(
            "</binaryDataArrayList>",
            &format!("{noise}</binaryDataArrayList>"),
        );
        let expected = serial(&content, true, Some(2)).unwrap();
        assert_eq!(expected.len(), 12);
        assert_eq!(expected[0].ms_level, 2);
        let sizes = Sizes {
            piece: 2048,
            chunk: Some(1),
        };
        assert!(same(
            &expected,
            &parallel("pg", &content, true, Some(2), sizes).unwrap()
        ));
    }

    #[test]
    fn a_spectrum_larger_than_a_piece_is_fine() {
        let content = synthetic(8, "\n");
        let expected = serial(&content, false, None).unwrap();
        let sizes = Sizes {
            piece: 64,
            chunk: Some(1),
        };
        assert!(same(
            &expected,
            &parallel("small", &content, false, None, sizes).unwrap()
        ));
    }

    #[test]
    fn unusual_files_go_to_the_serial_reader() {
        let sizes = Sizes {
            piece: 4096,
            chunk: Some(1),
        };
        let content = synthetic(10, "\n");
        let comment = content.replacen("</spectrum>", "</spectrum><!-- <spectrum id=\"x\"> -->", 3);
        let cdata = content.replacen("<precursorList", "<![CDATA[x]]><precursorList", 3);
        let pi = content.replacen("</spectrum>", "</spectrum><?pi x?>", 1);
        let not_after_end = content.replacen("</spectrum>", "</spectrum>text", 2);
        let no_list = content.replace(
            "<spectrumList count=\"10\" ",
            "<spectrumList2 count=\"10\" ",
        );
        let one = synthetic(1, "");
        let header_tag = content.replacen("<run ", "<!-- <spectrum id=\"x\"> --><run ", 1);
        for (name, c) in [
            ("comment", &comment),
            ("cdata", &cdata),
            ("pi", &pi),
            ("text", &not_after_end),
            ("list", &no_list),
            ("one", &one),
            ("header_tag", &header_tag),
        ] {
            assert!(parallel(name, c, false, None, sizes).is_err(), "{name}");
        }
        // a comment in the header is fine
        let header_comment = content.replacen("<run ", "<!-- converted for a test --><run ", 1);
        assert!(same(
            &serial(&header_comment, false, None).unwrap(),
            &parallel("hc", &header_comment, false, None, sizes).unwrap()
        ));
    }

    #[test]
    fn errors_go_to_the_serial_reader() {
        let sizes = Sizes {
            piece: 4096,
            chunk: Some(1),
        };
        let content = synthetic(10, "\n");
        // truncated inside a spectrum, after the spectrum list, and before </mzML>
        for cut in [
            content.len() / 2,
            content.find("</spectrumList>").unwrap(),
            content.find("</mzML>").unwrap(),
        ] {
            assert!(serial(&content[..cut], false, None).is_err());
            assert!(
                parallel("cut", &content[..cut], false, None, sizes).is_err(),
                "cut {cut}"
            );
        }
        let numpress = content.replace("MS:1000574", "MS:1002312");
        assert!(serial(&numpress, false, None).is_err());
        assert!(parallel("numpress", &numpress, false, None, sizes).is_err());
    }

    #[test]
    fn read_processed_takes_the_parallel_path_and_matches_the_serial_one() {
        let content = synthetic(40, "\n");
        let path = temp("util", &content);
        let url = crate::to_url(path.to_str().unwrap()).unwrap();
        let f = |s: RawSpectrum| (s.id.clone(), s.mz.len(), s.ms_level);
        assert!(read_processed(&path, 0, None, true, Sizes::default(), &f).is_ok());
        let parallel =
            crate::util::read_processed(&url, 0, None, Default::default(), false, f).unwrap();
        let serial: Vec<_> = crate::util::read_mzml_levels(&url, 0, None, true)
            .unwrap()
            .into_iter()
            .map(f)
            .collect();
        std::fs::remove_file(&path).ok();
        assert_eq!(parallel.len(), 32);
        assert_eq!(parallel, serial);
    }
}
