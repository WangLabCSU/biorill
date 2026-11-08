use std::{mem, ops::Range};

use leftover::LeftoverView;
use memchr::{memchr, memchr2};
use thiserror::Error;

use super::vo::{FastqParts, FastqRecordError};

/// Locate and validate FASTQ records in buffered input.
///
/// The caller keeps each input's progress in `ReadingState`. Parsing borrows
/// the view without copying or consuming its bytes; the caller decides when
/// to retain a chunk, discard a blank line, or take a completed record.
#[derive(Default)]
pub(crate) struct FastqParser;

/// Progress saved by the reader between chunks and calls.
pub(crate) struct ReadingState {
    pub(crate) state: RecordState,
    pub(crate) offset: usize,
}

impl ReadingState {
    pub(crate) fn new() -> Self {
        Self {
            state: RecordState::Empty,
            offset: 0,
        }
    }
}

/// The physical record lines located so far, excluding line endings.
pub(crate) enum RecordState {
    Empty,
    Header {
        header: Range<usize>,
    },
    Sequence {
        header: Range<usize>,
        sequence: Range<usize>,
    },
    Separator {
        header: Range<usize>,
        sequence: Range<usize>,
        separator: Range<usize>,
    },
}

pub(crate) enum ParsingStep {
    /// Discard the leftover and this many bytes of the current chunk.
    Skip { consumed: usize },
    /// Retain the current chunk, then supply another one with the same state.
    NeedMore,
    /// Take the leftover and `consumed` bytes of the current chunk.
    /// Field ranges exclude format markers and line endings.
    Complete {
        parts: FastqParts<Range<usize>>,
        consumed: usize,
    },
    /// Input ended with no record left to return. Keep the view's bytes.
    Eof,
}

impl FastqParser {
    /// Create a parser that can be used with independent input states.
    pub(crate) fn new() -> Self {
        Self
    }

    /// Scan buffered lines until a record completes or input is needed.
    ///
    /// Searches the current chunk for LF, strips CR only immediately before
    /// LF, and validates each complete line. The leftover contains bytes
    /// already searched for LF; retain the whole chunk after `NeedMore` before
    /// supplying the next view. Saved ranges address leftover + current chunk.
    /// `start` is the current line's start in those bytes. Keep it across chunks
    /// while reading one record; begin a new read with zero.
    /// `state` belongs to this input, so the same parser can process other
    /// inputs using their own states.
    ///
    /// Set `eof` when no more input is available. A nonempty remaining tail is
    /// then a final line, and an empty tail contributes no line. This allows
    /// a final nonempty quality line without LF. EOF itself consumes no bytes.
    /// The line count includes skipped blank lines and an invalid final line.
    ///
    /// # Errors
    /// Returns an error for invalid markers, unequal sequence/quality lengths,
    /// or an incomplete record at EOF. Stop processing this input after an error.
    ///
    /// # Panics
    /// Panics if a saved range does not address the supplied view.
    pub(crate) fn parse(
        &self,
        view: &LeftoverView<'_>,
        state: &mut ReadingState,
        start: &mut usize,
        eof: bool,
    ) -> Result<ParsingStep, ParseError> {
        let prefix = view.leftover().len();
        let chunk = view.buffer();
        // Earlier chunks have already been searched. Unlike start,
        // scanned is relative to this chunk and only lives for this call.
        let mut scanned = 0;

        loop {
            // SAFETY: scanned starts at zero. Its only continuing update adds
            // pos + 1 after memchr finds a byte in chunk[scanned..], so it stays
            // at most chunk.len(). The borrowed view cannot change this chunk.
            let buf = unsafe { chunk.get_unchecked(scanned..) };
            let (end, terminated) = if let Some(pos) = memchr(b'\n', buf) {
                let mut end = prefix + scanned + pos;
                // Only CR immediately before this LF is part of a line ending.
                // If pos == 0 and scanned > 0, the preceding byte is the LF
                // found in the previous iteration: this is an empty line.
                if pos > 0 {
                    // SAFETY: memchr returned pos < buf.len(), and pos > 0
                    // places the preceding byte inside this same slice.
                    if *unsafe { buf.get_unchecked(pos - 1) } == b'\r' {
                        end -= 1;
                    }
                } else if scanned == 0 && *start < prefix {
                    // LF begins this chunk and the current line has bytes in
                    // leftover, so CRLF may straddle the two slices. If start
                    // equals prefix, this line is empty; an earlier CR belongs
                    // to a previous line and must not be trimmed.
                    // SAFETY: start < prefix implies prefix > 0. prefix is the
                    // length of the unchanged leftover slice.
                    if *unsafe { view.leftover().get_unchecked(prefix - 1) } == b'\r' {
                        end -= 1;
                    }
                }
                // Trimming changes the field end only. The next line starts
                // after LF, including both bytes when the ending is CRLF.
                scanned += pos + 1;
                (end, true)
            } else if eof {
                let end = prefix + chunk.len();
                // Reader signals EOF with an empty chunk: only the retained
                // tail can contribute a final line. Also accept explicit EOF
                // with a final nonempty chunk, without inferring it from size.
                // No remaining bytes means no line, particularly no empty
                // quality line. Keep the saved state for incomplete input.
                if *start == end {
                    self.parse_eof(view, &state.state)?;
                    return Ok(ParsingStep::Eof);
                }
                scanned = chunk.len();
                // There is no LF, so a trailing bare CR stays in the field.
                (end, false)
            } else {
                return Ok(ParsingStep::NeedMore);
            };

            state.offset += 1;
            let range = *start..end;
            let line @ (line_leftover, line_chunk) = view.slice(range.clone());
            if matches!(state.state, RecordState::Empty)
                && line_leftover
                    .iter()
                    .chain(line_chunk)
                    .all(u8::is_ascii_whitespace)
            {
                if terminated {
                    *start = 0;
                    return Ok(ParsingStep::Skip { consumed: scanned });
                }
                // Keep final whitespace. A later read starts over and parses
                // this tail again, including counting it as another line.
                *start = end;
                return Ok(ParsingStep::Eof);
            }

            if let Some(parts) = self.parse_line(view, range, line, &mut state.state)? {
                *start = 0;
                return Ok(ParsingStep::Complete {
                    parts,
                    consumed: scanned,
                });
            }
            *start = prefix + scanned;
            if !terminated {
                // The nonempty tail did not complete a record. Report the
                // missing lines without discarding the saved record state.
                self.parse_eof(view, &state.state)?;
                return Ok(ParsingStep::Eof);
            }
        }
    }

    // Framing supplies one line without its ending. This method only advances
    // FASTQ record state; it neither searches for LF nor consumes buffer bytes.
    fn parse_line(
        &self,
        view: &LeftoverView<'_>,
        range: Range<usize>,
        line: (&[u8], &[u8]),
        state: &mut RecordState,
    ) -> Result<Option<FastqParts<Range<usize>>>, ParseError> {
        *state = match mem::replace(state, RecordState::Empty) {
            RecordState::Empty => {
                self.check_header(line)?;
                RecordState::Header { header: range }
            }
            RecordState::Header { header } => RecordState::Sequence {
                header,
                sequence: range,
            },
            RecordState::Sequence { header, sequence } => {
                self.check_separator(
                    view.slice(header.clone()),
                    view.slice(sequence.clone()),
                    line,
                )?;
                RecordState::Separator {
                    header,
                    sequence,
                    separator: range,
                }
            }
            RecordState::Separator {
                header,
                sequence,
                separator,
            } => {
                self.check_lengths(
                    view.slice(header.clone()),
                    view.slice(sequence.clone()),
                    view.slice(separator.clone()),
                    line,
                )?;
                let (id, description) = header_ranges(view.slice(header.clone()));
                let id = header.start + id.start..header.start + id.end;
                let description =
                    description.map(|range| header.start + range.start..header.start + range.end);
                let separator = separator.start + 1..separator.end;
                return Ok(Some(FastqParts::new(
                    id,
                    description,
                    sequence,
                    separator,
                    range,
                )));
            }
        };
        Ok(None)
    }

    /// Check that input ended between records, leaving `state` unchanged.
    /// The view must contain the preceding lines saved in `state`.
    ///
    /// # Errors
    /// Returns `IncompleteRecord` if any record lines are still pending.
    ///
    /// # Panics
    /// Panics if a saved line range is outside the view. May panic if a
    /// previously accepted header has changed.
    fn parse_eof(&self, view: &LeftoverView<'_>, state: &RecordState) -> Result<(), ParseError> {
        let source = match state {
            RecordState::Empty => return Ok(()),
            RecordState::Header { header } => {
                self.incomplete_record(view, header.clone(), None, None)
            }
            RecordState::Sequence { header, sequence } => {
                self.incomplete_record(view, header.clone(), Some(sequence.clone()), None)
            }
            RecordState::Separator {
                header,
                sequence,
                separator,
            } => self.incomplete_record(
                view,
                header.clone(),
                Some(sequence.clone()),
                Some(separator.clone()),
            ),
        };
        Err(source)
    }

    // These checks take constant time on valid input and report invalid lines
    // while framing, before another read can fail.
    /// Check the header's leading `@` marker.
    ///
    /// # Errors
    /// Returns `InvalidHead` if the marker is missing.
    fn check_header(&self, header: (&[u8], &[u8])) -> Result<(), ParseError> {
        let (leftover, chunk) = header;
        if leftover.first().or_else(|| chunk.first()) != Some(&b'@') {
            return Err(ParseError::InvalidHead {
                record: display_lines(&[header]),
            });
        }
        Ok(())
    }

    /// Check the separator's leading `+` marker.
    /// The preceding lines are included in the error message if needed.
    ///
    /// # Errors
    /// Returns `InvalidSep` if the marker is missing.
    fn check_separator(
        &self,
        header: (&[u8], &[u8]),
        sequence: (&[u8], &[u8]),
        separator: (&[u8], &[u8]),
    ) -> Result<(), ParseError> {
        let (leftover, chunk) = separator;
        if leftover.first().or_else(|| chunk.first()) != Some(&b'+') {
            return Err(ParseError::InvalidSep {
                record: display_lines(&[header, sequence, separator]),
            });
        }
        Ok(())
    }

    /// Check that sequence and quality have the same byte length.
    /// The header and separator are included in the error message if needed.
    ///
    /// # Errors
    /// Returns `UnequalLength` if the lengths differ.
    fn check_lengths(
        &self,
        header: (&[u8], &[u8]),
        sequence: (&[u8], &[u8]),
        separator: (&[u8], &[u8]),
        quality: (&[u8], &[u8]),
    ) -> Result<(), ParseError> {
        let (sequence_leftover, sequence_chunk) = sequence;
        let (quality_leftover, quality_chunk) = quality;
        let seq = sequence_leftover.len() + sequence_chunk.len();
        let qual = quality_leftover.len() + quality_chunk.len();
        if seq != qual {
            return Err(ParseError::UnequalLength {
                record: display_lines(&[header, sequence, separator, quality]),
                seq,
                qual,
            });
        }
        Ok(())
    }

    /// Format an incomplete-record error from the available record lines.
    /// Supply a header accepted by [`Self::check_header`], and include a
    /// separator only when a sequence line is also present. Ranges address
    /// the complete view.
    ///
    /// # Panics
    /// Panics if a range is outside the view. May panic if the supplied header
    /// has not passed the header check.
    fn incomplete_record(
        &self,
        view: &LeftoverView<'_>,
        header: Range<usize>,
        sequence: Option<Range<usize>>,
        separator: Option<Range<usize>>,
    ) -> ParseError {
        // Keep the existing diagnostic text for incomplete input.
        let record = match (sequence, separator) {
            (None, _) => {
                let (id, description) = header_ranges(view.slice(header.clone()));
                format!(
                    "{}{}",
                    display_lines(&[view.slice(header.start + id.start..header.start + id.end)]),
                    display_lines(&[description
                        .map(|range| {
                            view.slice(header.start + range.start..header.start + range.end)
                        })
                        .unwrap_or_default()]),
                )
            }
            (Some(sequence), None) => display_lines(&[view.slice(header), view.slice(sequence)]),
            (Some(sequence), Some(separator)) => {
                let (_, description) = header_ranges(view.slice(header.clone()));
                display_lines(&[
                    description
                        .map(|range| {
                            view.slice(header.start + range.start..header.start + range.end)
                        })
                        .unwrap_or_default(),
                    view.slice(sequence),
                    view.slice(separator),
                ])
            }
        };
        ParseError::IncompleteRecord { record }
    }
}

fn header_ranges((leftover, chunk): (&[u8], &[u8])) -> (Range<usize>, Option<Range<usize>>) {
    let len = leftover.len() + chunk.len();
    let split = memchr2(b' ', b'\t', leftover)
        .or_else(|| memchr2(b' ', b'\t', chunk).map(|at| leftover.len() + at));
    let id_end = split.unwrap_or(len);
    let description = split.filter(|&at| at + 1 < len).map(|at| at + 1..len);
    (1..id_end, description)
}

fn display_lines(lines: &[(&[u8], &[u8])]) -> String {
    // Only error diagnostics allocate here. Join bytes before decoding so a
    // UTF-8 character split across chunks is displayed as one character.
    let mut bytes = Vec::new();
    for (i, &(leftover, chunk)) in lines.iter().enumerate() {
        if i != 0 {
            bytes.push(b'\n');
        }
        bytes.extend_from_slice(leftover);
        bytes.extend_from_slice(chunk);
    }
    String::from_utf8_lossy(&bytes).into_owned()
}

/// A four-line FASTQ record could not be parsed.
#[derive(Debug, Error)]
pub enum ParseError {
    /// The first line does not start with `@`.
    #[error("expected '@' at record start\n{record}")]
    InvalidHead { record: String },
    /// The sequence and quality byte lengths differ.
    #[error("sequence and quality lengths do not match ({seq} vs {qual})\n{record}")]
    UnequalLength {
        record: String,
        seq: usize,
        qual: usize,
    },
    /// The third line does not start with `+`.
    #[error("expected '+' at separator line start\n{record}")]
    InvalidSep { record: String },
    /// The input ended before all four record lines were supplied.
    #[error("incomplete record\n{record}")]
    IncompleteRecord { record: String },
    /// The completed buffer violates the record's field or range constraints.
    #[error("invalid FASTQ record: {source}")]
    InvalidRecord {
        #[source]
        source: FastqRecordError,
    },
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use leftover::LeftoverReader;

    use super::{FastqParser, ParseError, ParsingStep, ReadingState, RecordState};
    use crate::{FastqFieldError, FastqParts, FastqRecord, FastqRecordError};

    #[test]
    fn parse_stops_at_the_first_record_without_consuming_the_view() {
        let first = b"@r note\nA\xff\n+r\n12\n";
        let second = b"@next\nT\n+\n3\n";
        let input = [first.as_slice(), second].concat();
        let parser = FastqParser::new();
        let mut state = ReadingState::new();
        let mut start = 0;
        let mut reader = LeftoverReader::new(input.as_slice());
        let view = reader.fill_view().unwrap();
        let ptr = view.buffer().as_ptr();

        let ParsingStep::Complete { parts, consumed } =
            parser.parse(&view, &mut state, &mut start, false).unwrap()
        else {
            panic!("expected the first complete record");
        };
        assert_eq!(consumed, first.len());
        assert_eq!(state.offset, 4);
        assert!(matches!(state.state, RecordState::Empty));
        assert_eq!(view.buffer(), input);
        assert_eq!(
            parts,
            FastqParts::new(1..2, Some(3..7), 8..10, 12..13, 14..16)
        );
        let record = FastqRecord::try_from((view.take_view(consumed), parts)).unwrap();
        assert_eq!(record.id(), b"r");
        assert_eq!(record.description(), Some(b"note".as_slice()));
        assert_eq!(record.sequence(), b"A\xff");
        assert_eq!(record.separator(), b"r");
        assert_eq!(record.quality(), b"12");
        let (raw, _) = record.into_parts();
        assert_eq!(raw.as_ptr(), ptr);
        assert_eq!(raw, first);

        let view = reader.fill_view().unwrap();
        assert_eq!(view.buffer(), second);
        assert!(matches!(
            parser.parse(&view, &mut state, &mut start, false),
            Ok(ParsingStep::Complete { .. })
        ));
        assert_eq!(state.offset, 8);
    }

    #[test]
    fn record_ranges_are_independent_of_every_buffer_split() {
        let parser = FastqParser::new();
        for input in [
            b"@r note\nA\xff\n+r\n12\n".as_slice(),
            b"@r note\r\nA\xff\r\n+r\r\n12\r\n",
            b"@r note\r\nA\xff\r\n+r\r\n12",
        ] {
            for split in 0..=input.len() {
                let (first, second) = input.split_at(split);
                let mut reader = LeftoverReader::new(first.chain(second));
                let mut state = ReadingState::new();
                let mut start = 0;
                let record = loop {
                    let view = reader.fill_view().unwrap();
                    let eof = view.buffer().is_empty();
                    match parser.parse(&view, &mut state, &mut start, eof).unwrap() {
                        ParsingStep::NeedMore => view.accumulate(),
                        ParsingStep::Complete { parts, consumed } => {
                            break FastqRecord::try_from((view.take_owned(consumed), parts))
                                .unwrap();
                        }
                        _ => panic!("expected a record at split {split}"),
                    }
                };
                assert_eq!(record.id(), b"r");
                assert_eq!(record.description(), Some(b"note".as_slice()));
                assert_eq!(record.sequence(), b"A\xff");
                assert_eq!(record.separator(), b"r");
                assert_eq!(record.quality(), b"12");
                assert_eq!(record.into_parts().0, input);
                assert_eq!(state.offset, 4);
            }
        }
    }

    #[test]
    fn explicit_eof_can_complete_a_nonempty_chunk_without_stripping_bare_cr() {
        let parser = FastqParser::new();
        let input = b"@r\nAC\n+\n1\r";
        let mut reader = LeftoverReader::new(input.as_slice());
        let mut state = ReadingState::new();
        let mut start = 0;
        let view = reader.fill_view().unwrap();
        let ParsingStep::Complete { parts, consumed } =
            parser.parse(&view, &mut state, &mut start, true).unwrap()
        else {
            panic!("the explicit EOF must complete the final quality line");
        };
        assert_eq!(consumed, input.len());
        assert_eq!(state.offset, 4);
        let record = FastqRecord::try_from((view.take_view(consumed), parts)).unwrap();
        assert_eq!(record.quality(), b"1\r");
    }

    #[test]
    fn eof_is_explicit_even_when_the_chunk_is_empty() {
        let parser = FastqParser::new();
        let mut reader = LeftoverReader::new(b"".as_slice());
        let view = reader.fill_view().unwrap();
        let mut state = ReadingState::new();
        let mut start = 0;
        assert!(matches!(
            parser.parse(&view, &mut state, &mut start, false),
            Ok(ParsingStep::NeedMore)
        ));
        assert!(matches!(
            parser.parse(&view, &mut state, &mut start, true),
            Ok(ParsingStep::Eof)
        ));
        assert_eq!(state.offset, 0);
    }

    #[test]
    fn eof_preserves_incomplete_state_and_does_not_invent_an_empty_quality_line() {
        let parser = FastqParser::new();
        for (input, lines) in [
            (b"@r\n".as_slice(), 1),
            (b"@r\nAC\n".as_slice(), 2),
            (b"@r\n\n+\n".as_slice(), 3),
        ] {
            let mut reader = LeftoverReader::new(input);
            let mut state = ReadingState::new();
            let mut start = 0;
            let view = reader.fill_view().unwrap();
            assert!(matches!(
                parser.parse(&view, &mut state, &mut start, false),
                Ok(ParsingStep::NeedMore)
            ));
            view.accumulate();
            let view = reader.fill_view().unwrap();
            for _ in 0..2 {
                assert!(matches!(
                    parser.parse(&view, &mut state, &mut start, true),
                    Err(ParseError::IncompleteRecord { .. })
                ));
                assert_eq!(state.offset, lines);
                assert_eq!(start, input.len());
                assert!(matches!(
                    (&state.state, lines),
                    (RecordState::Header { .. }, 1)
                        | (RecordState::Sequence { .. }, 2)
                        | (RecordState::Separator { .. }, 3)
                ));
                assert_eq!(view.leftover(), input);
            }
        }
    }

    #[test]
    fn one_parser_can_advance_independent_input_states() {
        let parser = FastqParser::new();
        let first = b"@a\nAC\n+\n12\n";
        let second = b"@b\nT\n+b\n3\n";
        let mut first_reader = LeftoverReader::new(first[..4].chain(&first[4..]));
        let mut second_reader = LeftoverReader::new(second[..7].chain(&second[7..]));
        let mut first_state = ReadingState::new();
        let mut first_start = 0;
        let mut second_state = ReadingState::new();
        let mut second_start = 0;
        for (reader, state, start) in [
            (&mut first_reader, &mut first_state, &mut first_start),
            (&mut second_reader, &mut second_state, &mut second_start),
        ] {
            let view = reader.fill_view().unwrap();
            assert!(matches!(
                parser.parse(&view, state, start, false),
                Ok(ParsingStep::NeedMore)
            ));
            view.accumulate();
        }
        assert_eq!(first_state.offset, 1);
        assert_eq!(second_state.offset, 2);
        // Resume the two streams in the opposite order using one parser.
        for (reader, state, start, id) in [
            (
                &mut second_reader,
                &mut second_state,
                &mut second_start,
                b"b",
            ),
            (&mut first_reader, &mut first_state, &mut first_start, b"a"),
        ] {
            let view = reader.fill_view().unwrap();
            let ParsingStep::Complete { parts, consumed } =
                parser.parse(&view, state, start, false).unwrap()
            else {
                panic!("expected a complete record");
            };
            let record = FastqRecord::try_from((view.take_view(consumed), parts)).unwrap();
            assert_eq!(record.id(), id);
            assert_eq!(state.offset, 4);
        }
    }

    #[test]
    fn whitespace_is_skipped_only_between_records() {
        let parser = FastqParser::new();
        let input = b" \t\r\n@r\n \t\n+\n \t\n \r";
        let mut reader = LeftoverReader::new(input.as_slice());
        let mut state = ReadingState::new();
        let mut start = 0;
        let view = reader.fill_view().unwrap();
        let ParsingStep::Skip { consumed } =
            parser.parse(&view, &mut state, &mut start, false).unwrap()
        else {
            panic!("expected leading whitespace to be skipped");
        };
        assert_eq!(consumed, 4);
        assert_eq!(view.buffer(), input);
        assert_eq!(state.offset, 1);
        view.consume(consumed);
        let view = reader.fill_view().unwrap();
        let ParsingStep::Complete { parts, consumed } =
            parser.parse(&view, &mut state, &mut start, false).unwrap()
        else {
            panic!("expected whitespace sequence and quality fields");
        };
        let record = FastqRecord::try_from((view.take_view(consumed), parts)).unwrap();
        assert_eq!(record.sequence(), b" \t");
        assert_eq!(record.quality(), b" \t");
        assert_eq!(state.offset, 5);
        let view = reader.fill_view().unwrap();
        assert!(matches!(
            parser.parse(&view, &mut state, &mut start, false),
            Ok(ParsingStep::NeedMore)
        ));
        view.accumulate();
        let view = reader.fill_view().unwrap();
        for lines in 6..=7 {
            // Each new read begins at the start of the retained tail.
            start = 0;
            assert!(matches!(
                parser.parse(&view, &mut state, &mut start, true),
                Ok(ParsingStep::Eof)
            ));
            assert_eq!(state.offset, lines);
            assert_eq!(view.leftover(), b" \r");
        }
    }

    #[test]
    fn errors_count_the_offending_line_across_buffer_splits() {
        let parser = FastqParser::new();
        for (input, expected_line) in [
            (b"#r\n".as_slice(), 1),
            (b"@r\nAC\n-\n".as_slice(), 3),
            (b"@r\nA\n+\n12\n".as_slice(), 4),
        ] {
            for split in 0..=input.len() {
                let (first, second) = input.split_at(split);
                let mut reader = LeftoverReader::new(first.chain(second));
                let mut state = ReadingState::new();
                let mut start = 0;
                let error = loop {
                    let view = reader.fill_view().unwrap();
                    match parser.parse(&view, &mut state, &mut start, view.buffer().is_empty()) {
                        Ok(ParsingStep::NeedMore) => view.accumulate(),
                        Err(error) => break error,
                        _ => panic!("invalid input must fail at split {split}"),
                    }
                };
                assert_eq!(state.offset, expected_line);
                assert!(matches!(
                    (error, expected_line),
                    (ParseError::InvalidHead { .. }, 1)
                        | (ParseError::InvalidSep { .. }, 3)
                        | (
                            ParseError::UnequalLength {
                                seq: 1,
                                qual: 2,
                                ..
                            },
                            4
                        )
                ));
            }
        }
    }

    #[test]
    fn record_construction_checks_field_ranges_and_contents() {
        let parts = FastqParts::new(1..2, None, 3..5, 7..7, 8..10);
        for input in [b"".as_slice(), b"@r\nAC\n+\n"] {
            let FastqRecordError::InvalidRange { field } =
                FastqRecord::try_from((input, parts.clone())).unwrap_err()
            else {
                panic!("expected an out-of-bounds field range");
            };
            assert_eq!(
                field,
                if input.is_empty() {
                    "identifier"
                } else {
                    "quality"
                }
            );
        }
        assert!(matches!(
            FastqRecord::try_from((b"@r\nA\n\n+\n12\n".as_slice(), parts)),
            Err(FastqRecordError::InvalidFields(FastqFieldError::LineFeed {
                field: "sequence"
            }))
        ));
    }
}
