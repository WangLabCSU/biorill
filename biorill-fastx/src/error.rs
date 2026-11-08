use std::io;

use thiserror::Error;

use crate::domain::ParseError;

/// A FASTQ source could not be read or a record could not be parsed.
/// Parse errors include the physical line number. I/O errors do not have a
/// line number; inspect [`crate::FastqReader::offset`] for input progress.
/// Inspect the underlying cause with `std::error::Error::source()`.
#[derive(Debug, Error)]
pub enum Error {
    /// The input failed while reading bytes for a record.
    #[error("Failed to read FASTQ: {source}")]
    Read {
        /// The underlying input failure.
        #[from]
        source: io::Error,
    },
    /// The domain parser rejected the record or found an incomplete record at EOF.
    #[error("FASTQ parse error (line: {line}): {source}")]
    Parse {
        /// The rejected line, or the last line supplied for an incomplete record.
        line: usize,
        /// The underlying record failure.
        #[source]
        source: ParseError,
    },
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use super::*;

    #[test]
    fn parse_errors_preserve_the_line_and_domain_cause() {
        let error = Error::Parse {
            line: 42,
            source: ParseError::InvalidHead {
                record: "SEQ_ID".into(),
            },
        };
        assert_eq!(
            error.to_string(),
            "FASTQ parse error (line: 42): expected '@' at record start\nSEQ_ID"
        );
        assert!(error.source().unwrap().is::<ParseError>());
        let error = Error::Parse {
            line: 100,
            source: ParseError::UnequalLength {
                record: "SEQ\nQUAL".into(),
                seq: 10,
                qual: 8,
            },
        };
        assert_eq!(error.to_string(), "FASTQ parse error (line: 100): sequence and quality lengths do not match (10 vs 8)\nSEQ\nQUAL");
    }

    #[test]
    fn read_errors_preserve_the_io_cause() {
        let error = Error::from(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated input",
        ));
        assert_eq!(error.to_string(), "Failed to read FASTQ: truncated input");
        assert_eq!(
            error
                .source()
                .unwrap()
                .downcast_ref::<io::Error>()
                .unwrap()
                .kind(),
            io::ErrorKind::UnexpectedEof
        );
    }
}
