//! Sequencing record models and buffered FASTQ reading.
//!
//! [`FastqReader`] handles LF and CRLF,
//! skips whitespace-only lines between records, and checks the header, separator
//! and sequence/quality lengths. It preserves sequence bytes, quality bytes and
//! separator text; it does not validate alphabets or quality encodings.
//! Wrapped sequences and qualities are not supported.
//! Separator fields contain the content after the leading `+`, which is empty
//! for a bare `+` line. Writing adds the marker before this content.
//!
//! [`FastqField`] and [`FastqRecord`] accept storage implementing
//! [`stringbytes::StringBytes`], including byte and text buffers, references,
//! and shared pointers. Custom storage can lend a byte slice for validation.
//!
//! Use [`FastqReader::read_view`] to borrow records. Records spanning buffer
//! refills are copied into reusable storage. [`FastqReader::read_owned`]
//! produces records that can outlive the reader by
//! copying from the input buffer and any collected record bytes.
//! All modes retain original line endings and header whitespace.
//! [`FastqRecord::write`] and [`FastqField::write`] write normalized FASTQ to any
//! [`std::io::Write`], using LF and a space before descriptions. Use
//! `encode_into_vec` or `encode_into_buf` to append encoded fields to a buffer.
//! Decompression, file handling and read selection belong to the caller.
//!
//! ```
//! use biorill_fastx::FastqReader;
//!
//! let input = b"@read1\tUMI:AC\r\nAC\r\n+\r\nII\r\n";
//! let mut reader = FastqReader::new(input.as_slice());
//! let record = reader.read_owned()?.unwrap();
//! assert_eq!(record.id(), b"read1");
//! assert_eq!(record.sequence(), b"AC");
//! assert_eq!(record.separator(), b"");
//! let mut encoded = Vec::new();
//! record.encode_into_vec(&mut encoded);
//! assert_eq!(encoded, b"@read1 UMI:AC\nAC\n+\nII\n");
//! assert_eq!(record.into_parts().0, input);
//! assert!(reader.read_owned()?.is_none());
//! # Ok::<(), biorill_fastx::Error>(())
//! ```
#![deny(unnameable_types, unreachable_pub)]

mod domain;
mod error;
mod reader;

pub use self::domain::{
    FastqField, FastqFieldError, FastqPairError, FastqParts, FastqRecord, FastqRecordError,
    LengthMismatch, ParseError,
};
pub use self::error::Error;
pub use self::reader::FastqReader;
