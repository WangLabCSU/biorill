//! Bioinformatics file formats and I/O.
//!
//! Enable the `fastx` feature for FASTQ records or the `kraken` feature for
//! Kraken reports. No features are enabled by default.
//!
//! The format crates, `biorill-fastx` and `biorill-kraken`, can also be used
//! independently. The `fastx` feature currently requires Rust nightly.

#[cfg(feature = "fastx")]
/// FASTQ records, reading, and writing.
///
/// ```
/// use biorill::fastx::FastqReader;
///
/// let mut reader = FastqReader::new(b"@read1\nAC\n+\n12\n".as_slice());
/// let record = reader.read_view()?.unwrap();
/// assert_eq!(record.sequence(), b"AC");
/// # Ok::<(), biorill::fastx::Error>(())
/// ```
pub use biorill_fastx as fastx;

#[cfg(feature = "kraken")]
/// Kraken reports and taxonomic selection.
///
/// ```
/// use biorill::kraken::KrakenReportReader;
///
/// let input = b"100\t2\t2\tD\t2\tBacteria\n";
/// let mut reader = KrakenReportReader::new(input.as_slice());
/// let entry = reader.read_entry()?.unwrap();
/// assert_eq!(entry.taxon().term(), "Bacteria");
/// # Ok::<(), biorill::kraken::Error>(())
/// ```
pub use biorill_kraken as kraken;
