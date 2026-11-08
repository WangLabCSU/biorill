use std::error::Error;
use std::fmt;

/// The identifiers of two records expected to form a read pair differ.
#[derive(Debug)]
pub struct FastqPairError {
    /// Identifier of the first read.
    pub read1_id: String,
    /// Identifier of the second read.
    pub read2_id: String,
    /// Optional physical line number for the first read.
    pub read1_pos: Option<usize>,
    /// Optional physical line number for the second read.
    pub read2_pos: Option<usize>,
}

impl fmt::Display for FastqPairError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "FASTQ pairing error: sequence IDs do not match\n  record1 ID"
        )?;
        if let Some(pos) = self.read1_pos {
            write!(f, " (line: {pos})")?;
        }
        write!(f, ": {}\n  record2 ID", self.read1_id)?;
        if let Some(pos) = self.read2_pos {
            write!(f, " (line: {pos})")?;
        }
        write!(f, ": {}", self.read2_id)
    }
}

impl Error for FastqPairError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_pair_error() {
        let err = FastqPairError {
            read1_id: "SEQ1".into(),
            read2_id: "SEQ2".into(),
            read1_pos: Some(3),
            read2_pos: Some(3),
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("sequence IDs do not match"),
            "Should mention pairing mismatch"
        );
        assert!(msg.contains("(line: 3): SEQ1"));
        assert!(msg.contains("(line: 3): SEQ2"));
    }
}
