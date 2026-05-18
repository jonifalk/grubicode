"""exercise_2_8819989_Falk.py"""
import argparse
import gzip
import itertools
import os
import textwrap

__author__ = "8819989, Falk"
__version__ = "1.0.0"


class FastaEntry:
    """A class to store and manipulate a single FASTA entry."""

    def __init__(self, header: str = "", sequence: str = ""):
        """
        Initialize the FastaEntry object.

        Args:
            header (str): The FASTA header line (starts with ">").
            sequence (str): The continuous nucleotide sequence.
        """
        self.header = header
        self.sequence = sequence

    def reverse_complement(self) -> str:
        """
        Converts the sequence to its reverse complement in-place.

        Returns:
            str: The newly computed reverse complement sequence.
        """
        complement_map = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
        self.sequence = self.sequence[::-1].translate(complement_map)
        return self.sequence

    def hamming_distance(self, other: 'FastaEntry') -> float:
        """
        Calculates the normalized Hamming distance between this sequence
        and another FastaEntry's sequence.

        If sequences are of different lengths, the length difference is
        counted as mismatches, and the maximum length is used as the
        denominator.

        Args:
            other (FastaEntry): The other FastaEntry object to compare
                                with.

        Returns:
            float: The normalized Hamming distance (0.0 to 1.0).
        """
        seq1 = self.sequence
        seq2 = other.sequence

        if not seq1 and not seq2:
            return 0.0

        max_len = max(len(seq1), len(seq2))

        # Base mismatch count starts with the difference in lengths
        mismatches = abs(len(seq1) - len(seq2))

        # Add mismatches from overlapping characters
        for c1, c2 in zip(seq1, seq2):
            if c1 != c2:
                mismatches += 1

        return mismatches / max_len

    def __str__(self) -> str:
        """
        Wraps the sequence to 60 characters per line to adhere to typical
        FASTA formatting guidelines.
        """
        wrapped_sequence = textwrap.fill(self.sequence, width=60)
        return f"{self.header}\n{wrapped_sequence}"


def is_gzipped(filepath: str) -> bool:
    """
    Checks if a file is gzip-compressed by reading its magic bytes.
    GZIP files always start with the magic number 1f 8b.

    Args:
        filepath (str): The path to the file.

    Returns:
        bool: True if the file is a gzip file, False otherwise.
    """
    if not os.path.exists(filepath):
        return False

    try:
        with open(filepath, 'rb') as f:
            magic_bytes = f.read(2)
            # \x1f\x8b are the standard magic bytes for a gzip file
            return magic_bytes == b'\x1f\x8b'
    except IOError:
        return False


def read_fasta(
        filepath: str,
        force_compressed: bool = False
) -> list[FastaEntry]:
    """
    Reads a Multi-FASTA file and parses it into a list of FastaEntry
    objects.

    Args:
        filepath (str): Path to the FASTA file.
        force_compressed (bool): Flag indicating if decompression should
                                 be forced.

    Returns:
        list: A list of parsed FastaEntry objects.
    """
    entries = []

    # Auto-detect if file is gzipped (Problem 4 bonus)
    is_compressed = is_gzipped(filepath) or force_compressed

    # Select appropriate file opener and mode based on compression state
    open_func = gzip.open if is_compressed else open
    mode = 'rt' if is_compressed else 'r'

    with open_func(filepath, mode) as f:
        current_header = ""
        current_sequence_lines = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # If we've already stored a header, append the
                # accumulated sequence to our list before starting a
                # new one.
                if current_header:
                    sequence = "".join(current_sequence_lines)
                    entries.append(FastaEntry(current_header, sequence))

                current_header = line
                current_sequence_lines = []
            else:
                # Accumulate nucleotide lines
                current_sequence_lines.append(line)

        # Make sure to append the final entry once the file is fully read
        if current_header:
            sequence = "".join(current_sequence_lines)
            entries.append(FastaEntry(current_header, sequence))

    return entries


def main(args):
    """
    Main execution logic:
    Reads a Multi-FASTA file, processes the entries, computes pairwise
    Hamming distances, and writes output to a TSV file.
    """
    # 1 & 3 & 4. Parse the Multi-FASTA file
    entries = read_fasta(args.input, force_compressed=args.compressed)

    if not entries:
        print("Error: Input file is empty or formatted incorrectly.")
        return

    # 2. Compute pairwise distances and write to TSV
    with open(args.output, 'w') as outfile:

        # Using itertools.combinations to generate all unique pairs
        for entry1, entry2 in itertools.combinations(entries, 2):
            dist = entry1.hamming_distance(entry2)

            # Output in Tab-Separated Format
            outfile.write(f"{entry1.header}\t{entry2.header}\t{dist:.4f}\n")

    # Calculate total distances safely below 79 chars
    total_dist = len(entries) * (len(entries) - 1) // 2

    print(
        f"Success! Calculated {total_dist} distances. "
        f"Saved to {args.output}"
    )


if __name__ == "__main__":
    # Declare a new parser
    parser = argparse.ArgumentParser(
        description=(
            "Read a Multi-FASTA file and output pairwise distances "
            "to a TSV file."
        )
    )

    # Add required arguments and their types
    parser.add_argument("input", type=str, help="Input FASTA file")
    parser.add_argument("output", type=str, help="Output TSV file")

    # The store_true action sets the variable to True if passed
    parser.add_argument(
        "--compressed",
        action='store_true',
        help="Force treat file as compressed"
    )

    # Read the arguments passed by the user and call the main function
    parsed_args = parser.parse_args()
    main(parsed_args)
