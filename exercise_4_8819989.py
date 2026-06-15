"""Exercise 4: Sequence Alignment"""
import argparse
import gzip
import os

__author__ = "Falk, 8819989"
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
        Generates the reverse complement of the sequence.

        Returns:
            str: The reverse complemented DNA sequence.
        """
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return "".join(
            comp.get(base, base) for base in reversed(self.sequence.upper())
        )

    def __str__(self) -> str:
        """Returns the basic string representation of the entry."""
        return f"{self.header}\n{self.sequence}"


def is_gzipped(filepath: str) -> bool:
    """
    Checks if a file is gzip-compressed by reading its magic bytes.

    Args:
        filepath: The path to the file.

    Returns:
        bool: True if the file is a gzip file, False otherwise.
    """
    if not os.path.exists(filepath):
        return False

    try:
        with open(filepath, 'rb') as f:
            magic_bytes = f.read(2)
            return magic_bytes == b'\x1f\x8b'
    except IOError:
        return False


def read_fasta(filepath: str) -> list[FastaEntry]:
    """
    Reads a Multi-FASTA file and parses it into a list of FastaEntry objects.

    Args:
        filepath (str): Path to the FASTA file.

    Returns:
        list: A list of parsed FastaEntry objects.
    """
    entries = []
    is_compressed = is_gzipped(filepath)
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
                if current_header:
                    sequence = "".join(current_sequence_lines)
                    entries.append(FastaEntry(current_header, sequence))

                current_header = line
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line)

        if current_header:
            sequence = "".join(current_sequence_lines)
            entries.append(FastaEntry(current_header, sequence))

    return entries


def prepare_matrix_needleman_wish(
        scoring_func: dict, seq1: str, seq2: str
) -> tuple[list[list[int]], list[list[str]]]:
    """
    Prepares the score matrix and traceback matrix for Needleman-Wunsch.

    Args:
        scoring_func (dict): Dictionary with match, mismatch, and gap scores.
        seq1: First sequence (rows in the matrix).
        seq2: Second sequence (columns in the matrix).

    Returns:
        tuple: (score_matrix, traceback_matrix)
    """
    match = scoring_func['match']
    mismatch = scoring_func['mismatch']
    gap = scoring_func['gap']

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score = [[0] * cols for _ in range(rows)]
    traceback = [['*'] * cols for _ in range(rows)]

    # Initialization of the matrix
    # The first row and column are initialized with gap penalties
    for i in range(1, rows):
        score[i][0] = i * gap
        traceback[i][0] = 'U'       # Up (Gap in sequence 2)
    for j in range(1, cols):
        score[0][j] = j * gap
        traceback[0][j] = 'L'       # Left (Gap in sequence 1)

    # Calculating Matrix entries
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                diag_score = score[i-1][j-1] + match
            else:
                diag_score = score[i-1][j-1] + mismatch

            up_score = score[i-1][j] + gap
            left_score = score[i][j-1] + gap

            # Determine the maximum score and record traceback direction.
            # Tie-breaking priority: Diagonal > Up > Left
            if diag_score >= up_score and diag_score >= left_score:
                score[i][j] = diag_score
                traceback[i][j] = 'D'
            elif up_score >= left_score:
                score[i][j] = up_score
                traceback[i][j] = 'U'
            else:
                score[i][j] = left_score
                traceback[i][j] = 'L'

    return score, traceback


def prepare_matrix_smith_waterman(
        scoring_func: dict, seq1: str, seq2: str
) -> tuple[list[list[int]], list[list[str]]]:
    """
    Prepares the score matrix and traceback matrix for Smith-Waterman.

    Args:
        scoring_func (dict): Dictionary with match, mismatch, and gap scores.
        seq1: First sequence (rows in the matrix).
        seq2: Second sequence (columns in the matrix).

    Returns:
        tuple: (score_matrix, traceback_matrix)
    """
    match = scoring_func['match']
    mismatch = scoring_func['mismatch']
    gap = scoring_func['gap']

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score = [[0] * cols for _ in range(rows)]
    traceback = [['*'] * cols for _ in range(rows)]

    # Initialization of the matrix
    # For Local Alignment, row 0 and col 0 start completely at 0
    # and traceback starts at stop symbol '*' (list comprehension)

    # Calculating Matrix entries
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                diag_score = score[i-1][j-1] + match
            else:
                diag_score = score[i-1][j-1] + mismatch

            up_score = score[i-1][j] + gap
            left_score = score[i][j-1] + gap

            # Local alignment allows paths to stop
            # completely by defaulting to 0
            max_score = max(0, diag_score, up_score, left_score)
            score[i][j] = max_score

            if max_score == 0:
                traceback[i][j] = '*'
            elif max_score == diag_score:
                traceback[i][j] = 'D'
            elif max_score == up_score:
                traceback[i][j] = 'U'
            else:
                traceback[i][j] = 'L'

    return score, traceback


def needleman_wish(scoring_func: dict, sequences: list) -> tuple[str, str]:
    """
    Performs global alignment using the Needleman-Wunsch algorithm.

    Args:
        scoring_func (dict): A dictionary with alignment scores.
        sequences (list): List containing FastaEntry objects.

    Returns:
        tuple: Two formatted sequences representing the alignment.
    """
    seq1 = sequences[0].sequence
    seq2 = sequences[1].sequence

    score, traceback = prepare_matrix_needleman_wish(
        scoring_func, seq1, seq2
    )

    aligned_seq1 = []
    aligned_seq2 = []

    # Backtracking starts in the bottom right corner
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if traceback[i][j] == 'D':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'U':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            break

    # The sequences were generated backwards, so we reverse them.
    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2))


def smith_waterman(scoring_func: dict, sequences: list) -> tuple[str, str]:
    """
    Performs local alignment using the Smith-Waterman algorithm.

    Args:
        scoring_func (dict): A dictionary with alignment scores.
        sequences (list): List containing FastaEntry objects.

    Returns:
        tuple: Two formatted sequences representing the alignment.
    """
    seq1 = sequences[0].sequence
    seq2 = sequences[1].sequence

    score, traceback = prepare_matrix_smith_waterman(
        scoring_func, seq1, seq2
    )

    max_val = -1
    max_i, max_j = 0, 0
    rows, cols = len(seq1) + 1, len(seq2) + 1

    # Locate the cell with the absolute highest score across the matrix
    for i in range(rows):
        for j in range(cols):
            if score[i][j] > max_val:
                max_val = score[i][j]
                max_i = i
                max_j = j

    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_i, max_j

    # Backtracking stops when a cell containing '0' (traceback '*') is hit.
    while traceback[i][j] != '*' and score[i][j] > 0:
        if traceback[i][j] == 'D':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'U':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            break

    # The sequences were generated backwards, so we reverse them.
    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2))


def main(args):
    """Main execution entrypoint for processing parameters and executing."""
    entries = read_fasta(args.input)

    if len(entries) < 2:
        print("Error: The FASTA file must contain at least two sequences.")
        return

    scoring_func = {
        'match': args.match,
        'mismatch': args.mismatch,
        'gap': args.gaps
    }

    # Decide alignment logic based on arguments
    if args.local:
        aligned_seq1, aligned_seq2 = smith_waterman(scoring_func, entries)
    else:
        aligned_seq1, aligned_seq2 = needleman_wish(scoring_func, entries)

    # Write output sequences to output file
    with open(args.output, 'w') as out_f:
        out_f.write(f"{aligned_seq1}\n{aligned_seq2}\n")

    print(f"Alignment completed. Saved format to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read a FASTA file and output aligned sequences."
    )

    # Required Positional Arguments
    parser.add_argument("input", type=str, help="Input FASTA file")
    parser.add_argument("output", type=str, help="Output text file")

    # The scores act positionally but fallback
    # using nargs='?' enabling default.
    parser.add_argument(
        "match", type=int, nargs="?", default=6,
        help="Match value for scoring function"
    )
    parser.add_argument(
        "mismatch", type=int, nargs="?", default=-2,
        help="Mismatch value for scoring function"
    )
    parser.add_argument(
        "gaps", type=int, nargs="?", default=-5,
        help="Gap value for scoring function"
    )

    # Optional Arguments
    parser.add_argument(
        "--local", action='store_true',
        help="Perform local alignment instead of global alignment"
    )

    parsed_args = parser.parse_args()
    main(parsed_args)
