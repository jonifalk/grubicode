"""Exercise 5: Tree Construction."""
import argparse
import math
import sys
import os
import gzip

__author__ = "Falk, 8819989"
__version__ = "1.0.0"

# Declare the parser at module level as given by the template
parser = argparse.ArgumentParser(
    description="Read a FASTA file detailing alignments and "
                "construct a tree using UPGMA."
)
parser.add_argument("input", type=str, help="Input FASTA Alignment")
parser.add_argument(
    "output", type=str,
    help="Output text file, contains the computed tree in Newick Format."
)


class FastaEntry:
    """A class to store and manipulate a single FASTA entry."""

    def __init__(self, header: str = "", sequence: str = ""):
        """
        Initialize the FastaEntry object.

        Args:
            header: The FASTA header line (starts with ">").
            sequence: The continuous nucleotide sequence.
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
        """
        Returns the basic string representation of the entry.

        Returns:
            str: Formatted string of header and sequence.
        """
        return f">{self.header}\n{self.sequence}"


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
            return f.read(2) == b'\x1f\x8b'
    except IOError:
        return False


def read_fasta(filepath: str) -> list[FastaEntry]:
    """
    Reads a Multi-FASTA file and parses it into FastaEntry objects.

    Args:
        filepath: Path to the FASTA file.

    Returns:
        list: A list of parsed FastaEntry objects.
    """
    entries = []
    is_comp = is_gzipped(filepath)
    open_func = gzip.open if is_comp else open
    mode = 'rt' if is_comp else 'r'

    with open_func(filepath, mode) as f:
        curr_head = ""
        curr_seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if curr_head:
                    seq = "".join(curr_seq)
                    entries.append(FastaEntry(curr_head, seq))

                # Save header without the leading '>'
                curr_head = line[1:]
                curr_seq = []
            else:
                curr_seq.append(line)

        if curr_head:
            seq = "".join(curr_seq)
            entries.append(FastaEntry(curr_head, seq))

    return entries


def compute_distance_matrix(
        entries: list[FastaEntry]
) -> list[list[float]]:
    """
    Computes a distance matrix based on Hamming distance and applies
    the J-Cantor correction model.

    Args:
        entries: List of FastaEntry objects.

    Returns:
        list[list[float]]: Symmetric distance matrix.
    """
    n = len(entries)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = entries[i].sequence
            seq2 = entries[j].sequence

            if len(seq1) != len(seq2):
                print("Error: Sequences must be aligned (equal length).")
                sys.exit(1)

            length = len(seq1)

            # Calculate absolute Hamming Distance (amount of mismatches)
            mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)

            # Proportion of differing sites
            p = mismatches / length

            # J-Cantor Correction
            if p >= 0.75:
                # Math limit: Sequence is too saturated to correct
                dist = float('inf')
            else:
                dist = -0.75 * math.log(1.0 - (4.0 / 3.0) * p)

            # Matrix is symmetric
            matrix[i][j] = dist
            matrix[j][i] = dist

    return matrix


def print_matrix(matrix: list[list[float]], labels: list[str]) -> None:
    """
    Prints the symmetric distance matrix to the console.

    Args:
        matrix: 2D float list representing distances.
        labels: Names of the taxa/sequences.
    """
    print("Distance Matrix (Hamming + J-Cantor):")
    max_label_len = max((len(lab) for lab in labels), default=0)

    # Print upper header for readability
    header = f"{'':<{max_label_len}} | "
    header += " | ".join(f"{lab[:10]:>10}" for lab in labels)
    print(header)
    print("-" * len(header))

    # Print row variables
    for i, row in enumerate(matrix):
        row_str = f"{labels[i]:<{max_label_len}} | "
        row_str += " | ".join(
            f"{val:>10.4f}" if val != float('inf') else f"{'inf':>10}"
            for val in row
        )
        print(row_str)
    print()


def run_upgma(matrix: list[list[float]], labels: list[str]) -> str:
    """
    Executes the UPGMA clustering algorithm on a given distance matrix.

    Args:
        matrix: 2D float distance matrix.
        labels: Names of the taxa.

    Returns:
        str: The constructed phylogenetic tree in Newick format.
    """
    clusters = list(labels)
    # Track the size of each cluster for distance weighting
    sizes = [1] * len(labels)

    # Deep copy !!! matrix so we don't mutate the original data
    curr_matrix = [row[:] for row in matrix]

    while len(clusters) > 1:
        n = len(clusters)
        min_dist = float('inf')
        min_i, min_j = -1, -1

        # Locate the pair of clusters with the minimum distance
        for i in range(n):
            for j in range(i + 1, n):
                if curr_matrix[i][j] < min_dist:
                    min_dist = curr_matrix[i][j]
                    min_i, min_j = i, j

        # Merge the two clusters into a new node format
        new_cluster = f"({clusters[min_i]},{clusters[min_j]})"
        new_size = sizes[min_i] + sizes[min_j]

        # Calculate distance from new cluster to all other clusters
        new_distances = []
        for k in range(n):
            if k != min_i and k != min_j:
                weight_i = sizes[min_i] * curr_matrix[min_i][k]
                weight_j = sizes[min_j] * curr_matrix[min_j][k]

                # UPGMA distances are scaled proportionally to cluster size
                avg_dist = (weight_i + weight_j) / new_size
                new_distances.append(avg_dist)

        # Filter out old (merged) clusters from matrix and list
        next_matrix = []
        kept = [k for k in range(n) if k not in (min_i, min_j)]

        for r in kept:
            row = [curr_matrix[r][c] for c in kept]
            next_matrix.append(row)

        # Append distances mapping from old clusters to the new cluster
        for r_idx, dist in enumerate(new_distances):
            next_matrix[r_idx].append(dist)

        # Add the new cluster row (distances to old + 0.0 to itself)
        next_matrix.append(new_distances + [0.0])

        # Overwrite references for the next iteration
        clusters = [clusters[k] for k in kept] + [new_cluster]
        sizes = [sizes[k] for k in kept] + [new_size]
        curr_matrix = next_matrix

    # Return the remaining single cluster node completed with a semicolon
    return f"{clusters[0]};"


def main(args):
    """
    Main function to read parameters and trigger matrix and tree workflows.
    """
    entries = read_fasta(args.input)
    if len(entries) < 2:
        print("Error: Require at least two sequences to construct a tree.")
        sys.exit(1)

    labels = [entry.header for entry in entries]

    # Output Problem 1: Build Distance Matrix and output to console
    dist_matrix = compute_distance_matrix(entries)
    print_matrix(dist_matrix, labels)

    # Output Problem 2: Perform UPGMA clustering to build tree topology
    newick_tree = run_upgma(dist_matrix, labels)
    print(f"Generated Newick Tree:\n{newick_tree}\n")

    # Output Problem 3: Store tree in correct format
    with open(args.output, 'w') as out_f:
        out_f.write(newick_tree + "\n")

    print(f"Success. Tree saved to {args.output}")


if __name__ == "__main__":
    parsed_args = parser.parse_args()
    main(parsed_args)
