"""Standalone Task Exercise 3"""
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

    def __str__(self) -> str:
        """Returns the basic string representation of the entry."""
        return f"{self.header}\n{self.sequence}"


def is_gzipped(filepath: str) -> bool:
    """
    Checks if a file is gzip-compressed by reading its magic bytes.

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


def naive_search(text: str, pattern: str, max_mismatches: int) -> list[int]:
    """
    Performs a naive pattern search supporting approximate matching.

    Args:
        text (str): The reference text to search in.
        pattern (str): The pattern to find.
        max_mismatches (int): The maximum allowed mismatches.

    Returns:
        list[int]: A list of starting indices where the pattern occurs.
    """
    matches = []
    n, m = len(text), len(pattern)

    if n < m:
        return matches

    for i in range(n - m + 1):
        mismatches = 0
        for j in range(m):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            matches.append(i)

    return matches


def rabin_karp_search(
        text: str, pattern: str, max_mismatches: int
) -> list[int]:
    """
    Uses the Rabin-Karp algorithm for pattern matching.

    A rolling hash effectively accelerates exact matching.
    If mismatches are allowed, the hash check is bypassed and evaluates
    windows directly to accommodate approximate matches.

    Args:
        text (str): The reference text to search in.
        pattern (str): The pattern to find.
        max_mismatches (int): The maximum allowed mismatches.

    Returns:
        list[int]: A list of starting indices where the pattern occurs.
    """
    matches = []
    n, m = len(text), len(pattern)

    if n < m:
        return matches

    # Parameters for the rolling hash
    d = 256
    q = 1000000007  # A large prime number
    h = 1

    # Calculate h = pow(d, m - 1) % q
    for _ in range(m - 1):
        h = (h * d) % q

    p_hash = 0
    t_hash = 0

    # Calculate initial hash values for pattern and the first text window
    for i in range(m):
        p_hash = (d * p_hash + ord(pattern[i])) % q
        t_hash = (d * t_hash + ord(text[i])) % q

    for i in range(n - m + 1):
        # If max_mismatches is 0 use the hash as a fast filter.
        # Otherwise, evaluate to handle the requested mismatches.
        if max_mismatches > 0 or p_hash == t_hash:
            mismatches = 0
            for j in range(m):
                if text[i + j] != pattern[j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                matches.append(i)

        # Calculate hash for the next window of text
        if i < n - m:
            t_hash = (d * (t_hash - ord(text[i]) * h) + ord(text[i + m])) % q
            # We might get a negative value, convert it to positive
            if t_hash < 0:
                t_hash += q

    return matches


class SuffixTrieNode:
    """A node belonging to the SuffixTrie structure."""

    def __init__(self):
        # Dictionary mapping a character to a child SuffixTrieNode
        self.children = {}
        # Stores the starting positions of all suffixes passing through node
        self.positions = []


class SuffixTrie:
    """A naive suffix trie implementation used for pattern searching."""

    def __init__(self, text: str):
        self.root = SuffixTrieNode()
        self.text = text
        self._build_suffix_trie()

    def _build_suffix_trie(self):
        """Builds a naive suffix trie from the text."""
        n = len(self.text)

        # Iterate over all suffixes
        for i in range(n):
            node = self.root
            node.positions.append(i)

            # Insert the current suffix into the trie
            for j in range(i, n):
                char = self.text[j]
                if char not in node.children:
                    node.children[char] = SuffixTrieNode()

                node = node.children[char]
                # Record the original starting position of the suffix
                node.positions.append(i)

    def search_with_mismatches(
            self, pattern: str, max_mismatches: int
    ) -> list[int]:
        """
        Searches the trie for occurrences of the pattern using DFS.

        Args:
            pattern (str): The sequence pattern to find.
            max_mismatches (int): The number of allowed mismatches.

        Returns:
            list[int]: A sorted list of match start positions.
        """
        matches = set()

        def dfs(node: SuffixTrieNode, pattern_idx: int, mismatches_left: int):
            # If we've reached the end of the pattern, it's a valid match
            if pattern_idx == len(pattern):
                for pos in node.positions:
                    matches.add(pos)
                return

            target_char = pattern[pattern_idx]

            # Explore child nodes using Depth-First-Search
            for child_char, child_node in node.children.items():
                if child_char == target_char:
                    # Exact match for this character
                    dfs(child_node, pattern_idx + 1, mismatches_left)
                elif mismatches_left > 0:
                    # Approximate match: consume one allowed mismatch
                    dfs(child_node, pattern_idx + 1, mismatches_left - 1)

        # Start DFS from the root node
        dfs(self.root, 0, max_mismatches)

        return sorted(list(matches))


def main(args):
    """Main execution logic for the pattern matching algorithms."""
    # Parse the FASTA file
    entries = read_fasta(args.input)

    if not entries:
        print("Warning: Input file is empty or formatted incorrectly.")
        return

    # Prepare for output
    with open(args.output, 'w') as outfile:
        # Write column headers
        outfile.write("Pattern\tHeader\tPosition\n")

        total_matches = 0

        for entry in entries:
            sequence = entry.sequence

            # Select the requested matching algorithm
            if args.use_tree:
                trie = SuffixTrie(sequence)
                matches = trie.search_with_mismatches(
                    args.pattern,
                    args.allowed_mismatches
                )
            elif args.use_kmer:
                matches = rabin_karp_search(
                    sequence,
                    args.pattern,
                    args.allowed_mismatches
                )
            else:
                matches = naive_search(
                    sequence,
                    args.pattern,
                    args.allowed_mismatches
                )

            # Write the matches to the TSV file
            for match_pos in matches:
                outfile.write(f"{args.pattern}\t{entry.header}\t{match_pos}\n")
                total_matches += 1

    print(
        f"Success! Found {total_matches} occurrences. "
        f"Saved to {args.output}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read a Multi-FASTA file to find all occurrences of a "
                    "pattern."
    )

    # Positional Arguments
    parser.add_argument(
        "input", type=str, help="Input FASTA file")
    parser.add_argument(
        "pattern", type=str, help="Pattern to search for"
    )
    parser.add_argument(
        "output", type=str, help="Output TSV")

    # Optional Arguments
    parser.add_argument(
        "-allowed_mismatches", "--allowed_mismatches",
        type=int, default=0,
        help="Allowed mismatches for approximate pattern matching."
    )

    # Flags for specific algorithms
    parser.add_argument(
        "--use_kmer", action="store_true",
        help="Use Rabin-Karp rolling hash for matching."
    )
    parser.add_argument(
        "--use_tree", action="store_true",
        help="Build and use a naive Suffix Trie for matching."
    )

    # Read the arguments passed by the user and call the main function
    parsed_args = parser.parse_args()
    main(parsed_args)
