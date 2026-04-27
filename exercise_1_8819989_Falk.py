import argparse
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
        # Create a translation table for standard DNA/RNA bases
        complement_map = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")

        # Reverse the sequence using slicing ([::-1])(reversion)
        # Translate the bases to their complements using the translation table
        self.sequence = self.sequence[::-1].translate(complement_map)

        return self.sequence

    def __str__(self) -> str:
        """
        The __str__ function determines the output of str(<FastaEntry> object).
        Wraps the sequence to 60 characters per line to adhere to typical
        FASTA formatting guidelines.
        """
        # Wrap the sequence into lines of up to 60 characters
        wrapped_sequence = textwrap.fill(self.sequence, width=60)
        return f"{self.header}\n{wrapped_sequence}"


def main(args):
    """
    Main execution logic:
    reads input, processes the FASTA entry, and writes output.
    """
    # Read the input FASTA file
    with open(args.input, 'r') as infile:
        lines = infile.read().splitlines()

    # Error Handle for empty File
    if not lines:
        print("Error: Input file is empty.")
        return

    # Extract the header (the first line, starts with >)
    header = lines[0].strip()

    # Concatenate the remaining lines to form the complete sequence.
    # By using strip(), we ignore the line breaks.
    sequence = "".join([line.strip() for line in lines[1:]])

    # Instantiate the class object
    fasta_entry = FastaEntry(header=header, sequence=sequence)

    # Convert the sequence to its reverse complement
    fasta_entry.reverse_complement()

    # Output the modified FASTA object to the specified path
    with open(args.output, 'w') as outfile:
        # Writing str(fasta_entry) automatically triggers the __str__ method
        outfile.write(f"{fasta_entry}\n")


if __name__ == "__main__":
    # Declare a new parser
    parser = argparse.ArgumentParser(
        description="Read a FASTA file and output "
                    "a FASTA file with the reverse complement"
    )

    # Add required arguments and their types
    parser.add_argument("input", type=str, help="Input FASTA file")
    parser.add_argument("output", type=str, help="Output FASTA file")

    # Read the arguments passed by the user and call the main function
    args = parser.parse_args()
    main(args)
