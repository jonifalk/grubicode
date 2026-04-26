import argparse

__author__ = "8819989, JonathanFalk"

# Declare a new parser
parser = argparse.ArgumentParser(description="Read a FASTA file and output a FASTA file with the reverse complement")
# Add required arguments and their types. The help messages will be shown when you run python exercise_1_Matrikelnr.py -h
parser.add_argument("input", type=str, help="Input FASTA file")
parser.add_argument("output", type=str, help="Output FASTA file")

class FastaEntry():
    def __init__(self):
        # You should the define the following attributes:
        # sequence (str)
        # header (str)
        pass

    def reverse_complement():
        # TODO: return the revcomp of the sequence
        pass


    def __str__(self):
        """
        The __str__ function determines the output of str(<FastaEntry> object). You can use this for shiny debugging or even when writing to files.
        """
        pass

def main(args):
    print(args.input)
    pass


if __name__ == "__main__":
    # Read the arguments passed by the user and call the main function
    args = parser.parse_args()
    main(args)
