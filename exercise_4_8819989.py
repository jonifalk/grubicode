import argparse

__author__ = "Your Name"

# Declare a new parser
parser = argparse.ArgumentParser(description="Read a FASTA file and output a text file with aligned sequences.")
parser.add_argument("input", type=str, help="Input FASTA file")
parser.add_argument("output", type=str, help="Output text file")
parser.add_argument("match", type=int, help="Match value for scoring function")
parser.add_argument("mismatch", type=int, help="Mismatch value for scoring function")
parser.add_argument("gaps", type=int, help="Gap value for scoring function")
parser.add_argument("--local", action='store_true', help="Perform local alignment instead of global alignment")

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


def prepare_matrix_needleman_wish(scoring_func: dict):

    # Store scores and traceback information
    score = []
    traceback = []
    
    # TODO add matrix code
    return score, traceback


def prepare_matrix_smith_waterman(scoring_func: dict):
        # Store scores and traceback information
    score = []
    traceback = []
    
    # TODO add matrix code
    return score, traceback


def needleman_wish(scoring_func: dict, sequences: list):
    # TODO: Implement Smith-Waterman
    pass

def smith_waterman(scoring_func: dict, sequences: list):
    # TODO: Implement Smith-Waterman
    pass


def main(args):
    pass


if __name__ == "__main__":
    # Read the arguments passed by the user and call the main function
    args = parser.parse_args()
    main(args)