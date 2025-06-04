from Bio import Phylo
from Bio.SeqIO import SeqIO
from io import StringIO

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_tree(filename):
    try:
        with open(filename, 'r') as f:
            content = f.read()
            tree = Phylo.read(StringIO(content), "newick")
        return True
    except Exception as e:
        return False

def is_text(filename):
    with open(filename, "r") as handle:
        lines = handle.readlines()
        if len(lines) <= 1 or len(lines) > 3 or lines[0][0]==">":
            return False

def assign_input(filename):
    match filename:
        case is_fasta(filename):
            return "fasta"
        case is_tree(filename):
            return "tree"
        case is_text(filename):
            return "text"
        case _:
            raise ValueError("Your input file is not a fasta, a newick tree, or a valid text file.")