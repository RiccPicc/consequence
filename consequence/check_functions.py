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
        if len(lines) != 2 or lines[0][0]==">":
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
        

def explicit(lines):
    good = False
    for line in lines:
        if line.strip(":")[0].lower() in ["taxa", "gene", "taxid", "taxon", "taxids", "genes"]:
            good = True
        else:
            good = False
    
    return good

def implicit(lines):
    return len(lines) == 2
        
def check_goodness(lines):
    if implicit(lines) or explicit(lines):
        return True
    
    return False

def extract_features(text, keyword, word_list):
    if text.startswith(f"{keyword}:"):
        text = text.replace("gene:", "").strip("\n").strip(";").strip(",").strip()
        word_list.extend(text.split(","))


def extract_genes(lines):
    genes = []

    for line in lines:
        extract_features(line, "gene", genes)
        extract_features(line, "genes", genes)

    return genes
    

def extract_taxa(lines):
    taxa = []

    for line in lines:
        extract_features(line, "taxon", taxa)
        extract_features(line, "taxa", taxa)
        extract_features(line, "taxid", taxa)
        extract_features(line, "taxids", taxa)

    return taxa