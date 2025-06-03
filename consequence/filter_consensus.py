from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import warnings

IUPAC_CODES = {
    frozenset(["a", "t"]): "w",
    frozenset(["a", "c"]): "m",
    frozenset(["a", "g"]): "r",
    frozenset(["t", "c"]): "y",
    frozenset(["t", "g"]): "k",
    frozenset(["c", "g"]): "s",
    frozenset(["a", "t", "g"]): "d",
    frozenset(["a", "t", "c"]): "h",
    frozenset(["g", "t", "c"]): "b",
    frozenset(["a", "c", "g"]): "v"
}

def extract_sequences(fasta):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    return sequences

def count_base(base, row):
    return row.count(base)

def get_position_dictionary(matrix):
    """
    Given a matrix, this function creates a dictionary containing a dictionary of all bases counted in that position.
    """
    positions = {}
    i = 0
    while i < len(matrix):
        bases = {}
        for base in list("atgc"):
            bases[base] = count_base(base, matrix[i])
        bases["gaps"] = count_base("-", matrix[i])
        bases["other"] = len(matrix[i]) - (bases["a"]+bases["c"]+bases["g"]+bases["t"]+bases["gaps"])
        bases["non-gaps"] = len(matrix[i]) - bases["gaps"]
        bases["length"] = len(matrix[i])
        positions[i] = bases
        i += 1
    return positions

def seqs_to_matrix(sequences):
    matrix = [list(str(seq.seq).lower()) for seq in sequences.values()]
    transposed_matrix = list(zip(*matrix))
    final_matrix = [list(row) for row in transposed_matrix]
    return final_matrix

def inspect_seqs(matrix):
    """
    This function inspect all sequences provided as input and returns a dictionary mapping the frequency of each bases present in each position.
    """
    pos_dict = get_position_dictionary(matrix)
    return pos_dict

def find_position(matrix, cutoff=.85, gap_count=10, gap_multiplier=None):
    pos_dict = inspect_seqs(matrix)

    target = 0
    counter = 0
    gap_counter = 0
    positions = list(pos_dict.keys())
    positions_details = {}

    gap_relative = False
    if isinstance(gap_multiplier, int) or isinstance(gap_multiplier, float):
        gap_relative = True
    
    for pos in positions:
        bases = pos_dict[pos]["non-gaps"] / pos_dict[pos]["length"]
        
        if bases > cutoff:
            if target == 0:
                target = pos
            counter += 1
            if gap_counter > 0:
                counter += gap_counter
                gap_counter = 0
                if gap_relative:
                    gap_count = counter * gap_multiplier
        elif counter > 0 and gap_counter < gap_count:
            gap_counter += 1
        else:
            positions_details[target] = counter
            target = 0
            counter = 0
            gap_counter = 0

    for pos in positions_details:
        if positions_details[pos] > gap_count:
            start = pos
            break

    for pos in list(positions_details.keys())[::-1]:
        if positions_details[pos] > gap_count:
            end = pos + positions_details[pos]
            break

    return start, end

def get_consensus(bases, nom_cutoff):
    if nom_cutoff == 0:
        raise ValueError("Misleading nomenclature cutoff inserted. You can't set it as 0. Please provide a greater number.")
    if nom_cutoff < .5:
        warnings.warn("Cutoff value provided is very low. Consensus result might be misleading.")

    count_bases = Counter(bases)
    sorted_bases = count_bases.most_common()
    highest_count = sorted_bases[0][1]
    
    top_bases = [nucleotide for nucleotide, count in sorted_bases if count >= highest_count*nom_cutoff]

    if len(top_bases) > 1 and "-" in top_bases:
        top_bases.remove("-")

    if len(top_bases) == 1:
        return top_bases[0]

    if len(top_bases) > 1:
        base_pair = frozenset(top_bases)
        return IUPAC_CODES.get(base_pair, "n")

    return "n"

def build_consensus(matrix, start=0, end=None, no_cut_consensus=False, nom_cutoff=.8):
    """
    nom_cutoff is the parameter used to establish if two bases have to be considered as one. For example if there are 8 "a" and 7 "t", consensus is "a" if nom_cutoff is 1, "w" if nom_cutoff is lower (like 0.8).
    """
    if not no_cut_consensus:
        if start == 0:
            warnings.warn("Start not provided for consensus cut, default is 0.")

        if end is None:
            warnings.warn("End not provided for consensus cut, default is the whole length of the MSA provided.")
            end = len(matrix)

    consensus = ""

    for row in matrix:
        consensus += get_consensus(row, nom_cutoff)

    return SeqRecord(Seq(consensus[start:end]), id="consensus", description="Consensus sequences from the original MSA")

def cut_msa(sequences, start=0, end=None, remove_cutoff=None):
    remove_low_informative = False

    if start == 0:
        if end is None:
            warnings.warn("You did not provide value to cut the MSA. Returning original MSA.")
            return sequences
        else:
            warnings.warn("Start not provided for MSA cut, default is 0.")

    if end is None:
        warnings.warn("End not provided for MSA cut, default is the whole length of the MSA provided.")
        end = len(next(iter(sequences.values())).seq)

    if remove_cutoff != None:
        if not isinstance(remove_cutoff, float) or remove_cutoff >= 1 or remove_cutoff <= 0:
            raise ValueError("Wrong remove_cutoff value %s. Insert a float number between 0 and 1 (for example 0.5)" % remove_cutoff)
        else:
            remove_low_informative = True

    cut_msa = {}
    removed = 0
    rem_seqs = {}
    for seq_id, record in sequences.items():
        cut_seq = str(record.seq)[start:end]
        gaps = cut_seq.count("-")
        if remove_low_informative:
            if gaps > len(cut_seq)*remove_cutoff:
                removed += 1
                cut_record = SeqRecord(Seq(cut_seq), id=record.id, description=record.description)
                rem_seqs[seq_id] = cut_record
                continue       
        cut_record = SeqRecord(Seq(cut_seq), id=record.id, description=record.description)
        cut_msa[seq_id] = cut_record

    if remove_low_informative:
        print("%s sequences removed because of low informative (too many gaps)" % removed)
    return cut_msa, rem_seqs

def save_fasta(sequences, output, separated=False):
    if isinstance(sequences, SeqRecord):
        with open(output, "w") as handle:
            SeqIO.write(sequences, handle, "fasta")
    elif isinstance(sequences, dict):
        with open(output, "w") as handle:
            SeqIO.write(sequences.values(), handle, "fasta")
        if separated:
            for record in sequences.values():
                with open(record.id+".fasta", "w") as handle:
                    SeqIO.write(record, handle, "fasta")