import argparse
from .check_functions import *
from .consensus_functions import *
from .download_seqs import *
from .manage_names import *

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a FASTA file to produce consensus sequence, cut alignment for best downstream phylogenetic analyses, and build phylogenetic tree. Read https://github.com/RiccPicc/consequence for more details.")

    # Mandatory arguments
    parser.add_argument('--input', type=str, required=True, help="Path to the input FASTA (sequences or alignment) or Newick (.treefile output of iqtree2) formats.")
    parser.add_argument('--output_path', type=str, required=True, help="Path where to place output files.")

    # Optional arguments with default values
    parser.add_argument('--gene_name', type=str, default="", help="gene name you want to downalod.")
    parser.add_argument('--taxon', type=str, default="", help="Taxon you are looking for (species, genus, family, etc.).")
    parser.add_argument('--find_best_position', type=lambda x: x.lower() == 'true', default=True, help="Whether to find the best position to cut alignment and consensus sequence (default: True). Does not work with unaligned files if no_msa=True. ")
    parser.add_argument('--cutoff_best_base', type=float, default=0.85, help="Cutoff value for selecting the best base per position (default: 0.85).")
    parser.add_argument('--gap_count', type=int, default=10, help="Maximum allowed gap open (default: 10).")
    parser.add_argument('--gap_multiplier', type=float, default=None, help="Multiplier for gap open, overrides gap_count (default: None).")
    parser.add_argument('--nomenclature_cutoff', type=float, default=0.8, help="Cutoff value for selecting consensus base(s) per position (default: 0.8).")
    parser.add_argument('--no_cut_consensus', type=lambda x: x.lower() == 'true', default=False, help="If set, do not cut consensus (default: False).")
    parser.add_argument('--remove_cutoff', type=float, default=None, help="Cutoff value for removal (default: None).")
    parser.add_argument('--cut_msa', type=lambda x: x.lower() == 'true', default=True, help="Whether to cut MSA (default: True). Does not work with unaligned files if no_msa=True.")
    parser.add_argument('--build_consensus', type=lambda x: x.lower() == 'true', default=True, help="Whether to build consensus (default: True).")
    parser.add_argument('--rename_leafs', type=lambda x: x.lower() == 'true', default=True, help="Whether to rename leaf nodes (default: True).")
    parser.add_argument('--save_intermediate', type=lambda x: x.lower() == 'true', default=False, help="Whether to save all intermediate files (default: False).")
    parser.add_argument('--save_individually', type=lambda x: x.lower() == 'true', default=False, help="Whether to save all sequences of msa individually (default: False).")
    parser.add_argument('--save_species', type=lambda x: x.lower() == 'true', default=False, help="Whether to save consensus of species individually (default: False).")
    parser.add_argument('--no_msa', type=lambda x: x.lower() == 'true', default=False, help="Does not perform MSA in case of sequences input (default: False).")

    return parser.parse_args()

def main():
    args = parse_arguments()

    # read input
    infiles = args.input.split(",") # add function to read input and make different things accoring to outputs. Consider that input file can be more than one. 
    infiles_type = {"text":[], "fasta":[], "tree":[]}
    for i in infiles:
        infiles_type[assign_input(i)].append(i)

    out = args.output_path
    if not out.endswith("/"):
        out += "/"

    taxa = args.taxon.split(",")
    genes = args.gene_name.split(",")
    # what to do if text?
    # if text or taxa/gene in input proceed. If just one of the two present pass but raise warning
    if len(taxa) > 0:
        print(f"Found {",".join(taxa)} as input taxa")
    if genes != "":
        print(f"Found {",".join(genes)} as input genes") 

    if len(infiles_type["text"]) > 0:
        for i in infiles_type["text"]:
            with  open(i, "r") as handle:
                lines = [line for line in handle.readlines() if not line.isspace() or line[0] in [";",",",".",":","_",";"]]
                if check_goodness(lines):
                    taxa.extend(extract_taxa(lines))
                    genes.extend(extract_genes(lines))
                else:
                    raise Warning(f"{i} not a valid text format. Please make sure you have taxon names or taxon ids in the first line and genes in the second line. Both lines are mandatory.")
        
    # proceed to execute fasta (msa, consensus, trim gaps), if tree skip this 
    # proceed tree (iqtree, renaming, picture)

    match infiles_type:
        case "text":
            "TODO: check how is the text and prepare to download; list of all genes, taxa, taxid produced"
        case "fasta":
            "TODO: check if MSA, in such case skip alignment"
            fasta = infile
        case "tree":
            "TODO: skip to tree file management"

    run_download(gene_name=gene_name_input,
                 taxon_name=taxon_name_input,
                 taxid=taxid_input)
    
    # read file
    print("Reading sequences file")
    sequences = extract_sequences(fasta)
    
    print("Check nature of the file")

    # check input nature
    # if fasta -> make msa
    # if msa -> make cut/consensus
    is_msa = True # this has to be set with check_input function in future issue

    if args.no_msa:
        is_msa = False
        print("No alignment will be performed as --no_msa=True")
    # else:
    #   align with subprocess prank
    

    # get matrix
    print("Building positions matrix...")
    matrix = seqs_to_matrix(sequences)  

    # get start, end
    if args.find_best_position and is_msa:
        print("Finding best positions of start and end...")
        start, end = find_position(matrix, 
                                   cutoff=args.cutoff_best_base, 
                                   gap_count=args.gap_count, 
                                   gap_multiplier=args.gap_multiplier)
    else:
        start, end = 0, len(matrix)-1

    id_species_map = map_id_species(fasta)

    # build consensus
    if args.build_consensus:
        print("Building consensus of all sequences...")
        consensus = build_consensus(matrix, 
                                    start=start, 
                                    end=end, 
                                    nom_cutoff=args.nomenclature_cutoff,
                                    no_cut_consensus=args.no_cut_consensus)
        # save consensus as fasta
        save_fasta(consensus, out+"consensus.fasta", separated=False)
        print("Consensus saved")
        if args.save_species:
            print("Building consensus for each species...")
            species_consensus = make_species_consensus(sequences,id_species_map, 
                                   cutoff=args.cutoff_best_base, 
                                   gap_count=args.gap_count, 
                                   gap_multiplier=args.gap_multiplier, 
                                    nom_cutoff=args.nomenclature_cutoff,
                                    no_cut_consensus=args.no_cut_consensus)
            save_fasta(species_consensus, out+"species-wise_consensus.fasta", separated=False)
            print("All species consensus saved")
            
    # cut msa
    if args.cut_msa and is_msa:
        print("Removing long, uninformative gaps from MSA...")
        msa_cut, rem_cut = cut_msa(sequences, start=start, end=end)
        # save msa as fasta
        save_fasta(msa_cut, out+"cut_msa.fasta", separated=False)
        print("Final MSA saved")
        if args.remove_cutoff != None:
            save_fasta(rem_cut, out+"removed.fasta", separated=False)
            print("Clean sequences used for MSA saved individually")
    
    # use iqtree2 on cut_msa, species_consensus, msa_not_cut
    # rename_leafs on .treefile
    if args.rename_leafs:
        print("Map produced")
        # rename tree leafs
        tree_path = out+"cut_msa.fasta.treefile"
        print("Renaming leafs...")
        renamed_tree = rename_leafs(tree_path, id_species_map)
        print("Leafs renamed")
        renamed_tree_output = out+"renamed_tree.treefile"
        print("Returning output %s..." % renamed_tree_output)
        with open(renamed_tree_output, "w") as handle:
            handle.write(renamed_tree)
        print("%s created with success!" % renamed_tree_output)
    
    print("Run finished with success!")



if __name__ == "__main__":
    main()
