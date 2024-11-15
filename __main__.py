import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a FASTA file to produce consensus sequence, cut alignment for best downstream phylogenetic analyses, and build phylogenetic tree. Read https://github.com/RiccPicc/consequence for more details.")

    # Mandatory arguments
    parser.add_argument('--input', type=str, required=True, help="Path to the input FASTA (sequences or alignment) or Newick (tree) file")
    parser.add_argument('--output', type=str, required=True, help="Path to the output file.")

    # Optional arguments with default values
    parser.add_argument('--find_best_position', type=bool, default=True, help="Whether to find the best position to cut alignment and consensus sequence (default: True).")
    parser.add_argument('--cutoff_best_base', type=float, default=0.85, help="Cutoff value for selecting the best base per position (default: 0.85).")
    parser.add_argument('--gap_count', type=int, default=10, help="Maximum allowed gap open (default: 10).")
    parser.add_argument('--gap_multiplier', type=float, default=None, help="Multiplier for gap open, overrides gap_count (default: None).")
    parser.add_argument('--nomenclature_cutoff', type=float, default=0.8, help="Cutoff value for selecting consensus base(s) per position (default: 0.8).")
    parser.add_argument('--no_cut_consensus', action='store_true', help="If set, do not cut consensus (default: False).")
    parser.add_argument('--remove_cutoff', type=float, default=None, help="Cutoff value for removal (default: None).")
    parser.add_argument('--cut_msa', type=bool, default=True, help="Whether to cut MSA (default: True).")
    parser.add_argument('--build_consensus', type=bool, default=True, help="Whether to build consensus (default: True).")
    parser.add_argument('--rename_leafs', type=bool, default=True, help="Whether to rename leaf nodes (default: True).")

    return parser.parse_args()

def main():
    args = parse_arguments()

    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Find best position: {args.find_best_position}")
    print(f"Cutoff best base: {args.cutoff_best_base}")
    print(f"Gap count: {args.gap_count}")
    print(f"Gap multiplier: {args.gap_multiplier}")
    print(f"Nomenclature cutoff: {args.nomenclature_cutoff}")
    print(f"No cut consensus: {args.no_cut_consensus}")
    print(f"Remove cutoff: {args.remove_cutoff}")
    print(f"Cut MSA: {args.cut_msa}")
    print(f"Build consensus: {args.build_consensus}")
    print(f"Rename leafs: {args.rename_leafs}")

if __name__ == "__main__":
    main()
