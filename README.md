# consequence

⚠️ **UNDER DEVELOPMENT** ⚠️

`consequence` is a `python3` tool that produces a consensus sequences based on a MSA fasta.

## Dependencies
- Python 3.10+
- Biopython

## Features
- [x] Produces a consensus sequence
- [x] Finds best position to cut the MSA and the consensus for better performance of downstream phylogenetic analyses
- [x] Removes from MSA uninformative sequences
- [x] Complete customisation of parameters to better fit your interests
- [ ] Useable from command line
- [ ] Recognition of fasta file and perform alignment if needed
- [ ] Possibility to remove internal long gaps
- [ ] Concatenate sequences based on metadata
- [ ] Build consensus species by species
- [ ] Concatenate species consesnsus sequences
- [ ] Build map sequence ID - species
- [ ] Build phylogenetic tree (iqtree2)
- [x] Rename phylogenetic tree leafs (from sequence ID to species name)
- [ ] Download data using NCBI queries
- [ ] Implement [barcoding gap analysis script](https://github.com/AleTatti/Barcoding-Analysis)

## Example
TODO

## Functions and customization parameters 
`consequence` let you customize paramaters for a full control of the dataset you need for your analysis. The aim is to produce a consensus sequence and/or an MSA ideal for your downstream analyses.
### `extract_sequences`
`extract_sequences` reads fasta files containing alignment data and store it as dictionary. 
### `seqs_to_matrix`
`seqs_to_matrix` translates the sequences into a matrix and transpose it. Needed to produce the consensus sequences and retrieve best position to cut off uninformative pieces of sequences.
### `find_position`
`find_position` aims to retrieve the most interesting positions that produce a valid, informative dataset for downstram phylogenetic analysis. It removes misleading, uninformative starting and ending pieces of sequences.
#### `cutoff`
defult: 0.85

The objective is that start and end of the sequence must mostly be sequences of bases in most sequences of the MSA. This parameter allows to assess if a position is considerable a gap or a base throughout the MSA.

Example (not real):

```
----ATGCATCAGC----CG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------GATCGAC-----
----------------TACG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------------------
-----------AGCAGT--G----ATG-----------------------------GCTAGCATCGAT------CAGTCAGAGAA------------------
------------------------ATG--------GTGATCGATCGACTAG-----------------------------GAGAA------GATCGAC-----
----------------TACG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------------------
------------------------ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT--------------GAA-----------AC-----
----------------T-------ATG--------GTGATCGATCGACTAG-----GCTA---TCGAT------CAGTCAGAGAA--------T---------
```
Start -> A at position 5

End -> A at position 24


#### `gap_count`
default: 10

This parameter allows to maintain the gap open. It is the number of gaps that the user tolerates after a base or a sequence of bases is retrieved.

