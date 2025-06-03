# consequence

⚠️ **UNDER DEVELOPMENT** ⚠️

`consequence` (pronunced /ˈkɒnˈsēkwəns/) is a bioinformatic tool that wraps other tools and scirpts to produce alignments, consensus sequences and phylogenetic trees.

## Dependencies
- Python 3.10+
- Biopython
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [IQ-TREE2](http://www.iqtree.org/#download)

For windows user is highly reccomanded to use [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and download Linux (Intel) versions of  [Mafft](https://mafft.cbrc.jp/alignment/software/) and [IQ-TREE2](http://www.iqtree.org/#download)

## Features
- [ ] Useable from command line
- [ ] Save all intermediate files
- [ ] Download data using NCBI queries
- [ ] Recognition of fasta file and perform alignment if needed
- [x] Build map sequence ID - species
- [x] Produces a consensus sequence
- [x] Finds best position to cut the MSA and the consensus for better performance of downstream phylogenetic analyses
- [x] Removes from MSA uninformative sequences
- [x] Customization parameters to better fit your interests
- [ ] Possibility to remove internal long gaps
- [ ] Concatenate sequences based on metadata
- [ ] Build consensus species by species
- [ ] Concatenate species consesnsus sequences
- [ ] Build phylogenetic tree (iqtree2)
- [x] Rename phylogenetic tree leafs (from sequence ID to species name)
- [ ] Implement [barcoding gap analysis script](https://github.com/AleTatti/Barcoding-Analysis)
- [ ] Interactive GUI

## Example
TODO

## Customization parameters 
`consequence` lets you customize paramaters for a full control of the dataset you need for your analysis. The aim is to produce a consensus sequence and/or an MSA ideal for your downstream analyses.

**REMINDER**: After each of these parameters you have to insert the value withuot "=", just a space. For example `python -m consequence --input test/my_sequences.fasta --output_path ./test/results/ --rename_leafs False` 

<details><summary><b><code>--find_best_position</code></b></summary>

default: True

This parameters allows to find best position for start and stop where to perform cuts and obtain meaningful sequences for following phylogenetic analyses.
  
</details>

<details><summary><b><code>--cutoff_best_base</code></b></summary>

defult: 0.85

The objective is that start and end of the sequence must mostly be sequences of bases in most sequences of the MSA. This parameter allows to assess if a position is considerable a gap or a base throughout the MSA.

Toy example:

```
----ATGCATCAGC----CG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------GATCGAC-----
----------------TACG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------------------
-----------AGCAGT--G----ATG-----------------------------GCTAGCATCGAT------CAGTCAGAGAA------------------
------------------------ATG--------GTGATCGATCGACTAG-----------------------------GAGAA------GATCGAC-----
----------------TACG----ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT------CAGTCAGAGAA------------------
------------------------ATG--------GTGATCGATCGACTAG-----GCTAGCATCGAT--------------GAA-----------AC-----
----------------T-------ATG--------GTGATCGATCGACTAG-----GCTA---TCGAT------CAGTCAGAGAA--------T---------
```
Start → A at position 5

End → A at position 24

</details>
  
<details><summary><b><code>--gap_count</code></b></summary>
  
default: 10

This parameter allows to maintain the gap open. It is the number of gaps that the user tolerates after a base or a sequence of bases is retrieved.

Toy example:

sequence = `----ATGCATCAGC-------CG-----CG----ATGTGCTGCTCGCTACATCATCGATCGAT--------`

setting `gap_count = 10` → `ATGCATCAGC-------CG-----CG----ATGTGCTGCTCGCTACATCATCGATCGAT`

setting `gap_count = 6` → `CG-----CG----ATGTGCTGCTCGCTACATCATCGATCGAT`

setting `gap_count = 3` → `ATGTGCTGCTCGCTACATCATCGATCGAT`
</details>
  
<details><summary><b><code>--gap_multiplier</code></b></summary>
  
default: None

This parameter accepts an integer or a float number. It overwrites `gap_count` as the gap open tolerance is modulated by the length of the preceeding sequence found. 

Toy example:
sequence = `----ATGCATCAGC----TGCCAA----CG----` (sequence composed by 4 gaps, 10 bases, 4 gaps, 6 bases, 4 gaps, 2 bases, 4 gaps)

setting `gap_multiplier = 1` → `ATGCATCAGC----TGCCAA----CG`

setting `gap_multiplier = .5` → `ATGCATCAGC----TGCCAA`

setting `gap_multiplier = .2` → `ATGCATCAGC`
</details>
  
<details><summary><b><code>--nomenclature_cutoff</code></b></summary>
  
default: 0.8

This parameter is used to select most representative bases per position and eventually translate in IUPAC nomenclature.

Toy example:
```
ATGCAT-AGT
AAGCATCAGT
ATGCAT-A-C
ATGCAT-A-T
ATGCAT-AGC
AGGA-T-A-T
ATGCATCAGC
AAG-AT-AGT
AGGCA-CAGC
AAG-AT-AGT
```
setting `nom_cutoff = 1` consensus result → `ATGCAT-AGT`

setting `nom_cutoff = .8` consensus result → `A`**`W`**`GCAT-AG`**`Y`**

setting `nom_cutoff = .2` consensus result → `A`**`D`**`GCAT`**`C`**`AG`**`Y`**

</details>

<details><summary><b><code>--no_cut_consensus</code></b></summary>

Default: False
This parameters does not allow the cut of the consensus, so it will be produced based on the whole MSA obtained from sequences or given as input.
  
</details>

<details><summary><b><code>--cut_msa</code></b></summary>

Default: False
This parameters allows to cut MSA. If not set or set as True, consequence will cut the MSA based on best positions found and perform the following phylogenetic analysis. If set False, following analyses and output produced will base on the whole MSA obtained
  
</details>

<details><summary><b><code>--build_consensus</code></b></summary>

Default: True
This parameters allows the building of consensus sequence based on input msa or sequences.
  
</details>
  
<details><summary><b><code>separated</code></b></summary>
  
default: False

If set True allows to save cut alignment sequences individually.

</details>
  
<details><summary><b><code>--remove_cutoff</code></b></summary>
  
default: None

This parameters is used to remove sequences from MSA with high percentage of gaps in the selected region. 

Toy example:

```
TACG----ATGGTGAT-----------
--------ATGGTGATCGATC-ACTAG
T-------ATGGTG---GATCGACTAG
--CG----ATGGATCGAC---------
TACG----ATG----------------
T--G----ATG----------------
--------ATGGATCGAC---------
TACG----ATGGTGAT-----------
--------ATGGTGATCGATCGACTAG
T-------ATGGTGATCGATCGACTAG
```
setting `remove_cutoff = 0.8` removes line 6 (`T--G----ATG----------------`)

setting `remove_cutoff = 0.5` removes 6 lines (1, 4, 5, 6, 7, 8)

setting `remove_cutoff = 0.3` removes 8 lines (last two lines remain)
</details>

<details><summary><b><code>--rename_leafs</code></b></summary>
  
default: True

This parameters allow the renaming of phylogenetic tree leafs from sequences ID into species scientific names.

</details>

<details><summary><b><code>--save_individually</code></b></summary>
  
default: False

This parameters allow to save all MSA sequences into single files. If --no_msa False, input sequences will be saved individually (unaligned).

</details>

<details><summary><b><code>--save_intermediate</code></b></summary>
  
default: False

This parameters allow to save all intermediate files created during the run.

</details>

<details><summary><b><code>--no_msa</code></b></summary>
  
default: False

This parameters does not allow to perform MSA in case input are uniligned sequences.

</details>
