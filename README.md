# In-Silico-PCR And Needleman-Wunsch

## Pre-requisites: <br>
Set up a new conda environment with python 3.11 , Blast, Seqtk using the following commands- <br>
```
conda create -n myenv
```
```
conda install python=3.8
```
```
conda install -c bioconda blast
```
```
conda install -c bioconda seqtk
```

About ispcr: <br>
isPCR is simply using a computer to simulate what would happen if you performed a real PCR reaction with a <br>
certain set of primers and template sequence. As isPCR tries to replicate the results you would get in a real  <br>
PCR, it follows the same rules that govern the results of a real PCR reaction.  <br>
isPCR is performed in three steps:  <br>
1. Identify locations where primers would anneal to the target sequence  <br>
2. Identify pairs of locations where two primers anneal close enough together and in the correct  <br>
orientation for amplification to occur  <br>
3. Extract the amplified sequence  <br>

## ispcr module (ispcr.py) <br>
This ispcr module consists of a function named ispcr. Provided with primers and an assembly file, it gives back predicted amplicons in a single step! <br>

This function takes as input the path to a primer file, the path to an assembly file, and a maximum amplicon. <br>

Check the implementation and usage of the module by following the given steps: <br>
(a) Activate the conda environment consisting of above-mentioned tools by running "conda activate myenv" <br>
(b) Run the q1.py script in your terminal.


## Needleman Wunch algorithm module (nw.py)
The Needleman-Wunsch algorithm is a dynamic programming algorithm used for global sequence alignment. It aligns two sequences by maximizing the total score based <br>
on match, mismatch, and gap penalties. This implementation is written in Python and uses NumPy for efficient matrix operations. <br>

Check the implementation and usage of the script by following the given steps: <br>
Run the script q2.py in your terminal. <br>

## amplicon.py
This script should have the following characteristics: <br>
1. It imports "magnumopus" package to perform isPCR and Needleman-Wunsch alignment. <br>
2. It accepts as command line input two assembly files, a primer file, a maximum amplicon size, and a match, mismatch, and gap score. <br>
3. It has a command line interface that includes named arguments, not simply positional arguments. <br>
4. It performs isPCR on two provided assembly files with the provided primer file (in the data folder). <br>
5. The amplicons from each assembly file would be aligned. There should only be one amplicon from each file. <br>
6. It will check which orientation the two sequences align best in and only return the best alignment (i.e., reverse complements one sequence to check) <br>
7. It prints the alignment followed by the alignment score to the terminal. <br>

Check the implementation and usage of the module by following the given steps: <br>
(a) Activate the conda environment consisting of the above-mentioned tools by running "conda activate myenv" <br>
(b)Run the amplicon_align.py in your terminal with following flags: <br>
![image](https://github.com/hinagaur/In-Silico-PCR-Needleman-Wunsch/assets/66309991/99f3a021-2260-43c2-94ab-1f231e81b4c7)









