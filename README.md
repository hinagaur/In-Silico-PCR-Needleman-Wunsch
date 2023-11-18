# In-Silico-PCR And Needleman-Wunsch

## Pre-requisites: <br>
Set up a new conda environment with python 3.11 , Blast, Seqtk using the following commands- <br>
conda create -n myenv
conda install python=3.8 <br>
conda install -c bioconda blast <br>
conda install -c bioconda seqtk <br>

## ispcr module (ispcr.py) <br>
This ispcr module consists of a function named ispcr. Provided with primers and an assembly file, it gives back predicted amplicons in a single step! <br>

This function takes as input the path to a primer file, the path to an assembly file, and a maximum amplicon. <br>

Check the implementation and usage of the script by following the given steps: <br>
Activate the conda environment consisting of above mentioned tools by running "conda activate myenv" <br>
Run the q1.py script <br>

About ispcr: <br>
isPCR is simply using a computer to simulate what would happen if you performed a real PCR reaction with a <br>
certain set of primers and template sequence. As isPCR tries to replicate the results you would get in a real  <br>
PCR, it follows the same rules that govern the results of a real PCR reaction.  <br>
isPCR is performed in three steps:  <br>
1. Identify locations where primers would anneal to the target sequence  <br>
2. Identify pairs of locations where two primers anneal close enough together and in the correct  <br>
orientation for amplification to occur  <br>
3. Extract the amplified sequence  <br>









