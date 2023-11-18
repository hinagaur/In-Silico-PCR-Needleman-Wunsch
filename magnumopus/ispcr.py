#!/usr/bin/env python3
import subprocess
import tempfile

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> list[tuple[list[str]]]:
    my_blast_hits= step_one(primer_file, assembly_file)
    primer_pair = step_two(my_blast_hits, max_amplicon_size)
    ampliseq = step_three(primer_pair, assembly_file)   
    return ampliseq

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    """
    This function performs primer annealing simulation by searching for sequence matches between
    the provided primer file and assembly file using BLASTN. It filters and retains full-length
    hits with a percent identity of at least 80%.

    Args:
        primer_file (str): Path to the primer file.
        assembly_file (str): Path to the assembly file.

    Returns:
        list[list[str]]: A list of lists containing information about the retained hits.
    """
    blast_run = subprocess.run(f"blastn -query {primer_file} -subject {assembly_file} -task blastn-short -outfmt '6 std qlen' | awk '{{ if($3 >= 80 && $4 == $13) {{ print }} }}'", stdout= subprocess.PIPE, text=True, shell=True)
    x = blast_run.stdout.strip().split("\n")

    blast_hits = []
    for i in x:
        k = i.split("\t")
        blast_hits.append(k)
    blast_hits.sort(key=lambda blast_hits: int(blast_hits[8]))

    return blast_hits

def step_two(sorted_hits: list[list[str]], max_amplicon_size: int) -> list[tuple[list[str]]]:
    """
    Identifies primer pairs that meet the criteria for valid PCR amplicons and returns them as tuples.

    This function takes a list of sorted BLAST hits, where each hit is represented as a list of strings.
    It also accepts the maximum allowed amplicon size (max_amplicon_size) as an integer.

    For each pair of hits in the sorted_hits list, this function checks if the hits meet the following criteria:
    1. Both primers are pointing towards each other.
    2. The distance between the 3' end of the first primer and the 3' end of the second primer is within the max_amplicon_size.

    If both criteria are met, the pair of hits is considered a valid PCR amplicon, and it is added to the result list.

    Parameters:
    sorted_hits (list[list[str]]): A list of sorted BLAST hits, where each hit is represented as a list of strings.
    max_amplicon_size (int): The maximum allowed amplicon size in base pairs.

    Returns:
    list[tuple[list[str]]]: A list of tuples, each containing two hits that form a valid PCR amplicon.
    """
    my_list = []

    for i in range(len(sorted_hits)):
        if int(sorted_hits[i][8]) < int(sorted_hits[i][9]):
            # Assigning names to the columns
            sstart_1, ssend_1 = int(sorted_hits[i][8]), int(sorted_hits[i][9])
            
            for k in range(i+1, len(sorted_hits)):
                if int(sorted_hits[k][8]) > int(sorted_hits[k][9]):
                    # Assigning names to the columns
                    sstart_2, ssend_2 = int(sorted_hits[k][9]), int(sorted_hits[k][8])
                    amp_len = sstart_2 - ssend_1
                    
                    if 0 < amp_len < max_amplicon_size:
                        my_list.append((sorted_hits[i], sorted_hits[k]))
    return my_list

def step_three(hit_pairs: list[tuple[list[str]]], assembly_f: str) -> str:
    """
    Extracts amplicon sequences from the provided assembly file based on the hit pairs.

    Parameters:
    hit_pairs (list[tuple[list[str]]): A list of hit pairs, where each hit pair is represented as a list of
                                       elements containing information about the amplicon location.
    assembly_file (str): The filepath to the assembly file containing the genomic sequences.

    Returns:
    str: A string containing the extracted amplicon sequences in FASTA format.

    The function extracts sequences from the assembly file for each hit pair specified in the hit_pairs list
    and returns the sequences in FASTA format.

    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.bed', delete=False) as tf:
        for i in hit_pairs:
            id = i[0][1]
            start = i[0][9] #start = i[0][8]
            end = int(i[1][9]) - 1
            tf.write(f"{id}\t{start}\t{end}\n")

    seq_tk = subprocess.run(f'seqtk subseq {assembly_f} {tf.name}', stdout=subprocess.PIPE, text=True, shell=True)

    return seq_tk.stdout



