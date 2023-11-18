#!/usr/bin/env python3

import numpy as np

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int):
    """
    Perform global sequence alignment using the Needleman-Wunsch algorithm.

    Args:
    - seq_a (str): The first input sequence.
    - seq_b (str): The second input sequence.
    - match (int): Score assigned for a match between corresponding characters in the sequences.
    - mismatch (int): Score assigned for a mismatch between corresponding characters in the sequences.
    - gap (int): Score penalty for introducing a gap in one of the sequences.

    Returns:
    - tuple: A tuple containing two elements:
        - A tuple of two strings representing the aligned sequences.
        - An integer representing the alignment score.

    The function aligns the input sequences `seq_a` and `seq_b` based on the specified match, mismatch, and gap scores, and returns the aligned sequences along with the alignment score.
    """
    final_matrix = np.zeros((len(seq_a) + 1, len(seq_b) + 1))
    match_mismatch_matrix = np.zeros((len(seq_a), len(seq_b)))

    for i in range(len(seq_a)):
        for j in range(len(seq_b)):
            if seq_a[i] == seq_b[j]:
                match_mismatch_matrix[i][j] = match
            else:
                match_mismatch_matrix[i][j] = mismatch

    for i in range(len(seq_a) + 1):
        final_matrix[i][0] = i * gap
    for j in range(len(seq_b) + 1):
        final_matrix[0][j] = j * gap

    for i in range(1, len(seq_a) + 1):
        for j in range(1, len(seq_b) + 1):
            final_matrix[i][j] = max(final_matrix[i-1][j-1] + match_mismatch_matrix[i-1][j-1],
                                    final_matrix[i-1][j] + gap,
                                    final_matrix[i][j-1] + gap)

    aligned_seq_a = ""
    aligned_seq_b = ""

    ti = len(seq_a)
    tj = len(seq_b)

    while ti > 0 or tj > 0:
        current_score = final_matrix[ti][tj]
        
        if ti > 0 and tj > 0 and current_score == final_matrix[ti-1][tj-1] + match_mismatch_matrix[ti-1][tj-1]:
            aligned_seq_a = seq_a[ti-1] + aligned_seq_a
            aligned_seq_b = seq_b[tj-1] + aligned_seq_b
            ti -= 1
            tj -= 1
        elif ti > 0 and current_score == final_matrix[ti-1][tj] + gap:
            aligned_seq_a = seq_a[ti-1] + aligned_seq_a
            aligned_seq_b = "-" + aligned_seq_b
            ti -= 1
        else:
            aligned_seq_a = "-" + aligned_seq_a
            aligned_seq_b = seq_b[tj-1] + aligned_seq_b
            tj -= 1

    alignment_score = int(final_matrix[-1][-1])
    return (aligned_seq_a, aligned_seq_b), alignment_score


# Example usage:
if __name__ == "__main__":
    seq_a = "CATCGAGGAAGGCATCCGCGAAGTGATGAGCGCCATCGCCCAGTTCCCGGGCACGGTGGACAGCATCCTGGCCGACTACAATCGCATCGTCGCCGAAGGCGGTCGCCTCTCCGACGTCCTCAGCGGCTATATCGATCCCGATGACGGCAGCCTGCCCGCCGAAGAGGTGGAGCCGGTCAACCTGAAGGACGATTCCGCCGACTCGAAAGAGAAGGACGACGAGGAAGAAGAAAGCGACGACAGCAGCGACAGCGACGACGAAGGCGACGGCGGTCCGGATCCGGAAGAAGCCCGCCTGCGTTTCACCGCGGTCTCCGAGCAGCTCGACAAGGCCAAGAAGGCCCTGAAGAAGCACGGTCGCGGCAGCAAGCAGGCCACCGCCGAACTCACCGGCCTGGCCGAGCTGTTCATGCCGATCAAGCTGGTGCCCAAGCAGTTCGACGCCCTGGTCGCCCGCGTGCGCTCCGCCCTGGAAGGCGTGCGCGCCCAGGAACGCGCCATCATGCAGCTCTGCGTGCGTGACGCGCGCATGCCGCGTGCCGACTTCCTGCGCCTGTTCCCGAACCACGAGACCGACGAGAAGTGGGTCGACAGCGTCCTGAAGAGCAAGCCGAAGTACGCCGAGGCCATCGAGCGCCTGCGCGACGACATCCTGCGCAACCAGCAGAAGCTGGCGGCCCTGGAAAGCGAGGTCGAGCTGACCGTCGCCGAGA"
    seq_b = "TATCGAAGAGGGTATCCGTGAAGTGATGGGCGCAATCGCGCACTTCCCTGGCACGGTTGATCACATCCTCTCCGAGTACACTCGCGTCACCACCGAAGGTGGACGCCTGTCCGACGTTCTCAGCGGGTACATCGACCCGGATGACGGCATTACGCCGCCTGCCGCCGAAGTGCCGCCGCCTGTCGACACCAAGACCGCGAAAGCGGATGACGACAGCGAAGACGAAGAAGCGGAAGCGACCGAAGACGAAGAAGAAGCCGAAAGCGGTCCGGATCCGGTCATCGCCGCACAGCGCTTTGGCGCCGTCGCCGACCAGATGGAAGTCACCCGCAAGGCGCTGAAGAAACACGGCCGCGAGAACAAGCAAGCCATCGCCGAAATGCTGGCCCTGGCTGAACTGTTCATGCCGATCAAACTGGTTCCGAAGCAATTCGAAGGCCTGGTTGAACGTGTACGCAGTGCCCTGGATCGCCTGCGCCAGCAGGAACGCGCCATCATGCAGCTCTGCGTGCGTGATGCCCGCATGCCACGTGCCGACTTCCTGCGCCAGTTCCCGGGCAACGAAGTGGATGAAAGCTGGACCGACGCCCTGGCCAAAGGCAAGAGCAAGTACGCCGAAGCCATCGCCCGCCTGCAGCCGGACATCATCCGTTGCCAGCAGAAGCTGACCGCTCTTGAAGTCGAAACAGGCTTGAAGATCGCCGAGA"
    match = 1
    mismatch = -1
    gap = -1
    aligned_sequences, score = needleman_wunsch(seq_a, seq_b, match, mismatch, gap)
    print(aligned_sequences[0])
    print(aligned_sequences[1])
    print(score)


