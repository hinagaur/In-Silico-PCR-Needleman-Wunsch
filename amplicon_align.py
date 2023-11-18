#!/usr/bin/env python3
from magnumopus import *

import subprocess
import tempfile
import argparse

parser = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")
parser.add_argument("-1", dest="assembly1", required=True, help="Path to the first assembly file")
parser.add_argument("-2", dest="assembly2", required=True, help="Path to the second assembly file")
parser.add_argument("-p", dest="primers", required=True, help="Path to the primer file")
parser.add_argument("-m", dest="max_amplicon_size", type=int, required=True, help="Maximum amplicon size for isPCR")
parser.add_argument("--match", type=int, required=True, help="Match score to use in alignment")
parser.add_argument("--mismatch", type=int, required=True, help="Mismatch penalty to use in alignment")
parser.add_argument("--gap", type=int, required=True, help="Gap penalty to use in alignment")
    
args = parser.parse_args()

ASSEMBLY1 = args.assembly1
ASSEMBLY2 = args.assembly2
PRIMERS = args.primers
MAX_AMPLICON_SIZE = args.max_amplicon_size
MATCH = args.match
MISMATCH = args.mismatch
GAP = args.gap


ispcr_1 = ispcr(PRIMERS, ASSEMBLY1, MAX_AMPLICON_SIZE)
seq_a_ = ispcr_1.strip().split("\n")
seq_a = seq_a_[1]

ispcr_2 = ispcr(PRIMERS, ASSEMBLY2, MAX_AMPLICON_SIZE)
seq_b_ = ispcr_2.strip().split("\n")
seq_b = seq_b_[1]

# Reverse complementing seq_a
seq_a_rc = seq_a[::-1].translate(str.maketrans("ATGC", "TACG"))

# Reverse complementing seq_b
seq_b_rc = seq_b[::-1].translate(str.maketrans("ATGC", "TACG"))

# Alignment without reverse complementing
alignment_1, score_1 = needleman_wunsch(seq_a, seq_b, MATCH, MISMATCH, GAP)

# Alignment when reverse complementing seq_a
alignment_2, score_2 = needleman_wunsch(seq_a_rc, seq_b, MATCH, MISMATCH, GAP)

# Alignment when reverse complementing seq_b
alignment_3, score_3 = needleman_wunsch(seq_a, seq_b_rc, MATCH, MISMATCH, GAP)

# Alignment when reverse complementing seq_b and seq_a
alignment_4, score_4 = needleman_wunsch(seq_a_rc, seq_b_rc, MATCH, MISMATCH, GAP)

# Compare scores and print the alignment with the highest score
if score_1 > score_2 and score_1 > score_3 and score_1 > score_4:
    print("\n".join(alignment_1) + "\n")
    print(score_1)
elif score_2 > score_3 and score_2 > score_4:
    print("\n".join(alignment_2) + "\n")
    print(score_2)
elif score_3 > score_4:
    print("\n".join(alignment_3) + "\n")
    print(score_3)
else:
    print("\n".join(alignment_4) + "\n")
    print(score_4)



























# a = np.array([1,2,3]) # initializing 1d array

# # print(a)

# b = np.array([[9,0,8,0,7,0], [6,0,5,0,4,0]]) # initializing 2d array
# # print(b)

# print(b.shape)



# def init_matrix(num1: int, num2: int) -> list[list[list[int]]]:
#     string1= str(num1)
#     string2= str(num2)

#     matrix = []

#     for _ in range(len(string1)):
#         row = []
#         for _ in range(len(string2)):
#             row.append([])
            
#         matrix.append(row)

#     return matrix

# def main():
#     num1= 5433
#     num2= 789

#     init_matrix(num1, num2)

# if __name__ == '__main__': 
#     main()

# a= 12345
# b= 98765

# # making a matrix
# matrix = []
# for _ in range(len(str(a))):
#     row= []
#     for _  in range(len(str(b))):
#         row.append([])
#     matrix.append(row)

# # for x in matrix:
# #     print(x)

# for i in range(len(str(a))):
#     n_i = str(a)[i]
#     # print(n_i)
#     for j in range(len(str(b))):
#         n_j = str(b)[j]
#         # print(n_j)
#         prod = int(n_i)* int(n_j)
        
#         matrix[i][j] = [prod]
#         print(matrix[i][j])

# # for l in matrix:
# #     print(l)

