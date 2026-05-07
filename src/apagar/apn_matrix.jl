# todo: run this code https://github.com/cbe90/self_equivalent_apn/blob/v1.1/8bit/8bit_class1/README.txt
# todo: investigate this dataset https://zenodo.org/records/6983500

#= This code extracts the APN (Almost Perfect Nonlinear) matrix from a given S-box.
The function `extract_apn_matrix` takes an S-box represented as a vector of integers and
the number of bits `n` as input. It computes the APN matrix by iterating through pairs
of basis elements and applying the formula to calculate the entries of the matrix.
The resulting compact matrix is returned as a vector of integers.
The code also includes an example S-box and prints the extracted APN matrix. =#

# Reference: Weng et al. - On Quadratic Almost Perfect Nonlinear Functions and Their Related Algebraic Object


# todo: the tuple generation code is available now; the next step is to
# implement Algorithm 1 from the paper. They first run preprocessing with Propositions 4 and 5
# to reduce the search space, so this can be done first to validate the tuple counts reported in the paper.
#
# For n=6: these filters remove 8 of the 17 tuples immediately.
# For n=7: these filters remove 13 of the 27 tuples.
# For n=8: these filters remove 15 of the 32 tuples.
#
#
# According to Table 1 from the paper, running the 17 tuples through the filters should give exactly this behavior:
#
#    Proposition 4 checks whether affine subspace dimensions fall into the forbidden sizes 2, 4, or n - 1.
#    It should automatically discard Classes 6, 9, 13, 16, and 17.
#    Proposition 5 checks whether there is a quadrinomial that is a multiple of the minimum polynomials
#    of matrices A and B. It should automatically discard Classes 4, 8, and 12.

# todo: try rendering the matrices in Gray code to visualize their differences more clearly


# todo: https://github.com/cbe90/Supplementary-code-to-Trims-and-extensions-of-quadratic-APN-functions/blob/main/search_8bit_extensions/prefixes7.h
#https://eprint.iacr.org/2020/1515