# todo: rodar esse codigo https://github.com/cbe90/self_equivalent_apn/blob/v1.1/8bit/8bit_class1/README.txt
# todo: investigar essa bomba https://zenodo.org/records/6983500

#= This code extracts the APN (Almost Perfect Nonlinear) matrix from a given S-box.
The function `extract_apn_matrix` takes an S-box represented as a vector of integers and
the number of bits `n` as input. It computes the APN matrix by iterating through pairs
of basis elements and applying the formula to calculate the entries of the matrix.
The resulting compact matrix is returned as a vector of integers.
The code also includes an example S-box and prints the extracted APN matrix. =#

# Reference: Weng et al. - On Quadratic Almost Perfect Nonlinear Functions and Their Related Algebraic Object

