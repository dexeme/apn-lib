# todo: rodar esse codigo https://github.com/cbe90/self_equivalent_apn/blob/v1.1/8bit/8bit_class1/README.txt
# todo: investigar essa bomba https://zenodo.org/records/6983500

#= This code extracts the APN (Almost Perfect Nonlinear) matrix from a given S-box.
The function `extract_apn_matrix` takes an S-box represented as a vector of integers and
the number of bits `n` as input. It computes the APN matrix by iterating through pairs
of basis elements and applying the formula to calculate the entries of the matrix.
The resulting compact matrix is returned as a vector of integers.
The code also includes an example S-box and prints the extracted APN matrix. =#

# Reference: Weng et al. - On Quadratic Almost Perfect Nonlinear Functions and Their Related Algebraic Object


# todo: basicamente agor aeu ja tenho o codigo para gerar as tuplas; o que eu preciso fazer agora é
# implementar o algoritmo 1 do artigo, mas primeiro eles fazem um pre processamento, ou seja, eles aplicam as proposicoes 4 e 5
# para reduzir o espaço de busca, entao daria pra fazer isso primeiro, assim eu ja validaria pq eles dizem que eliminou x tuplas
#
#Para n=6: Elas cortam 8 das 17 tuplas instantaneamente.
#Para n=7: Elas excluem 13 das 27 tuplas.
#Para n=8: Elas eliminam 15 das 32 tuplas.
#
#
#De acordo com a Tabela 1 do artigo, ao passar as suas 17 tuplas pelos filtros, seu código deve apresentar exatamente este comportamento:
#
#    A Proposição 4 (que verifica se as dimensões dos subespaços afins caem nos tamanhos proibidos de 2, 4 ou n−1)
#   deve descartar automaticamente as Classes 6, 9, 13, 16 e 17.
#    A Proposição 5 (que verifica a existência de um quadrinômio que seja múltiplo dos polinômios mínimos
# das matrizes A e B) deve descartar automaticamente as Classes 4, 8 e 12.

# todo: tentar colocar em gray code as matrizes pra visualizar melhor as diferenças entre elas