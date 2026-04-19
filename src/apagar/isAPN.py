# Consider a 2-bit mapping to illustrate the distinction.
#
# **The S-box**
# The S-box is the data structure. It is the lookup table that executes the mapping $F: \mathbb{F}_2^2 \to \mathbb{F}_2^2$. Let us define an S-box $S$:
#
# $S = [0, 1, 3, 2]$
#
# This table maps inputs to outputs:
# *   Input 0 $\to$ Output 0
# *   Input 1 $\to$ Output 1
# *   Input 2 $\to$ Output 3
# *   Input 3 $\to$ Output 2
#
# **The APN Property**
# The APN property is a mathematical test applied to that structure. To check if $S$ is APN, we compute the discrete derivative $S(x \oplus a) \oplus S(x)$ for a non-zero input difference $a$.
#
# Let $a = 1$. We evaluate the derivative for all possible values of $x$:
# *   $x = 0: S(0 \oplus 1) \oplus S(0) = S(1) \oplus S(0) = 1 \oplus 0 = 1$
# *   $x = 1: S(1 \oplus 1) \oplus S(1) = S(0) \oplus S(1) = 0 \oplus 1 = 1$
# *   $x = 2: S(2 \oplus 1) \oplus S(2) = S(3) \oplus S(2) = 2 \oplus 3 = 1$
# *   $x = 3: S(3 \oplus 1) \oplus S(3) = S(2) \oplus S(3) = 3 \oplus 2 = 1$
#
# For the input difference $a = 1$, the resulting output difference $b = 1$ occurs 4 times.
#
# A function is only APN if the maximum number of solutions for any derivative equation $S(x \oplus a) \oplus S(x) = b$ is exactly 2. Because this S-box produced 4 solutions for a single differential pair, $S$ is a valid S-box, but it is not an APN function.
#
# def is_apn(func_table):
#     field_size = len(func_table)
#     for a in range(1, field_size):
#         diff_values = set()
#         for x in range(field_size):
#             res = func_table[x ^ a] ^ func_table[x]
#             if res in diff_values:
#                 continue
#             diff_values.add(res)
#
#         if len(diff_values) != (field_size // 2):
#             return False
#     return True

# f_table = [0, 1, 3, 6, 7, 4, 5, 2]
# f1_table = [2, 5, 3, 1, 0, 4, 7, 6]
# is_valid_apn = is_apn(f_table)
# print(is_valid_apn)
# is_valid_apn2 = is_apn(f1_table)
# print(is_valid_apn2)

# f_3_table = [0,0,0,4,0,10,30,16,0,36,78,110,26,52,74,96,0,228,2,226,98,140,126,148,13,205,65,133,117,191,39,233,0,48,112,68,17,43,127,65,18,6,44,60,25,7,57,35,40,252,90,138,91,133,55,237,55,199,11,255,94,164,124,130,0,56,237,209,76,126,191,137,180,168,23,15,226,244,95,77,236,48,3,219,194,20,51,225,85,173,244,8,97,147,222,40,132,140,25,21,217,219,90,92,34,14,241,217,101,67,168,138,64,172,223,55,127,153,254,28,235,35,58,246,206,12,1,199,0,151,198,85,190,35,102,255,68,247,204,123,224,89,118,203,212,167,16,103,8,113,210,175,157,202,23,68,91,6,207,150,164,3,18,177,11,166,163,10,242,113,10,141,71,206,161,44,88,27,236,171,149,220,63,114,3,100,249,154,212,185,48,89,246,89,221,118,4,161,49,144,6,141,99,236,238,111,149,16,206,133,231,168,94,31,105,44,51,92,84,63,185,220,192,161,214,73,141,22,53,160,112,225,52,143,33,158,205,124,198,115,198,189,159,224,71,54,0,117,41,118,62,101,178,231,187,234]
# is_valid_apn_3 = is_apn(f_3_table)
# print(is_valid_apn_3)

# TODO: OBJETIVO: https://zenodo.org/records/16752428
# PEGAR AS ENTRADAS DESSE SITE, E CONVERTER PRA FORMA MATRICIAL COMPACTA DESCRITA EM
# Weng et al. - On Quadratic Almost Perfect Nonlinear Functions and Their Related Algebraic Object
# NESSE SITE, TEM A QUANTIDADE DE FUNCOES EM CADA TRABALHO, OBJETIVO É SEPARAR QUAL É QUAL E CATALOGAR EM CODIGO
# A FORMA MATRICIAL (SE POSSIVEL), ASSIM TENDO UM DATASET MUITO MAIS COMPACTO

