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
ALL_TUPLES_7_1_TUPLE = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127]
ALL_TUPLES_7_1_SEARCH = [0, 1, 2, 5, 4, 7, 10, 6, 8, 21, 14, 25, 20, 34, 12, 47, 16, 82, 42, 59, 28, 23, 50, 90, 40, 91, 68, 102, 24, 52, 94, 39, 32, 97, 37, 72, 84, 101, 118, 13, 56, 57, 46, 62, 100, 121, 53, 125, 80, 19, 55, 108, 9, 103, 77, 60, 48, 85, 104, 99, 61, 119, 78, 49, 64, 66, 67, 3, 74, 76, 17, 87, 41, 93, 75, 45, 109, 51, 26, 83, 112, 36, 114, 70, 92, 31, 124, 126, 73, 54, 115, 30, 106, 113, 123, 88, 33, 65, 38, 107, 110, 86, 89, 105, 18, 35, 79, 63, 27, 15, 120, 44, 96, 117, 43, 116, 81, 95, 71, 22, 122, 58, 111, 11, 29, 69, 98, 127]
ALL_TUPLES_7_4_SEARCH = [0, 1, 2, 107, 4, 39, 105, 126, 8, 84, 78, 82, 109, 34, 67, 58, 16, 44, 23, 59, 35, 72, 27, 60, 101, 31, 68, 13, 57, 9, 116, 19, 32, 119, 88, 74, 46, 104, 118, 14, 70, 98, 47, 127, 54, 99, 120, 125, 117, 10, 62, 93, 55, 71, 26, 50, 114, 40, 18, 11, 87, 21, 38, 112, 64, 53, 81, 75, 15, 33, 43, 115, 92, 77, 111, 48, 83, 103, 28, 45, 51, 61, 123, 12, 94, 3, 65, 6, 108, 73, 121, 24, 79, 102, 69, 122, 85, 56, 20, 90, 124, 25, 5, 113, 110, 97, 49, 96, 52, 7, 100, 37, 91, 86, 80, 89, 36, 30, 22, 66, 17, 29, 42, 41, 76, 63, 95, 106]

def is_apn(func_table):
    field_size = len(func_table)
    for a in range(1, field_size):
        diff_values = set()
        for x in range(field_size):
            res = func_table[x ^ a] ^ func_table[x]
            if res in diff_values:
                continue
            diff_values.add(res)

        if len(diff_values) != (field_size // 2):
            return False
    return True

is_valid_apn = is_apn(ALL_TUPLES_7_4_SEARCH)
print(is_valid_apn)
is_valid_apn2 = is_apn(ALL_TUPLES_7_1_SEARCH)
print(is_valid_apn2)

# f_3_table = [0,0,0,4,0,10,30,16,0,36,78,110,26,52,74,96,0,228,2,226,98,140,126,148,13,205,65,133,117,191,39,233,0,48,112,68,17,43,127,65,18,6,44,60,25,7,57,35,40,252,90,138,91,133,55,237,55,199,11,255,94,164,124,130,0,56,237,209,76,126,191,137,180,168,23,15,226,244,95,77,236,48,3,219,194,20,51,225,85,173,244,8,97,147,222,40,132,140,25,21,217,219,90,92,34,14,241,217,101,67,168,138,64,172,223,55,127,153,254,28,235,35,58,246,206,12,1,199,0,151,198,85,190,35,102,255,68,247,204,123,224,89,118,203,212,167,16,103,8,113,210,175,157,202,23,68,91,6,207,150,164,3,18,177,11,166,163,10,242,113,10,141,71,206,161,44,88,27,236,171,149,220,63,114,3,100,249,154,212,185,48,89,246,89,221,118,4,161,49,144,6,141,99,236,238,111,149,16,206,133,231,168,94,31,105,44,51,92,84,63,185,220,192,161,214,73,141,22,53,160,112,225,52,143,33,158,205,124,198,115,198,189,159,224,71,54,0,117,41,118,62,101,178,231,187,234]
# is_valid_apn_3 = is_apn(f_3_table)
# print(is_valid_apn_3)

# TODO: target dataset: https://zenodo.org/records/16752428
# Download the entries from this site and convert them to the compact matrix form described in
# Weng et al. - On Quadratic Almost Perfect Nonlinear Functions and Their Related Algebraic Object
# The site lists the number of functions in each work; the goal is to separate and catalog each source in code,
# preferably in matrix form, so the dataset is much more compact.
