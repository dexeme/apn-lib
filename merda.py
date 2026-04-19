from sage.all import GF, PolynomialRing, matrix, vector
import sboxU

F = GF(2 ** 7, name='g')
g = F.gen()
R = PolynomialRing(F, name='x')
x = R.gen()

A_gold = matrix(F, 7, 7)
for i in range(7):
    for j in range(7):
        if i != j:
            A_gold[i, j] = g ** (2 * i + j) + g ** (2 * j + i)


def int_to_elem(n):
    return sum([F((n >> k) & 1) * (g ** k) for k in range(7)])


def elem_to_int(elem):
    try:
        return int(elem.integer_representation())
    except AttributeError:
        return sum([int(c) * (1 << k) for k, c in enumerate(list(elem))])


iteracao = 0

while True:
    iteracao += 1
    if iteracao % 100 == 0:
        print(f"Buscando... {iteracao} matrizes testadas.")

    B = matrix(F, 7, 7)
    for i in range(7):
        for j in range(i + 1, 7):
            if i < 2 or j < 2:
                B[i, j] = A_gold[i, j]
            else:
                B[i, j] = F.random_element()

    A = B + B.transpose()
    matriz_valida = True

    for linha in range(2, 7):
        elementos_linha = [A[linha, c] for c in range(7) if c != linha]

        matriz_F2 = matrix(GF(2), [
            [(elem_to_int(e) >> k) & 1 for k in range(7)]
            for e in elementos_linha
        ])

        if matriz_F2.rank() != 6:
            matriz_valida = False
            break

    if not matriz_valida:
        continue

    lut = []
    for i in range(128):
        v_vec = vector(F, [F((i >> k) & 1) for k in range(7)])
        val = v_vec * B * v_vec
        lut.append(elem_to_int(val))

    if sboxU.differential_uniformity(lut) == 2:
        print(f"\nFuncao APN encontrada na iteracao {iteracao}!")
        print("Formatando o polinomio univariado...\n")

        pontos = [(int_to_elem(i), int_to_elem(lut[i])) for i in range(128)]
        poly = R.lagrange_polynomial(pontos)

        terms = []
        for i, c in enumerate(poly.list()):
            if c != 0:
                if c == 1:
                    terms.append(f"x^{i}")
                else:
                    pot = c.log(g)
                    terms.append(f"g{pot}x^{i}")

        print(" + ".join(terms))
        break