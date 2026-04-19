using Nemo

function is_prime_int(n)
    n = Int(n)

    if n <= 1
        return false
    end

    if n == 2
        return true
    end

    if n % 2 == 0
        return false
    end

    d = 3

    while d * d <= n
        if n % d == 0
            return false
        end

        d += 2
    end

    return true
end

function matrix_multiplicative_order(A::FqMatrix)
    check_square(A)

    F = base_ring(A)
    n = nrows(A)
    I = identity_matrix(F, n)

    power = I

    # Safe enough for n = 6, 7, 8 in this context, but you can increase if needed.
    max_iterations = 10_000_000

    for k in 1:max_iterations
        power = power * A

        if power == I
            return k
        end
    end

    error("Could not determine multiplicative order within $max_iterations iterations")
end

function matrix_minimal_polynomial(A::FqMatrix)
    # First try the Sage-like name, if exported by your environment.
    if isdefined(Main, :minimal_polynomial)
        return getfield(Main, :minimal_polynomial)(A)
    end

    # Then try a shorter common CAS-style name.
    if isdefined(Main, :minpoly)
        return getfield(Main, :minpoly)(A)
    end

    error("""
    No minimal polynomial function was found.

    Try in the Julia REPL:

        methods(minpoly)
        methods(minimal_polynomial)

    Then replace matrix_minimal_polynomial(A) with the correct function for your Nemo version.
    """)
end

function matrix_is_similar(A::FqMatrix, B::FqMatrix)
    check_compatible_pair(A, B)

    if isdefined(Main, :is_similar)
        return getfield(Main, :is_similar)(A, B)
    end

    error("""
    No matrix similarity function was found.

    Try in the Julia REPL:

        methods(is_similar)

    If Nemo does not provide it directly for your matrix type, we need to implement similarity via rational canonical form invariants.
    """)
end