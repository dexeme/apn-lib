tuple_search_representative(n::Int, index::Int)::APNRepresentative =
    APNRepresentative(Symbol("tuple_search_$index"),
                      :PrecomputedTupleSearch,
                      n,
                      () -> precomputed_tuple_search(n, index),
                      metadata = Dict(:index => index))

const APN_REPRESENTATIVES_BY_DIMENSION = Dict{Int, Vector{APNRepresentative}}(
    3 => APNRepresentative[
        tuple_search_representative(3, 1),
        tuple_search_representative(3, 3),
        APNRepresentative(:gold_i1, :Gold, 3, APNFunction(3, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 3, family_c4(3)),
        APNRepresentative(:c5, :C5, 3, family_c5(3)),
        APNRepresentative(:c6, :C6, 3, family_c6(3)),
    ],
    4 => APNRepresentative[
        APNRepresentative(:gold_i1, :Gold, 4, APNFunction(4, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 4, family_c4(4)),
    ],
    5 => APNRepresentative[
        APNRepresentative(:gold_i1, :Gold, 5, APNFunction(5, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 5, family_c4(5)),
    ],
    6 => APNRepresentative[
        tuple_search_representative(6, 5),
        APNRepresentative(:gold_i1, :Gold, 6, APNFunction(6, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 6, family_c4(6)),
        APNRepresentative(:c5, :C5, 6, family_c5(6)),
        APNRepresentative(:c6, :C6, 6, family_c6(6)),
    ],
    7 => APNRepresentative[
        tuple_search_representative(7, 4),
        tuple_search_representative(7, 5),
        tuple_search_representative(7, 7),
        tuple_search_representative(7, 8),
        tuple_search_representative(7, 9),
        tuple_search_representative(7, 10),
        APNRepresentative(:gold_i1, :Gold, 7, APNFunction(7, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 7, family_c4(7)),
    ],
    8 => APNRepresentative[
        APNRepresentative(:gold_i1, :Gold, 8, APNFunction(8, monomial_expr(3)), metadata = Dict(:i => 1)),
        APNRepresentative(:c4, :C4, 8, family_c4(8)),
    ],
)
