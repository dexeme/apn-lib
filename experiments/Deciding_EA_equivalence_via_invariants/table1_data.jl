# Data from Table 1 of Kaleyski, "Deciding EA-equivalence via invariants",
# with test functions taken from Edel and Pott's APN catalogue.
#@article{Edel:200900,
# doi = {10.3934/amc.2009.3.59},
# author = {Yves Edel and Alexander Pott},
# title = {{A new almost perfect nonlinear function which is not quadratic}},
# year = 2009,
# journal = {{Advances in Mathematics of Communications}},
# volume = 3,
# number = 1,
# pages = {59--81},
#}
# Keep this file as the readable source of truth. The LUTs used by tests are
# generated pointwise from these definitions by
# experiments/Deciding_EA_equivalence_via_invariants/generate_luts.jl. In particular, reference(...) means
# "use the already generated LUT value for this input and add the perturbation".

const KALEYSKI_TABLE1_EXPECTED = Dict{Tuple{Int, String}, Int}(
    (6, "1.1") => 1008,
    (6, "1.2") => 336,
    (6, "2.1") => 10,
    (6, "2.2") => 336,
    (6, "2.3") => 1008,
    (6, "2.4") => 8,
    (6, "2.5") => 60,
    (6, "2.6") => 8,
    (6, "2.7") => 10,
    (6, "2.8") => 8,
    (6, "2.9") => 8,
    (6, "2.10") => 8,
    (6, "2.11") => 8,
    (6, "2.12") => 48,
    (8, "1.1") => 680,
    (8, "1.2") => 8,    # Kaleyski alocou na linha 1.3
    (8, "1.3") => 4,    # Kaleyski alocou na linha 1.5
    (8, "1.4") => 1,    # A 1ª 'new', Kaleyski alocou na linha 1.7
    (8, "1.5") => 4,    # A 2ª 'new', Kaleyski alocou na linha 1.8
    (8, "1.6") => 4,    # A 3ª 'new', Kaleyski alocou na linha 1.9
    (8, "1.7") => 2,    # A 4ª 'new', Kaleyski alocou na linha 1.10
    (8, "1.8") => 4,    # A 5ª 'new', Kaleyski alocou na linha 1.11
    (8, "1.9") => 4,    # A 6ª 'new', Kaleyski alocou na linha 1.12
    (8, "1.10") => 2,   # A 7ª 'new', Kaleyski alocou na linha 1.13
    (8, "1.11") => 8,   # Kaleyski alocou na linha 1.4
    (8, "1.12") => 4,   # Kaleyski alocou na linha 1.6
    (8, "1.13") => 2,   # A 8ª 'new', Kaleyski alocou na linha 1.14
    (8, "1.14") => 1,   # A 9ª 'new', Kaleyski alocou na linha 1.15
    (8, "1.15") => 680, # Kaleyski alocou na linha 1.2
    (8, "1.16") => 2,   # A 10ª 'new', Kaleyski alocou na linha 1.16
    (8, "1.17") => 2,   # A 11ª 'new', Kaleyski alocou na linha 1.17
    (8, "2.1") => 360,
    (8, "3.1") => 4,
    (8, "4.1") => 16,
    (8, "5.1") => 8,
    (8, "6.1") => 8,
    (8, "7.1") => 680,
)

const KALEYSKI_TABLE1_CATALOGUE = Catalogue(
    APNFunction(6, "1.1",
        monomial_expr(3)),
    APNFunction(6, "1.2",
        reference("1.1"),
        absolute_trace(trace_term(56, 3); scale = 1),
        relative_trace(3, trace_term(18, 9); scale = 1)),
    APNFunction(6, "2.3",
        monomial_expr(3), monomial_expr(1, 24; base = :u), monomial_expr(10)),
    APNFunction(6, "2.4",
        reference("2.1"),
        absolute_trace(trace_term(10, 3), trace_term(53, 5); scale = 3),
        relative_trace(3, trace_term(36, 9); scale = 3)),
    APNFunction(6, "2.7",
        reference("2.1"),
        absolute_trace(trace_term(34, 3), trace_term(48, 5)),
        relative_trace(3, trace_term(9, 9))),
    APNFunction(6, "2.10",
        reference("2.1"),
        absolute_trace(trace_term(24, 3), trace_term(28, 5); scale = 2),
        relative_trace(3, trace_term(0, 9); scale = 2)),
    APNFunction(6, "2.1",
        reference("2.3"),
        absolute_trace(trace_term(10, 3), trace_term(51, 5); scale = 42),
        relative_trace(3, trace_term(9, 9); scale = 42)),
    APNFunction(6, "2.5",
        reference("2.3"),
        absolute_trace(trace_term(31, 3), trace_term(49, 5); scale = 23),
        relative_trace(3, trace_term(9, 9); scale = 23)),
    APNFunction(6, "2.6",
        reference("2.3"),
        absolute_trace(trace_term(42, 3), trace_term(13, 5); scale = 12),
        relative_trace(3, trace_term(54, 9); scale = 12)),
    APNFunction(6, "2.8",
        reference("2.3"),
        absolute_trace(trace_term(51, 3), trace_term(60, 5); scale = 1),
        relative_trace(3, trace_term(18, 9); scale = 1)),
    APNFunction(6, "2.9",
        reference("2.3"),
        absolute_trace(trace_term(18, 3), trace_term(61, 5); scale = 14),
        relative_trace(3, trace_term(18, 9); scale = 14)),
    APNFunction(6, "2.11",
        reference("2.3"),
        absolute_trace(trace_term(50, 3), trace_term(56, 5); scale = 17)),
    APNFunction(6, "2.12",
        reference("2.3"),
        absolute_trace(trace_term(11, 3), trace_term(7, 5), trace_term(38, 7), trace_term(61, 11), trace_term(23, 13); scale = 19),
        relative_trace(3, trace_term(54, 9); scale = 19),
        relative_trace(2, trace_term(42, 21); scale = 19)),
    APNFunction(6, "2.2",
        reference("2.4"),
        absolute_trace(trace_term(54, 3), trace_term(47, 5); scale = 1),
        relative_trace(3, trace_term(9, 9); scale = 1)),

    APNFunction(8, "1.1",
        monomial_expr(3)),
    APNFunction(8, "1.2",
        reference("1.1"), absolute_trace(trace_term(48, 3), trace_term(0, 9))),
    APNFunction(8, "1.3",
        reference("1.1"), absolute_trace(trace_term(63, 3), trace_term(252, 9); scale = 1)),
    APNFunction(8, "1.4",
        reference("1.2"), absolute_trace(trace_term(84, 3), trace_term(213, 9); scale = 38)),
    APNFunction(8, "1.5",
        reference("1.2"), absolute_trace(trace_term(253, 3), trace_term(102, 9); scale = 51)),
    APNFunction(8, "1.6",
        reference("1.3"), absolute_trace(trace_term(68, 3), trace_term(235, 9); scale = 154)),
    APNFunction(8, "1.7",
        reference("1.4"), absolute_trace(trace_term(147, 3), trace_term(20, 9); scale = 69)),
    APNFunction(8, "1.8",
        reference("1.5"), absolute_trace(trace_term(153, 3), trace_term(51, 9); scale = 68)),
    APNFunction(8, "1.9",
        reference("1.6"), absolute_trace(trace_term(216, 3), trace_term(116, 9); scale = 35)),
    APNFunction(8, "1.10",
        reference("1.7"), absolute_trace(trace_term(232, 3), trace_term(195, 9); scale = 22)),
    APNFunction(8, "1.11",
        reference("1.8"), absolute_trace(trace_term(243, 3), trace_term(170, 9); scale = 85)),
    APNFunction(8, "1.12",
        reference("1.9"), absolute_trace(trace_term(172, 3), trace_term(31, 9); scale = 103)),
    APNFunction(8, "1.13",
        reference("1.10"), absolute_trace(trace_term(87, 3), trace_term(141, 5), trace_term(20, 9); scale = 90),
        relative_trace(4, trace_term(51, 17); scale = 90)),
    APNFunction(8, "1.14",
        reference("1.11"), absolute_trace(trace_term(160, 3), trace_term(250, 9); scale = 5)),
    APNFunction(8, "1.15",
        reference("1.11"), absolute_trace(trace_term(6, 3), trace_term(119, 9); scale = 102)),
    APNFunction(8, "1.16",
        reference("1.14"), absolute_trace(trace_term(133, 3), trace_term(30, 9); scale = 64)),
    APNFunction(8, "1.17",
        reference("1.16"), absolute_trace(trace_term(235, 3), trace_term(146, 9); scale = 78)),
    APNFunction(8, "2.1",
        monomial_expr(3), monomial_expr(17), monomial_expr(16, 18; base = :u), monomial_expr(16, 33; base = :u), monomial_expr(15, 48; base = :u)),
    APNFunction(8, "3.1",
        monomial_expr(3), monomial_expr(24, 6; base = :u), monomial_expr(182, 132; base = :u), monomial_expr(67, 192; base = :u)),
    APNFunction(8, "4.1",
        monomial_expr(3), monomial_expr(6), monomial_expr(68), monomial_expr(80), monomial_expr(132), monomial_expr(160)),
    APNFunction(8, "5.1",
        monomial_expr(3), monomial_expr(5), monomial_expr(18), monomial_expr(40), monomial_expr(66)),
    APNFunction(8, "6.1",
        monomial_expr(3), monomial_expr(12), monomial_expr(40), monomial_expr(66), monomial_expr(130)),
    APNFunction(8, "7.1",
        monomial_expr(57)),
)

const KALEYSKI_TABLE1_FUNCTIONS = KALEYSKI_TABLE1_CATALOGUE.functions

const KALEYSKI_TABLE1_FUNCTION_BY_KEY = Dict((function_.n, function_.id) => function_
                                             for function_ in KALEYSKI_TABLE1_FUNCTIONS)

const KALEYSKI_TABLE1_EQUATION_IDS = Dict{Tuple{Int, String}, String}(
    (6, "1.1") => "1.1",
    (6, "1.2") => "1.2",
    (6, "2.3") => "2.1",
    (6, "2.4") => "2.2",
    (6, "2.7") => "2.3",
    (6, "2.10") => "2.4",
    (6, "2.1") => "2.5",
    (6, "2.5") => "2.6",
    (6, "2.6") => "2.7",
    (6, "2.8") => "2.8",
    (6, "2.9") => "2.9",
    (6, "2.11") => "2.10",
    (6, "2.12") => "2.11",
    (6, "2.2") => "2.12",
)

const KALEYSKI_TABLE1_FUNCTION_BY_EQUATION_KEY = Dict((function_.n, get(KALEYSKI_TABLE1_EQUATION_IDS, (function_.n, function_.id), function_.id)) => function_
                                                      for function_ in KALEYSKI_TABLE1_FUNCTIONS)

const KALEYSKI_TABLE1_GENERATION_CASES = [(n = function_.n, id = function_.id)
                                          for function_ in KALEYSKI_TABLE1_FUNCTIONS]

const KALEYSKI_TABLE1_CASES = [
    (n = 6, id = "1.1"),
    (n = 6, id = "1.2"),
    (n = 6, id = "2.1"),
    (n = 6, id = "2.2"),
    (n = 6, id = "2.3"),
    (n = 6, id = "2.4"),
    (n = 6, id = "2.5"),
    (n = 6, id = "2.6"),
    (n = 6, id = "2.7"),
    (n = 6, id = "2.8"),
    (n = 6, id = "2.9"),
    (n = 6, id = "2.10"),
    (n = 6, id = "2.11"),
    (n = 6, id = "2.12"),
    (n = 8, id = "1.1"),
    (n = 8, id = "1.2"),
    (n = 8, id = "1.3"),
    (n = 8, id = "1.4"),
    (n = 8, id = "1.5"),
    (n = 8, id = "1.6"),
    (n = 8, id = "1.7"),
    (n = 8, id = "1.8"),
    (n = 8, id = "1.9"),
    (n = 8, id = "1.10"),
    (n = 8, id = "1.11"),
    (n = 8, id = "1.12"),
    (n = 8, id = "1.13"),
    (n = 8, id = "1.14"),
    (n = 8, id = "1.15"),
    (n = 8, id = "1.16"),
    (n = 8, id = "1.17"),
    (n = 8, id = "2.1"),
    (n = 8, id = "3.1"),
    (n = 8, id = "4.1"),
    (n = 8, id = "5.1"),
    (n = 8, id = "6.1"),
    (n = 8, id = "7.1"),
]

function kaleyski_table1_expected_permutations(n::Int, id::AbstractString)
    return KALEYSKI_TABLE1_EXPECTED[(n, String(id))]
end

function kaleyski_table1_formula(n::Int, id::AbstractString)
    return string(KALEYSKI_TABLE1_FUNCTION_BY_KEY[(n, String(id))])
end

function kaleyski_table1_equation_id(n::Int, id::AbstractString)
    key = (n, String(id))
    return get(KALEYSKI_TABLE1_EQUATION_IDS, key, key[2])
end

function kaleyski_coefficient(field, exponent::Int)
    iszero(exponent) && return one(field)
    return gen(field)^exponent
end
