# Data from Table 1 of Kaleyski, "Deciding EA-equivalence via invariants",
# with test functions taken from Edel and Pott's APN catalogue.
#
# Nemo's GF(2, 6, "a") and GF(2, 8, "a") use the primitive polynomials
# required by Edel-Pott:
#   n = 6: x^6 + x^4 + x^3 + x + 1
#   n = 8: x^8 + x^4 + x^3 + x^2 + 1
#
# Keep this file as the readable source of truth. The LUTs used by tests are
# generated pointwise from these definitions by
# scripts/generate_kaleyski_table1_luts.jl. In particular, reference(...) means
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
    (8, "1.2") => 680,
    (8, "1.3") => 8,
    (8, "1.4") => 8,
    (8, "1.5") => 4,
    (8, "1.6") => 4,
    (8, "1.7") => 1,
    (8, "1.8") => 4,
    (8, "1.9") => 4,
    (8, "1.10") => 2,
    (8, "1.11") => 4,
    (8, "1.12") => 4,
    (8, "1.13") => 2,
    (8, "1.14") => 2,
    (8, "1.15") => 1,
    (8, "1.16") => 2,
    (8, "1.17") => 2,
    (8, "2.1") => 360,
    (8, "3.1") => 4,
    (8, "4.1") => 16,
    (8, "5.1") => 8,
    (8, "6.1") => 8,
    (8, "7.1") => 680,
)

abstract type KaleyskiTable1Component end

struct KaleyskiMonomial <: KaleyskiTable1Component
    coefficient_power::Int
    exponent::Int
end

struct KaleyskiTraceTerm
    coefficient_power::Int
    exponent::Int
end

struct KaleyskiAbsoluteTrace <: KaleyskiTable1Component
    scale_power::Int
    terms::Vector{KaleyskiTraceTerm}
end

struct KaleyskiRelativeTrace <: KaleyskiTable1Component
    scale_power::Int
    extension_degree::Int
    terms::Vector{KaleyskiTraceTerm}
end

struct KaleyskiReference <: KaleyskiTable1Component
    id::String
end

struct KaleyskiTable1Definition
    n::Int
    table_id::String
    equation_id::String
    formula::String
    components::Vector{KaleyskiTable1Component}
end

term(coefficient_power::Int, exponent::Int) = KaleyskiTraceTerm(coefficient_power, exponent)
xpow(exponent::Int) = KaleyskiMonomial(0, exponent)
monomial(coefficient_power::Int, exponent::Int) = KaleyskiMonomial(coefficient_power, exponent)
reference(id::AbstractString) = KaleyskiReference(String(id))

function absolute_trace(terms::KaleyskiTraceTerm...; scale::Int = 0)
    return KaleyskiAbsoluteTrace(scale, collect(terms))
end

function relative_trace(extension_degree::Int, terms::KaleyskiTraceTerm...; scale::Int = 0)
    return KaleyskiRelativeTrace(scale, extension_degree, collect(terms))
end

function entry(n::Int, id::AbstractString, formula::AbstractString, components::KaleyskiTable1Component...)
    return KaleyskiTable1Definition(n, String(id), String(id), String(formula), collect(components))
end

function entry(n::Int, table_id::AbstractString, equation_id::AbstractString, formula::AbstractString,
               components::KaleyskiTable1Component...)
    return KaleyskiTable1Definition(n, String(table_id), String(equation_id), String(formula), collect(components))
end

const KALEYSKI_TABLE1_DEFINITIONS = [
    entry(6, "1.1", "x^3",
        xpow(3)),
    entry(6, "1.2", "(No. 1.1) + u*(tr(u^56*x^3) + tr_8/2(u^18*x^9))",
        reference("1.1"),
        absolute_trace(term(56, 3); scale = 1),
        relative_trace(3, term(18, 9); scale = 1)),
    entry(6, "2.3", "2.1", "x^3 + u*x^24 + x^10",
        xpow(3), monomial(1, 24), xpow(10)),
    entry(6, "2.4", "2.2", "(No. 2.1) + u^3*(tr(u^10*x^3 + u^53*x^5) + tr_8/2(u^36*x^9))",
        reference("2.1"),
        absolute_trace(term(10, 3), term(53, 5); scale = 3),
        relative_trace(3, term(36, 9); scale = 3)),
    entry(6, "2.7", "2.3", "(No. 2.1) + tr(u^34*x^3 + u^48*x^5) + tr_8/2(u^9*x^9)",
        reference("2.1"),
        absolute_trace(term(34, 3), term(48, 5)),
        relative_trace(3, term(9, 9))),
    entry(6, "2.10", "2.4", "(No. 2.1) + u^2*(tr(u^24*x^3 + u^28*x^5) + tr_8/2(x^9))",
        reference("2.1"),
        absolute_trace(term(24, 3), term(28, 5); scale = 2),
        relative_trace(3, term(0, 9); scale = 2)),
    entry(6, "2.1", "2.5", "(No. 2.3) + u^42*(tr(u^10*x^3 + u^51*x^5) + tr_8/2(u^9*x^9))",
        reference("2.3"),
        absolute_trace(term(10, 3), term(51, 5); scale = 42),
        relative_trace(3, term(9, 9); scale = 42)),
    entry(6, "2.5", "2.6", "(No. 2.3) + u^23*(tr(u^31*x^3 + u^49*x^5) + tr_8/2(u^9*x^9))",
        reference("2.3"),
        absolute_trace(term(31, 3), term(49, 5); scale = 23),
        relative_trace(3, term(9, 9); scale = 23)),
    entry(6, "2.6", "2.7", "(No. 2.3) + u^12*(tr(u^42*x^3 + u^13*x^5) + tr_8/2(u^54*x^9))",
        reference("2.3"),
        absolute_trace(term(42, 3), term(13, 5); scale = 12),
        relative_trace(3, term(54, 9); scale = 12)),
    entry(6, "2.8", "(No. 2.3) + u*(tr(u^51*x^3 + u^60*x^5) + tr_8/2(u^18*x^9))",
        reference("2.3"),
        absolute_trace(term(51, 3), term(60, 5); scale = 1),
        relative_trace(3, term(18, 9); scale = 1)),
    entry(6, "2.9", "(No. 2.3) + u^14*(tr(u^18*x^3 + u^61*x^5) + tr_8/2(u^18*x^9))",
        reference("2.3"),
        absolute_trace(term(18, 3), term(61, 5); scale = 14),
        relative_trace(3, term(18, 9); scale = 14)),
    entry(6, "2.11", "2.10", "(No. 2.3) + u^17*tr(u^50*x^3 + u^56*x^5)",
        reference("2.3"),
        absolute_trace(term(50, 3), term(56, 5); scale = 17)),
    entry(6, "2.12", "2.11", "(No. 2.3) + u^19*(tr(u^11*x^3 + u^7*x^5 + u^38*x^7 + u^61*x^11 + u^23*x^13) + tr_8/2(u^54*x^9) + tr_4/2(u^42*x^21))",
        reference("2.3"),
        absolute_trace(term(11, 3), term(7, 5), term(38, 7), term(61, 11), term(23, 13); scale = 19),
        relative_trace(3, term(54, 9); scale = 19),
        relative_trace(2, term(42, 21); scale = 19)),
    entry(6, "2.2", "2.12", "(No. 2.4) + u*(tr(u^54*x^3 + u^47*x^5) + tr_8/2(u^9*x^9))",
        reference("2.4"),
        absolute_trace(term(54, 3), term(47, 5); scale = 1),
        relative_trace(3, term(9, 9); scale = 1)),

    entry(8, "1.1", "x^3",
        xpow(3)),
    entry(8, "1.2", "(No. 1.1) + tr(u^48*x^3 + x^9)",
        reference("1.1"), absolute_trace(term(48, 3), term(0, 9))),
    entry(8, "1.3", "(No. 1.1) + u*tr(u^63*x^3 + u^252*x^9)",
        reference("1.1"), absolute_trace(term(63, 3), term(252, 9); scale = 1)),
    entry(8, "1.4", "(No. 1.2) + u^38*tr(u^84*x^3 + u^213*x^9)",
        reference("1.2"), absolute_trace(term(84, 3), term(213, 9); scale = 38)),
    entry(8, "1.5", "(No. 1.2) + u^51*tr(u^253*x^3 + u^102*x^9)",
        reference("1.2"), absolute_trace(term(253, 3), term(102, 9); scale = 51)),
    entry(8, "1.6", "(No. 1.3) + u^154*tr(u^68*x^3 + u^235*x^9)",
        reference("1.3"), absolute_trace(term(68, 3), term(235, 9); scale = 154)),
    entry(8, "1.7", "(No. 1.4) + u^69*tr(u^147*x^3 + u^20*x^9)",
        reference("1.4"), absolute_trace(term(147, 3), term(20, 9); scale = 69)),
    entry(8, "1.8", "(No. 1.5) + u^68*tr(u^153*x^3 + u^51*x^9)",
        reference("1.5"), absolute_trace(term(153, 3), term(51, 9); scale = 68)),
    entry(8, "1.9", "(No. 1.6) + u^35*tr(u^216*x^3 + u^116*x^9)",
        reference("1.6"), absolute_trace(term(216, 3), term(116, 9); scale = 35)),
    entry(8, "1.10", "(No. 1.7) + u^22*tr(u^232*x^3 + u^195*x^9)",
        reference("1.7"), absolute_trace(term(232, 3), term(195, 9); scale = 22)),
    entry(8, "1.11", "(No. 1.8) + u^85*tr(u^243*x^3 + u^170*x^9)",
        reference("1.8"), absolute_trace(term(243, 3), term(170, 9); scale = 85)),
    entry(8, "1.12", "(No. 1.9) + u^103*tr(u^172*x^3 + u^31*x^9)",
        reference("1.9"), absolute_trace(term(172, 3), term(31, 9); scale = 103)),
    entry(8, "1.13", "(No. 1.10) + u^90*(tr(u^87*x^3 + u^141*x^5 + u^20*x^9) + tr_16/2(u^51*x^17))",
        reference("1.10"),
        absolute_trace(term(87, 3), term(141, 5), term(20, 9); scale = 90),
        relative_trace(4, term(51, 17); scale = 90)),
    entry(8, "1.14", "(No. 1.11) + u^5*tr(u^160*x^3 + u^250*x^9)",
        reference("1.11"), absolute_trace(term(160, 3), term(250, 9); scale = 5)),
    entry(8, "1.15", "(No. 1.11) + u^102*tr(u^6*x^3 + u^119*x^9)",
        reference("1.11"), absolute_trace(term(6, 3), term(119, 9); scale = 102)),
    entry(8, "1.16", "(No. 1.14) + u^64*tr(u^133*x^3 + u^30*x^9)",
        reference("1.14"), absolute_trace(term(133, 3), term(30, 9); scale = 64)),
    entry(8, "1.17", "(No. 1.16) + u^78*tr(u^235*x^3 + u^146*x^9)",
        reference("1.16"), absolute_trace(term(235, 3), term(146, 9); scale = 78)),
    entry(8, "2.1", "x^3 + x^17 + u^16*(x^18 + x^33) + u^15*x^48",
        xpow(3), xpow(17), monomial(16, 18), monomial(16, 33), monomial(15, 48)),
    entry(8, "3.1", "x^3 + u^24*x^6 + u^182*x^132 + u^67*x^192",
        xpow(3), monomial(24, 6), monomial(182, 132), monomial(67, 192)),
    entry(8, "4.1", "x^3 + x^6 + x^68 + x^80 + x^132 + x^160",
        xpow(3), xpow(6), xpow(68), xpow(80), xpow(132), xpow(160)),
    entry(8, "5.1", "x^3 + x^5 + x^18 + x^40 + x^66",
        xpow(3), xpow(5), xpow(18), xpow(40), xpow(66)),
    entry(8, "6.1", "x^3 + x^12 + x^40 + x^66 + x^130",
        xpow(3), xpow(12), xpow(40), xpow(66), xpow(130)),
    entry(8, "7.1", "x^57",
        xpow(57)),
]

const KALEYSKI_TABLE1_DEFINITION_BY_KEY = Dict((definition.n, definition.table_id) => definition
                                               for definition in KALEYSKI_TABLE1_DEFINITIONS)
const KALEYSKI_TABLE1_DEFINITION_BY_EQUATION_KEY = Dict((definition.n, definition.equation_id) => definition
                                                        for definition in KALEYSKI_TABLE1_DEFINITIONS)

const KALEYSKI_TABLE1_GENERATION_CASES = [(n = definition.n, id = definition.table_id)
                                          for definition in KALEYSKI_TABLE1_DEFINITIONS]

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
    return KALEYSKI_TABLE1_DEFINITION_BY_KEY[(n, String(id))].formula
end

function kaleyski_table1_equation_id(n::Int, id::AbstractString)
    return KALEYSKI_TABLE1_DEFINITION_BY_KEY[(n, String(id))].equation_id
end

function kaleyski_coefficient(field, exponent::Int)
    iszero(exponent) && return one(field)
    return gen(field)^exponent
end
