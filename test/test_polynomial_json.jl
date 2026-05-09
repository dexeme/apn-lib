using Test
using Nemo
using APNLib

@testset "Polynomial JSON representation" begin
    json_str = """
    {
      "field": "GF(2^7)",
      "basis": "power",
      "modulus": "t^7+t+1",
      "terms": [
        {"exp": 3, "coef": 22},
        {"exp": 9, "coef": 5}
      ]
    }
    """

    polynomial = build_polynomial_from_json(json_str)
    field = base_ring(parent(polynomial))

    @test coeff(polynomial, 3) == int_to_field_element(22, field, 7)
    @test coeff(polynomial, 9) == int_to_field_element(5, field, 7)
    @test coeff(polynomial, 4) == zero(field)
end
