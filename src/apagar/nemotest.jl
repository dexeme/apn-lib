#= Finite Fields
Nemo package
https://github.com/Nemocas/Nemo.jl/tree/master/docs/src
Nemo.jl is a computer algebra package for Julia.
It provides a wide range of functionality for working with algebraic structures, including finite fields.

Usage:

finite_field(p::Int, n::Int, name::String) -> (R, x) # Create a finite field of order p^n with a generator named name
finite_field(3, 2, "a") -> (K, a) # Finite field of degree 2 and characteristic 3, a
finite_field(9, "a") -> (K, a) # Finite field of degree 2 and characteristic 3, a
=#

using Nemo

R, x = finite_field(7, 11, "x")
println("Finite field: ", R)
println("Generator: ", x)

# todo: completar com o que ta na wiki deles
a = R(3)
b = R(5)
println("a: ", a)
println("b: ", b)
println("a + b: ", add(a, b))
println("a * b: ", mul(a, b))
println("Inverse of a: ", inv(a))
println("a - b: ", sub(a, b))

