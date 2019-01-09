using AlphaCertified
using DynamicPolynomials
using Test

@test AlphaCertified begin
    @polyvar x y z
    f = (0.5 + 3im)*x*y^2*z
    @test matrix_form(f) == [1 2 1 1 2 3 1]
end
