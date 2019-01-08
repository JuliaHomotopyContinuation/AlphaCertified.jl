using AlphaCertified
using DynamicPolynomials
using Test
using HomotopyContinuation
using PolynomialTestSystems


@polyvar x y z

f = (0.5 + 3im)*x*y^2*z
@test matrix_form(f) == [1 2 1 1 2 3 1]

g = [x^2-3*x*y+z^4,
    (0.5 + 3im)*x*y^2*z - 9z + 2x + 7,
    z^10-im]

points = [[2/3+2im, -7, 3im], [1/8, -7+2im, -1-0.4im]]


g = [x^2-1, y^2-1]
points = [[0.9999, 0.99999]]
certify(g, points)
# @test matrix_form(g) ==

AlphaCertified.matrix_form!(stdout, g)

exponents(f)


rationalize(2.3+5im)


F = equations(katsura(5))
r = solve(F)

certify(F, solutions(r),
    rationalize=false,
    NEWTONONLY=1, ARITHMETICTYPE=1,
    NUMITERATIONS=5)

using DelimitedFiles


entries = readdlm("/var/folders/d8/9j36bh6j5190wk69s6qwpf140000gp/T/tmpkMn1Cq/refinedPoints", ' ', String, skipstart=1, skipblanks=false)
