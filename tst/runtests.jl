using AlphaCertified
using DynamicPolynomials
using Test


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
