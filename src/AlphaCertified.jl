module AlphaCertified

export certify

import MultivariatePolynomials
const MP = MultivariatePolynomials

Base.rationalize(i::Int64) = i // 1
function matrix_form(f::MP.AbstractPolynomialLike, variables=MP.variables(f))
    nvars = length(variables)
    nterms = MP.nterms(f)
    M = zeros(Int, nterms, nvars + 4)
    for (i, term) in enumerate(MP.terms(f))
        for (j, v) in enumerate(variables)
            M[i, j] = MP.degree(term, v)
        end
        re, im = rationalize.(reim(MP.coefficient(term)))
        M[i, nvars+1] = numerator(re)
        M[i, nvars+2] = denominator(re)
        M[i, nvars+3] = numerator(im)
        M[i, nvars+4] = denominator(im)
    end
    M
end

function matrix_form(io::IO, F::Vector{<:MP.AbstractPolynomialLike})
    vars = MP.variables(F)
    npolys = length(F)
    nvars = length(vars)
    # header
    println(io, nvars, " ", npolys)
    # systems block
    println(io)
    for f in F
        M = matrix_form(f, vars)
        nterms = size(M, 1)
        println(io, nterms)
        for i ∈ 1:nterms
            for j ∈ 1:nvars
                print(io, M[i, j], " ")
            end
            # real coeff part
            print(io, M[i, nvars+1], "/", M[i, nvars+2], " ")
            print(io, M[i, nvars+3], "/", M[i, nvars+4], "\n")
        end
    end
    nothing
end

function input_points(io::IO, points::Vector{<:Vector{<:Number}})
    println(io, length(points))
    for p in points
        println(io, "\n")
        for xᵢ in p
            re, im = rationalize.(reim(xᵢ))
            print(io, numerator(re), "/", denominator(re))
            print(io, " ")
            print(io, numerator(im), "/", denominator(im))
            print(io, "\n")
        end
    end
    nothing
end


"""
    certify(F, solutions; dir=mktempdir())

"""
function certify(F::Vector{<:MP.AbstractPolynomialLike}, solutions;
    dir = mktempdir())
    polySys = sprint(matrix_form, F)
    println("File directory: ", dir)
    write(joinpath(dir, "polySys"), polySys)

    points = sprint(input_points, solutions)
    write(joinpath(dir, "points"), points)

    old_dir = pwd()
    cd(dir)
    run(`alphaCertified polySys points`)
    cd(old_dir)
end

end # module
