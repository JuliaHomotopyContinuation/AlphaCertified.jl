module AlphaCertified

export certify, refine_points

import MultivariatePolynomials
using DelimitedFiles
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

function matrix_form(io::IO, F::AbstractVector{<:MP.AbstractPolynomialLike})
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

function input_points(io::IO, points::Vector{<:AbstractVector{<:Number}}, convert_to_rational=true)
    println(io, length(points))
    for p in points
        println(io, "\n")
        for xᵢ in p
            if convert_to_rational
                re, im = rationalize_re_im(xᵢ)
                print(io, numerator(re), "/", denominator(re))
                print(io, " ")
                print(io, numerator(im), "/", denominator(im))
                print(io, "\n")
            else
                re, im = reim(xᵢ)
                print(io, re)
                print(io, " ")
                print(io, im)
                print(io, "\n")
            end
        end
    end
    nothing
end

rationalize_re_im(x::Complex{BigFloat}) = rationalize.(BigInt, reim(x))
rationalize_re_im(x::BigFloat) = rationalize.(BigInt, reim(x))
rationalize_re_im(x) = rationalize.(reim(x))

function write_settings(io::IO, settings)
    for (k, v) in settings
        println(io, k, ": ", v)
    end
end

"""
    certify(F, solutions; dir=mktempdir())

"""
function certify(F::AbstractVector{<:MP.AbstractPolynomialLike}, solutions;
    system_file=nothing,
    points_file=nothing,
    settings_file=nothing,
    rationalize=true,
    dir = mktempdir(), kwargs...)

    if system_file === nothing
        polySys = sprint(matrix_form, F)
        println("File directory: ", dir)
        write(joinpath(dir, "polySys"), polySys)
        system_file = "polySys"
    end

    if points_file === nothing
        points = sprint(input_points, solutions, rationalize)
        write(joinpath(dir, "points"), points)
        points_file = "points"
    end

    if !isempty(kwargs)
        settings = sprint(write_settings, kwargs)
        write(joinpath(dir, "customSettings"), settings)
        settings_file = "customSettings"
    end
    old_dir = pwd()
    try
        cd(dir)
        if settings_file === nothing
            run(`alphaCertified polySys $(points_file)`)
        else
            run(`alphaCertified polySys $(points_file) $(settings_file)`)
        end
    finally
        cd(old_dir)
    end
    dir
end

function refine_points(F::AbstractVector{<:MP.AbstractPolynomialLike}, solutions; numiterations=3)
    dir = certify(F, solutions,
        rationalize = false,
        NEWTONONLY = 1,
        ARITHMETICTYPE = 1,
        NUMITERATIONS = numiterations)
    parse_refined_points(dir)
end

function parse_refined_points(dir)
    R = readdlm(joinpath(dir, "refinedPoints"), ' ', String, skipstart = 0)
    l = size(R,1) - 1
    s = parse(Int, R[1,1])
    r = Int(l/s)
    map(1:s) do i
        map((i-1)*r+2:i*r+1) do k
            complex(parse(BigFloat, R[k,1]), parse(BigFloat, R[k,2]))
         end
    end
end



end # module
