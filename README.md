#AlphaCertified

This is a simple Julia wrapper for [alphaCertified](https://arxiv.org/abs/1011.1091). Loading the package provides two functions `certify` and `refine_points`. Both functions work with `BigFloat`. To set the precision, use `setprecision`. For instance, for working with 96 bits use `setprecision(96)`.


`certify` applies alphaCertified's certification algorithm to the system `F` and the solutions `solutions`.
```julia
certify(F, solutions;
    system_file=nothing,
    points_file=nothing,
    settings_file=nothing,
    rationalize=false,
    dir = mktempdir(), kwargs...)
```
Here, `F` is an element or vector of type `<:MultivariatePolynomials.AbstractPolynomialLike`.

`refine_points` applies the Newton operator to `F` and `solutions` to refine the solutions.
```julia
refine_points(F, solutions; numiterations=3)
```
