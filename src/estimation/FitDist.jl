"""
    FitDist(samples, resolution)

Create an object that generates samples according to the distribution
(histogram) of `samples`. Samples are generated using inverse transform
sampling.

# Arguments
- `samples`: Data points to which to fit a distribution
- `resolution::Real`: Spacing of consecutive points in the fitted distribution

# Return
- `f::FitDist`: FitDist object representing the distribution given by `samples`

# Examples
```
julia> f = FitDist(randn(10000), 0.01);

julia> v = rand(f, 10000); # Histogram of v is approximately standard normal
```
"""
struct FitDist{T}
    x::Array{T,1}   # Values in support
    cdf::Array{T,1} # Cumulative density at each value of x

    FitDist(samples, resolution) = begin
        # Generate the values in the support
        x   = collect(minimum(samples):resolution:maximum(samples))
        # Estimate the cdf from the samples
        cdf = ecdf(samples)
        # Create the FitDist object
        new{eltype(x)}(x, cdf.(x))
    end
end
Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{FitDist{T}}) where {T} =
begin
    d = s[]
    tmp = findlast(d.cdf .<= rand())
    d.x[isnothing(tmp) ? 1 : tmp]
end
Base.eltype(::Type{FitDist{T}}) where {T} = T
