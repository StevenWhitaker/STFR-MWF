"""
    stfrblochsimperk_bias(y, noise, scanDesigns; T, H, ρ, λ, σ, M0DistTrain,
        ffDistTrain, fmDistTrain, T1fDistTrain, T1sDistTrain, T1mDistTrain,
        T2fDistTrain, T2sDistTrain, T2mDistTrain, τfsDistTrain, τfmDistTrain,
        ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δω, κ, useFitDist, resetRNG)

Estimate unknown parameters using STFR scans and training with a
three-compartment exchanging tissue model.

# Arguments
- `y::AbstractArray{<:Real,2}`: STFR scan data \\[D,N\\] (D is number of scans,
  and N is number of voxels or data points per scan)
- `noise::AbstractArray{<:Real,1}`: Noise/background voxels (for estimating the
  noise level to add for training)
- `scanDesigns::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1}`:
  STFR scan designs

# Options
Note: Typically, these will not need to be modified.
- `T::Integer = 20000`: Number of training data points to simulate
- `H::Integer = 4000`: Size of approximation of training Gram matrix
- `ρ::Real = 2^-60`: Ridge regression regularization parameter
- `λ::Real = 2^3.5`: Length scales regularization parameter
- `σ::Real = ...`: Standard deviation of noise to add to training data; default
  estimates it from `noise`
- `M0DistTrain = nothing`: M0 training distribution; if set to `nothing`
  (default) it will be determined using the STFR signal model and the observed
  data `y`
- `ffDistTrain = Uniform(0.03, 0.31)`: ff (MWF) training distribution
- `fmDistTrain = Uniform(0.03, 0.31)`: fm (volume fraction of macromolecules)
  training distribution
- `T1fDistTrain = Uniform(320, 480)`: T1f (T1 of myelin water) training
  distribution (ms)
- `T1sDistTrain = Uniform(800, 1200)`: T1s (T1 of non-myelin water) training
  distribution (ms)
- `T1mDistTrain = Uniform(800, 3000)`: T1m (T1 of macromolecules) training
  distribution (ms)
- `T2fDistTrain = Uniform(16, 24)`: T2f (T2 of myelin water) training
  distribution (ms)
- `T2sDistTrain = Uniform(64, 96)`: T2s (T2 of non-myelin water) training
  distribution (ms)
- `T2mDistTrain = Uniform(0.01, 0.1)`: T2m (T2 of macromolecules) training
  distribution (ms)
- `τfsDistTrain = Uniform(80, 150)`: τfs (residence time of exchange from myelin
  water to non-myelin water) training distribution (ms)
- `τfmDistTrain = Uniform(40, 75)`: τfm (residence time of exchange from myelin
  water to macromolecules) training distribution (ms)
- `ΔωfDistTrain = Uniform(0, 35 * 2π)`: Δωf (additional off-resonance
  experienced by myelin water) training distribution (rad/s)
- `ΔωDistTrain = Uniform(-50 * 2π, 50 * 2π)`: Default training distribution
  for Δω (bulk off-resonance) (rad/s); used if `useFitDist` is `false` or `Δω`
  is unknown
- `κDistTrain = Uniform(0.8, 1.2)`: Default training distribution for κ
  (flip-angle scaling constant); used if `useFitDist` is `false` or `κ`
  is unknown
- `Δω = nothing`: Vector of length N (i.e., `size(y,2)`) containing the known
  Δω values for each voxel, if known (unknown by default)
- `κ = nothing`: Vector of length N (i.e., `size(y,2)`) containing the known κ
  values for each voxel, if known (unknown by default)
- `useFitDist::Bool = true`: Whether to fit training distributions for the
  known parameters to the known parameter values; if `false`, use the default
  training distribution
- `resetRNG::Bool = true`: Whether to reset the random number generator

# Return
- `xhat::Array{Float64,2}`: Estimates of all unknown parameters \\[L,N\\] (L is
  number of unknown parameters)
- `trainData::TrainingData`: Training data object that can be passed to `perk`
  directly; typically this should be discarded because right now the training
  data uses the actual data `y`, so is subject dependent
- `ttrain::Float64`: Amount of time training took (s)
- `ttest::Float64`: Amount of time testing took (s)
"""
function stfrblochsimperk_bias(
    y::AbstractArray{<:Real,2}, # [D,N]
    noise::AbstractArray{<:Real,1}, # noise voxels (i.e., background)
    scanDesigns::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    T::Integer = 20000,
    H::Integer = 4000,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    σ::Real = sqrt(sum(noise.^2) / (2 * length(noise))),
    M0DistTrain = nothing, # Replace with distribution if sure of possible M0 values
    ffDistTrain = Uniform(0.03, 0.31),
    fmDistTrain = Uniform(0.03, 0.31),
    T1fDistTrain = Uniform(320, 480),
    T1sDistTrain = Uniform(800, 1200),
    T1mDistTrain = Uniform(800, 3000),
    T2fDistTrain = Uniform(16, 24),
    T2sDistTrain = Uniform(64, 96),
    T2mDistTrain = Uniform(0.01, 0.1),
    τfsDistTrain = Uniform(80, 150), # Ballpark value
    τfmDistTrain = Uniform(40, 75), # About twice as fast as τfs
    ΔωfDistTrain = Uniform(0 * 2π, 35 * 2π),
    ΔωDistTrain = Uniform(-50 * 2π, 50 * 2π),
    κDistTrain = Uniform(0.8, 1.2),
    Δω = nothing, # Set to array of values if known
    κ = nothing, # Set to array of values if known
    useFitDist::Bool = true,
    resetRNG::Bool = true
)

    # Reset the RNG for reproducible results
    if resetRNG
        Random.seed!(20191007)
    end

    # Check whether Δω and κ are known or unknown
    Δωknown = Δω isa AbstractArray
    κknown  = κ  isa AbstractArray

    # Set up the signal model
    signalModels = getstfrblochsimsignalModels_bias(scanDesigns, Δωknown, κknown)

    # Create distributions for the known parameters for training
    ΔωDistTrain = (Δωknown && useFitDist) ? FitDist(Δω, 0.1) : ΔωDistTrain
    κDistTrain  =  (κknown && useFitDist) ? FitDist(κ, 0.01) : κDistTrain

    # Determine M0DistTrain, unless provided by user
    if isnothing(M0DistTrain)

        params = randparams(N = 10000, [1], ffDistTrain, fmDistTrain,
            T1fDistTrain, T1sDistTrain, T1mDistTrain, T2fDistTrain,
            T2sDistTrain, T2mDistTrain, τfsDistTrain, τfmDistTrain,
            ΔωfDistTrain, ΔωDistTrain, κDistTrain)[2:end] # Exclude M0

        M0DistTrain = getstfrblochsimM0DistTrain_bias(params..., Δωknown, κknown,
            maximum(y), signalModels)

    end

    # Collect training priors
    (xDists, νDists) = collectdists(M0DistTrain, ffDistTrain, fmDistTrain,
        T1fDistTrain, T1sDistTrain, T1mDistTrain, T2fDistTrain, T2sDistTrain,
        T2mDistTrain, τfsDistTrain, τfmDistTrain, ΔωfDistTrain, ΔωDistTrain,
        κDistTrain, Δωknown, κknown)

    # Generate training data
    xtrain = rand.(xDists, T) # [L][T]
    νtrain = rand.(νDists, T) # [K][T]

    # Specify training noise
    noiseDist = Normal(0, σ)
    noisetrain = complex.(rand(noiseDist, size(y, 1), T),
                          rand(noiseDist, size(y, 1), T))

    # Random Fourier features
    f = randn(H, size(y, 1) + length(νtrain))
    ph = rand(H)

    # Set Δω and κ to nothing if unknown
    Δω = Δωknown ? Δω : nothing
    κ  = κknown  ? κ  : nothing

    # Figure out ΔωTrain and κTrain
    if Δωknown
        ΔωTrain = νtrain[1]
    elseif κknown
        ΔωTrain = xtrain[end]
    else
        ΔωTrain = xtrain[end-1]
    end
    κTrain = κknown ? νtrain[end] : xtrain[end]

    return stfrblochsimperk_bias(y, scanDesigns, Δω, κ, noisetrain, nothing,
        xtrain[1:12]..., ΔωTrain, κTrain, f, ph, ρ = ρ, λ = λ)

end

function stfrblochsimperk_bias(
    y::AbstractArray{<:Real,2}, # [D,N]
    scanDesigns::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1},
    Δω::Union{<:AbstractArray{<:Real,1},Nothing}, # [N]
    κ::Union{<:AbstractArray{<:Real,1},Nothing}, # [N]
    noisetrain::AbstractArray{<:Complex}, # [D,T]
    M0params, # Output of randparams (excluding first entry), or nothing
    M0Train::AbstractArray{<:Real,1}, # [T]
    ffTrain::AbstractArray{<:Real,1}, # [T]
    fmTrain::AbstractArray{<:Real,1}, # [T]
    T1fTrain::AbstractArray{<:Real,1}, # [T]
    T1sTrain::AbstractArray{<:Real,1}, # [T]
    T1mTrain::AbstractArray{<:Real,1}, # [T]
    T2fTrain::AbstractArray{<:Real,1}, # [T]
    T2sTrain::AbstractArray{<:Real,1}, # [T]
    T2mTrain::AbstractArray{<:Real,1}, # [T]
    τfsTrain::AbstractArray{<:Real,1}, # [T]
    τfmTrain::AbstractArray{<:Real,1}, # [T]
    ΔωfTrain::AbstractArray{<:Real,1}, # [T]
    ΔωTrain::AbstractArray{<:Real,1}, # [T]
    κTrain::AbstractArray{<:Real,1}, # [T]
    f::AbstractArray{<:Real,2}, # RFF frequency (unscaled) [H,D+K]
    ph::AbstractArray{<:Real,1}; # RFF phase [H]
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5
)

    # Check whether Δω and κ are known or unknown
    Δωknown = !isnothing(Δω)
    κknown  = !isnothing(κ)

    # Set up the signal model
    signalModels = getstfrblochsimsignalModels_bias(scanDesigns, Δωknown, κknown)

    # Scale M0Train to use updated M0DistTrain without actually regenerating
    # random samples
    if !isnothing(M0params)

        M0DistTrain = getstfrblochsimM0DistTrain_bias(M0params..., Δωknown, κknown,
            maximum(y), signalModels)

        M0Train = M0Train .* (mean(M0DistTrain) / mean(M0Train))

    end

    # Collect known parameters and specify minimum length scales
    # The purpose for specifying minimum length scales is to avoid scaling by
    # eps() when testing with Δω fixed to (close to) 0
    # I suppose that using similar logic would mandate that I avoid eps() for
    # all the length scales, but the other parameters won't ever be close
    # enough to 0 to matter
    if Δωknown
        if κknown
            ν = [Δω, κ]
            minscales = [fill(eps(), size(y, 1)); 10 * 2π; eps()]
        else
            ν = [Δω]
            minscales = [fill(eps(), size(y, 1)); 10 * 2π]
        end
    elseif κknown
        ν = [κ]
        minscales = fill(eps(), size(y, 1) + 1)
    else
        ν = Array{Array{Float64,1},1}()
        minscales = fill(eps(), size(y, 1))
    end

    # Generate length scales
    if isempty(ν)
        q = y # [D+K,N], K = 0
    else
        q = [y; transpose(hcat(ν...))] # [D+K,N]
    end
    Λ = λ * max.([mean(abs.(q[d,:])) for d = 1:size(q, 1)], minscales)

    # Specify the kernel
    kernel = GaussianRFF(length(ph), Λ)

    # Generate training data
    (xtrain, νtrain) = collectdists(M0Train, ffTrain, fmTrain, T1fTrain,
        T1sTrain, T1mTrain, T2fTrain, T2sTrain, T2mTrain, τfsTrain, τfmTrain,
        ΔωfTrain, ΔωTrain, κTrain, Δωknown, κknown)
    ytrain = vcat([hcat(signalModels[i].(xtrain..., νtrain...)...)
        for i = 1:length(signalModels)]...) # [D,T]

    # Add noise and take magnitude
    ytrain = abs.(ytrain .+ noisetrain)

    # Train
    ttrain = @elapsed trainData = train(ytrain, xtrain, νtrain, kernel, f, ph)

    # Run PERK
    (xhat, ttest) = perk(y, ν, trainData, kernel, ρ)

    # Return PERK estimates, training data, and timing info
    return (xhat, trainData, ttrain, ttest)

end

"""
    getstfrblochsimsignalModels_bias(scanDesigns, Δωknown, κknown)

Set up the signal models to be passed to `perk`.
"""
function getstfrblochsimsignalModels_bias(scanDesigns, Δωknown, κknown)

    (Tfree, Tg, α, β, ϕ, TE) = scanDesigns[1]

    if Δωknown && !κknown

        signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, κ, Δω) ->
            stfrblochsim.(M0, [[ff, 1-ff-fm, fm]],
                          [[T1f, T1s, T1m]], [[T2f, T2s, T2m]],
                          [[Δω+Δωf, Δω, Δω]],
                          [[τfs, τfm, τfs*(1-ff-fm)/ff, Inf, Inf, Inf]],
                          κ, Tfree, Tg, α, β, ϕ, TE)]

    else

        signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ) ->
            stfrblochsim.(M0, [[ff, 1-ff-fm, fm]],
                          [[T1f, T1s, T1m]], [[T2f, T2s, T2m]],
                          [[Δω+Δωf, Δω, Δω]],
                          [[τfs, τfm, τfs*(1-ff-fm)/ff, Inf, Inf, Inf]],
                          κ, Tfree, Tg, α, β, ϕ, TE)]

    end

    return signalModels

end

"""
    getstfrblochsimM0DistTrain_bias(ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm,
        Δωf, Δω, κ, Δωknown, κknown, ymax, signalModels)

Determine an appropriate training distribution for M0.

An appropriate distribution is chosen by simulating several noiseless
evaluations of the signal models over a variety of randomly chosen parameters,
but fixing M0 to 1 (the other parameters should be randomly chosen and then
passed to this function). Then, the maximum training M0 value is equal to the
maximum observed signal value divided by the mean simulated signal value,
while the minimum training M0 value is set to 0.
"""
function getstfrblochsimM0DistTrain_bias(ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs,
    τfm, Δωf, Δω, κ, Δωknown, κknown, ymax, signalModels)

    M0 = 1
    (unknown, known) = collectdists(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m,
                                    τfs, τfm, Δωf, Δω, κ, Δωknown, κknown)

    y = vcat([hcat(signalModels[i].(unknown..., known...)...)
              for i = 1:length(signalModels)]...)
    y = abs.(y)

    M0DistTrain = Uniform(eps(), ymax / mean(y))

    return M0DistTrain

end

"""
    randparams(M0Dist, ffDist, fmDist, T1fDist, T1sDist, T1mDist, T2fDist,
        T2sDist, T2mDist, τfsDist, τfmDist, ΔωfDist, ΔωDist, κDist; N)

Randomly generate parameters for a three-compartment, exchanging model.
"""
function randparams(M0Dist, ffDist, fmDist, T1fDist, T1sDist, T1mDist, T2fDist,
    T2sDist, T2mDist, τfsDist, τfmDist, ΔωfDist, ΔωDist, κDist; N = 10000)

    M0  = rand(M0Dist,  N)
    ff  = rand(ffDist,  N)
    fm  = rand(fmDist,  N)
    T1f = rand(T1fDist, N)
    T1s = rand(T1sDist, N)
    T1m = rand(T1mDist, N)
    T2f = rand(T2fDist, N)
    T2s = rand(T2sDist, N)
    T2m = rand(T2mDist, N)
    τfs = rand(τfsDist, N)
    τfm = rand(τfmDist, N)
    Δωf = rand(ΔωfDist, N)
    Δω  = rand(ΔωDist,  N)
    κ   = rand(κDist,   N)

    return (M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ)

end

"""
    collectdists(M0DistTrain, ffDistTrain, fmDistTrain, T1fDistTrain,
        T1sDistTrain, T1mDistTrain, T2fDistTrain, T2sDistTrain, T2mDistTrain,
        τfsDistTrain, τfmDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain,
        Δωknown, κknown)

Collect training distributions into lists of known and unknown parameters.
"""
function collectdists(M0DistTrain, ffDistTrain, fmDistTrain, T1fDistTrain,
    T1sDistTrain, T1mDistTrain, T2fDistTrain, T2sDistTrain, T2mDistTrain,
    τfsDistTrain, τfmDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δωknown,
    κknown)

    if Δωknown
        if κknown
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain]
            νDists = [ΔωDistTrain, κDistTrain]
        else
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      κDistTrain]
            νDists = [ΔωDistTrain]
        end
    else
        if κknown
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      ΔωDistTrain]
            νDists = [κDistTrain]
        else
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      ΔωDistTrain, κDistTrain]
            νDists = Array{Array{Float64,1},1}()
        end
    end

    return (xDists, νDists)

end
