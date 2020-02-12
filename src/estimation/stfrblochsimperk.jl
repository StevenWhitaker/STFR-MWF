function stfrblochsimperk(
    y::AbstractArray{<:Real,2}, # [D,N]
    noise::AbstractArray{<:Real,1}; # noise voxels (i.e., background)
    T::Integer = 20000,
    H::Integer = 4000,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    σ::Real = sqrt(sum(noise.^2) / (2 * length(noise))),
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    M0DistTrain = nothing,
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
    ΔωDefaultDist = Uniform(-50 * 2π, 50 * 2π),
    κDefaultDist = Uniform(0.8, 1.2),
    Δω = ΔωDefaultDist,
    κ = κDefaultDist,
    resetRNG::Bool = true
)

    # Reset the RNG for reproducible results
    if resetRNG
        Random.seed!(20190823)
    end

    # Check whether Δω and κ are known or unknown
    Δωknown = Δω isa AbstractArray
    κknown  = κ  isa AbstractArray

    # Load the scan design
    scanDesigns = load(scanDesignFile, "P")

    # Set up the signal model
    signalModels = getstfrblochsimsignalModels(scanDesigns, Δωknown, κknown)

    # Create distributions for the known parameters for training
    ΔωDistTrain = Δωknown ? FitDist(Δω, 0.1) : Δω
    κDistTrain  =  κknown ? FitDist(κ, 0.01) : κ

    # Determine M0DistTrain, unless provided by user
    if isnothing(M0DistTrain)

        M0DistTrain = getstfrblochsimM0DistTrain(ffDistTrain, fmDistTrain,
            T1fDistTrain, T1sDistTrain, T1mDistTrain, T2fDistTrain,
            T2sDistTrain, T2mDistTrain, τfsDistTrain, τfmDistTrain,
            ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δωknown, κknown, maximum(y),
            signalModels)

    end

    # Collect known parameters and training priors
    if Δωknown
        if κknown
            ν = [Δω, κ]
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain]
            νDists = [ΔωDistTrain, κDistTrain]
        else
            ν = [Δω]
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      κDistTrain]
            νDists = [ΔωDistTrain]
        end
    else
        if κknown
            ν = [κ]
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      ΔωDistTrain]
            νDists = [κDistTrain]
        else
            ν = Array{Array{Real,1},1}()
            xDists = [M0DistTrain, ffDistTrain, fmDistTrain,
                      T1fDistTrain, T1sDistTrain, T1mDistTrain,
                      T2fDistTrain, T2sDistTrain, T2mDistTrain,
                      τfsDistTrain, τfmDistTrain, ΔωfDistTrain,
                      ΔωDistTrain, κDistTrain]
            νDists = Array{Array{Float64,1},1}()
        end
    end

    # Generate length scales
    Λ = λ * max.(dropdims(mean(abs.(y), dims = 2), dims = 2), eps())
    if Δωknown
        # Use 10 Hz instead of eps() to avoid scaling issues when testing with
        # Δω fixed to 0 (or close to 0)
        push!(Λ, λ * max(mean(abs.(Δω)), 10 * 2π))
    end
    if κknown
        push!(Λ, λ * max(mean(abs.(κ)), eps()))
    end

    # Specify the kernel
    kernel = GaussianRFF(H, Λ)

    noiseDist = Normal(0, σ)

    # Run PERK
    (xhat, trainData, ttrain, ttest) = perk(y, ν, T, xDists, νDists, noiseDist,
                                            signalModels, kernel, ρ)

    # Return PERK estimates, training data, and timing info
    return (xhat, trainData, ttrain, ttest)

end

function getstfrblochsimsignalModels(scanDesigns, Δωknown, κknown)

    (Tfree, Tg, α, β, ϕ, TE) = scanDesigns[1]

	signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ) ->
		stfrblochsim.(M0, [[ff, 1-ff-fm, fm]],
					  [[T1f, T1s, T1m]], [[T2f, T2s, T2m]],
					  [[Δω+Δωf, Δω, Δω]],
					  [[τfs, τfm, τfs*(1-ff-fm)/ff, Inf, Inf, Inf]],
					  κ, Tfree, Tg, α, β, ϕ, TE)]
    if Δωknown && !κknown
		tmp = signalModels[1]
		signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, κ, Δω) ->
						tmp(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ)]
    end

    return signalModels

end

# See Gopal's MATLAB code here (lines 65-73)
# https://github.com/gopal-nataraj/perk/blob/master/map/t1-t2/perk_train.m
function getstfrblochsimM0DistTrain(ffDistTrain, fmDistTrain, T1fDistTrain,
    T1sDistTrain, T1mDistTrain, T2fDistTrain, T2sDistTrain, T2mDistTrain,
    τfsDistTrain, τfmDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δωknown,
    κknown, ymax, signalModels)

    N   = 10000
    M0  = 1
    ff  = rand(ffDistTrain,  N)
    fm  = rand(fmDistTrain,  N)
    T1f = rand(T1fDistTrain, N)
    T1s = rand(T1sDistTrain, N)
    T1m = rand(T1mDistTrain, N)
    T2f = rand(T2fDistTrain, N)
    T2s = rand(T2sDistTrain, N)
    T2m = rand(T2mDistTrain, N)
    τfs = rand(τfsDistTrain, N)
    τfm = rand(τfmDistTrain, N)
    Δωf = rand(ΔωfDistTrain, N)
    Δω  = rand(ΔωDistTrain,  N)
    κ   = rand(κDistTrain,   N)
    if Δωknown
        unknown = [M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf]
        known = [Δω]
    else
        unknown = [M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω]
        known = []
    end
    if κknown
        push!(known, κ)
    else
        push!(unknown, κ)
    end

    y = reduce(vcat, [reduce(hcat, signalModels[i].(unknown..., known...))
                        for i = 1:length(signalModels)])
    y = abs.(y)

    M0DistTrain = Uniform(eps(), ymax / mean(y))

    return M0DistTrain

end
