function stfrperk(
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
    T1fDistTrain = Uniform(320, 480),
    T1sDistTrain = Uniform(800, 1200),
    T2fDistTrain = Uniform(16, 24),
    T2sDistTrain = Uniform(64, 96),
    ΔωfDistTrain = Uniform(0 * 2π, 35 * 2π),
    ΔωDefaultDist = Uniform(-50 * 2π, 50 * 2π),
    κDefaultDist = Uniform(0.8, 1.2),
    Δω = ΔωDefaultDist,
    κ = κDefaultDist,
    ignoreΔωf::Bool = false,
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
    signalModels = getstfrsignalModels(scanDesigns, Δωknown, κknown, ignoreΔωf)

    # Create distributions for the known parameters for training
    ΔωDistTrain = Δωknown ? FitDist(Δω, 0.1) : Δω
    κDistTrain  =  κknown ? FitDist(κ, 0.01) : κ

    # Determine M0DistTrain, unless provided by user
    if isnothing(M0DistTrain)

        M0DistTrain = getstfrM0DistTrain(ffDistTrain, T1fDistTrain,
            T1sDistTrain, T2fDistTrain, T2sDistTrain, ΔωfDistTrain,
            ΔωDistTrain, κDistTrain, Δωknown, κknown, maximum(y), signalModels,
            ignoreΔωf)

    end

    # Collect known parameters and training priors
    if ignoreΔωf
        if Δωknown
            if κknown
                ν = [Δω, κ]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain]
                νDists = [ΔωDistTrain, κDistTrain]
            else
                ν = [Δω]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, κDistTrain]
                νDists = [ΔωDistTrain]
            end
        else
            if κknown
                ν = [κ]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωDistTrain]
                νDists = [κDistTrain]
            else
                ν = Array{Array{Real,1},1}()
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωDistTrain, κDistTrain]
                νDists = []
            end
        end
    else
        if Δωknown
            if κknown
                ν = [Δω, κ]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωfDistTrain]
                νDists = [ΔωDistTrain, κDistTrain]
            else
                ν = [Δω]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωfDistTrain,
                         κDistTrain]
                νDists = [ΔωDistTrain]
            end
        else
            if κknown
                ν = [κ]
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωfDistTrain, ΔωDistTrain]
                νDists = [κDistTrain]
            else
                ν = Array{Array{Real,1},1}()
                xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                         T2fDistTrain, T2sDistTrain, ΔωfDistTrain, ΔωDistTrain,
                         κDistTrain]
                νDists = Array{Array{Float64,1},1}()
            end
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

function getstfrsignalModels(scanDesigns, Δωknown, κknown, ignoreΔωf)

    (Tfree, Tg, α, β, ϕ, TE) = scanDesigns[1]

    if ignoreΔωf
		signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δω, κ) ->
						stfr.(M0, ff, T1f, T1s, T2f, T2s, 0, Δω, κ,
							  Tfree, Tg, α, β, ϕ, TE)]
        if Δωknown && !κknown
			tmp = signalModels[1]
			signalModels = [(M0, ff, T1f, T1s, T2f, T2s, κ, Δω) ->
							tmp(M0, ff, T1f, T1s, T2f, T2s, Δω, κ)]
        end
    else
		signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ) ->
						stfr.(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ,
							  Tfree, Tg, α, β, ϕ, TE)]
        if Δωknown && !κknown
			tmp = signalModels[1]
			signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, κ, Δω) ->
							tmp(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ)]
        end
    end

    return signalModels

end

# See Gopal's MATLAB code here (lines 65-73)
# https://github.com/gopal-nataraj/perk/blob/master/map/t1-t2/perk_train.m
function getstfrM0DistTrain(ffDistTrain, T1fDistTrain, T1sDistTrain,
	T2fDistTrain, T2sDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δωknown,
	κknown, ymax, signalModels, ignoreΔωf)

    N   = 10000
    M0  = 1
    ff  = rand(ffDistTrain,  N)
    T1f = rand(T1fDistTrain, N)
    T1s = rand(T1sDistTrain, N)
    T2f = rand(T2fDistTrain, N)
    T2s = rand(T2sDistTrain, N)
    Δωf = rand(ΔωfDistTrain, N)
    Δω  = rand(ΔωDistTrain,  N)
    κ   = rand(κDistTrain,   N)
    if Δωknown
        if ignoreΔωf
            unknown = [M0, ff, T1f, T1s, T2f, T2s]
        else
            unknown = [M0, ff, T1f, T1s, T2f, T2s, Δωf]
        end
        known = [Δω]
    else
        if ignoreΔωf
            unknown = [M0, ff, T1f, T1s, T2f, T2s, Δω]
        else
            unknown = [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω]
        end
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
