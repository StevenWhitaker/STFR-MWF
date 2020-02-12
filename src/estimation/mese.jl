function meseblochsimperk(
    y::AbstractArray{<:Real,2}, # [D,N]
    noise::AbstractArray{<:Real,1}; # noise voxels (i.e., background)
    T::Integer = 20000,
    H::Integer = 4000,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    σ = sqrt(sum(noise.^2) / (2 * length(noise))),
    TR::Real = 1200,
    αex::Real = π/2,
    αinv::Real = π,
    TE::Real = 10,
    nechoes::Integer = 32,
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
    nspins::Integer = 30, # Set to 1 to use ideal spoiling
    ncycles::Real = 1,
    resetRNG::Bool = true
)

    # Reset the RNG for reproducible results
    if resetRNG
        Random.seed!(20190823)
    end

    # Check whether Δω and κ are known or unknown
    Δωknown = Δω isa AbstractArray
    κknown  = κ  isa AbstractArray

    # Set up the scan design
    scanDesigns = [Real[TR, αex, αinv, TE, nechoes]]

    # Set up the signal model
    if nspins == 1
        signalModels = getmesesignalModels(scanDesigns, Δωknown, κknown)
    else
        signalModels = getmesesignalModels(scanDesigns, Δωknown, κknown, nspins,
            ncycles)
    end

    # Create distributions for the known parameters for training
    ΔωDistTrain = Δωknown ? FitDist(Δω, 0.1) : Δω
    κDistTrain  = κknown  ? FitDist(κ, 0.01) : κ

    # Determine M0DistTrain, unless provided by user
    if isnothing(M0DistTrain)

        M0DistTrain = getmeseM0DistTrain(ffDistTrain, fmDistTrain,
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
        push!(Λ, λ * max(mean(abs.(Δω)), eps()))
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

function getmesesignalModels(scanDesigns, Δωknown, κknown)

    (TR, αex, αinv, TE, nechoes) = scanDesigns[1]

    signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ) ->
        meseblochsim(M0, [ff, 1-ff-fm, fm], [T1f, T1s, T1m],
                     [T2f, T2s, T2m], [Δω+Δωf, Δω, Δω],
                     [τfs, τfm, τfs*(1-ff-fm)/ff, Inf, Inf, Inf], κ, TR,
                     αex, αinv, TE, nechoes)]
    if Δωknown && !κknown
        tmp = signalModels[1]
        signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, κ, Δω) ->
                        tmp(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ)]
    end

    return signalModels

end

function getmesesignalModels(scanDesigns, Δωknown, κknown, nspins, ncycles)

    (TR, αex, αinv, TE, nechoes) = scanDesigns[1]

    signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ) ->
        meseblochsim(M0, [ff, 1-ff-fm, fm], [T1f, T1s, T1m],
                     [T2f, T2s, T2m], [Δω+Δωf, Δω, Δω],
                     [τfs, τfm, τfs*(1-ff-fm)/ff, Inf, Inf, Inf], κ, TR,
                     αex, αinv, TE, nechoes, nspins, ncycles)]
    if Δωknown && !κknown
        tmp = signalModels[1]
        signalModels = [(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, κ, Δω) ->
                        tmp(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ)]
    end

    return signalModels

end

# See Gopal's MATLAB code here (lines 65-73)
# https://github.com/gopal-nataraj/perk/blob/master/map/t1-t2/perk_train.m
function getmeseM0DistTrain(ffDistTrain, fmDistTrain, T1fDistTrain,
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

function meseblochsim(M0, frac, T1, T2, Δω, τ, κ, TR, αex, αinv, TE, nechoes)

    spin = SpinMC(M0, frac, T1, T2, Δω / 2π, τ)
    return mese!(spin, TR, TE, nechoes, αex = κ * αex, αinv = κ * αinv)

end

function meseblochsim(M0, frac, T1, T2, Δω, τ, κ, TR, αex, αinv, TE, nechoes,
                      nspins, ncycles)

    spin = SpinMC(M0, frac, T1, T2, Δω / 2π, τ)
    return mese(spin, TR, TE, nechoes, αex = κ * αex, αinv = κ * αinv,
                nspins = nspins, ncycles = ncycles)

end

function mese2compblochsimperk(
    y::AbstractArray{<:Real,2}, # [D,N]
    noise::AbstractArray{<:Real,1}; # noise voxels (i.e., background)
    T::Integer = 20000,
    H::Integer = 4000,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    σ = sqrt(sum(noise.^2) / (2 * length(noise))),
    TR::Real = 1200,
    αex::Real = π/2,
    αinv::Real = π,
    TE::Real = 10,
    nechoes::Integer = 32,
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
    nspins::Integer = 30, # Set to 1 to use ideal spoiling
    ncycles::Real = 1,
    resetRNG::Bool = true
)

    # Reset the RNG for reproducible results
    if resetRNG
        Random.seed!(20190823)
    end

    # Check whether Δω and κ are known or unknown
    Δωknown = Δω isa AbstractArray
    κknown  = κ  isa AbstractArray

    # Set up the scan design
    scanDesigns = [Real[TR, αex, αinv, TE, nechoes]]

    # Set up the signal model
    if nspins == 1
        signalModels = getmese2compsignalModels(scanDesigns, Δωknown, κknown)
    else
        signalModels = getmese2compsignalModels(scanDesigns, Δωknown, κknown,
            nspins, ncycles)
    end

    # Create distributions for the known parameters for training
    ΔωDistTrain = Δωknown ? FitDist(Δω, 0.1) : Δω
    κDistTrain  = κknown  ? FitDist(κ, 0.01) : κ

    # Determine M0DistTrain, unless provided by user
    if isnothing(M0DistTrain)

        M0DistTrain = getmese2compM0DistTrain(ffDistTrain, T1fDistTrain,
            T1sDistTrain, T2fDistTrain, T2sDistTrain, ΔωfDistTrain, ΔωDistTrain,
            κDistTrain, Δωknown, κknown, maximum(y), signalModels)

    end

    # Collect known parameters and training priors
    if Δωknown
        if κknown
            ν = [Δω, κ]
            xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                      T2fDistTrain, T2sDistTrain, ΔωfDistTrain]
            νDists = [ΔωDistTrain, κDistTrain]
        else
            ν = [Δω]
            xDists = [M0DistTrain, ffDistTrain, T1fDistTrain, T1sDistTrain,
                      T2fDistTrain, T2sDistTrain, ΔωfDistTrain, κDistTrain]
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

    # Generate length scales
    Λ = λ * max.(dropdims(mean(abs.(y), dims = 2), dims = 2), eps())
    if Δωknown
        push!(Λ, λ * max(mean(abs.(Δω)), eps()))
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

function getmese2compsignalModels(scanDesigns, Δωknown, κknown)

    (TR, αex, αinv, TE, nechoes) = scanDesigns[1]

    signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ) ->
        meseblochsim(M0, [ff, 1-ff], [T1f, T1s],
                     [T2f, T2s], [Δω+Δωf, Δω],
                     [Inf, Inf], κ, TR,
                     αex, αinv, TE, nechoes)]
    if Δωknown && !κknown
        tmp = signalModels[1]
        signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, κ, Δω) ->
                        tmp(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ)]
    end

    return signalModels

end

function getmese2compsignalModels(scanDesigns, Δωknown, κknown, nspins, ncycles)

    (TR, αex, αinv, TE, nechoes) = scanDesigns[1]

    signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ) ->
        meseblochsim(M0, [ff, 1-ff], [T1f, T1s],
                     [T2f, T2s], [Δω+Δωf, Δω],
                     [Inf, Inf], κ, TR,
                     αex, αinv, TE, nechoes, nspins, ncycles)]
    if Δωknown && !κknown
        tmp = signalModels[1]
        signalModels = [(M0, ff, T1f, T1s, T2f, T2s, Δωf, κ, Δω) ->
                        tmp(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ)]
    end

    return signalModels

end

# See Gopal's MATLAB code here (lines 65-73)
# https://github.com/gopal-nataraj/perk/blob/master/map/t1-t2/perk_train.m
function getmese2compM0DistTrain(ffDistTrain, T1fDistTrain, T1sDistTrain,
    T2fDistTrain, T2sDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain, Δωknown,
    κknown, ymax, signalModels)

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
        unknown = [M0, ff, T1f, T1s, T2f, T2s, Δωf]
        known = [Δω]
    else
        unknown = [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω]
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

# Estimate MWF for many voxels and jointly estimate B1
function mesennls(y;
    κrange = LinRange(0.28, 1, 8),
    T1 = 1000,
    T2 = exp10.(LinRange(log10(15), log10(2000), 40)),
    TR = 1200,
    αex = π/2,
    αinv = π,
    TE = 10,
    nechoes = 32,
    β = 0
)

    if αinv isa AbstractArray
        κ = estimateb1.(y, [κrange], T1, [T2], TR, αex, [αinv], TE, nechoes)
    else
        κ = estimateb1.(y, [κrange], T1, [T2], TR, αex, αinv, TE, nechoes)
    end

    MWF = mesennls(y, κ, T1 = T1, T2 = T2, TR = TR, αex = αex, αinv = αinv,
                   TE = TE, nechoes = nechoes, β = β)

    return (MWF, κ)

end

# Estimate B1
function estimateb1(y, κrange, T1, T2, TR, αex, αinv, TE, nechoes)

    Nk = length(κrange)

    if αinv isa AbstractArray
        result = mesennls.([y], T1, [T2], κrange, TR, αex, [αinv], TE, nechoes)
    else
        result = mesennls.([y], T1, [T2], κrange, TR, αex, αinv, TE, nechoes)
    end
    err = [result[k][2] for k = 1:Nk]

    # Return the κ that minimizes the error
    upsamplefactor = 100
    errinterp = CubicSplineInterpolation(1:Nk, err)(1:1/upsamplefactor:Nk)
    κinterp = LinearInterpolation(1:Nk, κrange)(1:1/upsamplefactor:Nk)
    return κinterp[argmin(errinterp)]

end

# Estimate MWF for many voxels
# y is [N][nechoes]
function mesennls(y, κ;
    T1 = 1000,
    T2 = exp10.(LinRange(log10(15), log10(2000), 40)),
    TR = 1200,
    αex = π/2,
    αinv = π,
    TE = 10,
    nechoes = 32,
    β = 0
)

    N = length(y)

    if αinv isa AbstractArray
        result = mesennls.(y, T1, [T2], κ, TR, αex, [αinv], TE, nechoes, β) # [N]([D],[1])
    else
        result = mesennls.(y, T1, [T2], κ, TR, αex, αinv, TE, nechoes, β) # [N]([D],[1])
    end
    T2dist = [result[n][1] for n = 1:N] # [N][D]

    idx = T2 .<= 40
    MWF = [sum(T2dist[n][idx]) / sum(T2dist[n]) for n = 1:N]

    return MWF

end

# Estimate T2 distribution for one voxel
function mesennls(y, T1, T2, κ, TR, αex, αinv, TE, nechoes, β = 0)

    D = length(T2)

    # Generate data matrix A for unity M0
    A = createA(T1, T2, κ, TR, αex, αinv, TE, nechoes) # [nechoes,D]

    # Include regularization
    if β != 0
        A = [A; diagm(0 => fill(sqrt(β), D))]
        y = [y; zeros(eltype(y), D)]
    end

    # Run NNLS
    T2dist = nnls(A, y) # [D]

    # Compute the unregularized residual
    residual = A[1:nechoes,:] * T2dist .- y[1:nechoes] # [nechoes]

    return (T2dist, residual' * residual)

end

function createA(T1, T2, κ, TR, αex, αinv::AbstractArray{<:Real,1}, TE, nechoes)

    A = meseepg.(1, T1, T2, κ, TR, αex, [αinv], TE, nechoes) # [D][nechoes]
    A = reduce(hcat, A) # [nechoes,D]

end

function createA(T1, T2, κ, TR, αex, αinv::Real, TE, nechoes)

    A = meseepg.(1, T1, T2, κ, TR, αex, αinv, TE, nechoes) # [D][nechoes]
    A = reduce(hcat, A) # [nechoes,D]

end

# MESE via Extended Phase Graph
# Returns magnitude signal
# See Weigel et al. Extended Phase Graphs. JMRI 2015.
# See http://web.stanford.edu/class/rad229/Notes/1b-ExtendedPhaseGraphs.pdf
# See https://github.com/gopal-nataraj/mwf/blob/master/model/mese/mese.m
# See Prasloski et al. Stimulated Echo Correction in T2 Analysis. MRM 2012.
function meseepg(M0, T1, T2, κ, TR, αex, αinv::AbstractArray{<:Real,1}, TE, nechoes)

    length(αinv) == nechoes ||
        throw(DimensionMismatch("number of refocussing pulses = $(length(αinv)), nechoes = $nechoes"))

    # Set up rotation matrices
    Rex = rotateepg(π/2, κ * αex)
    Rinv = rotateepg.(0, κ * αinv)

    # Set up transition operator
    T! = x -> begin
        x[1,2:end] = x[1,1:end-1] # F_n transitions to F_{n+1}
        x[2,1:end-1] = x[2,2:end] # F_n* transitions to F_{n-1}*
        x[1,1] = conj(x[2,1]) # Use new F_1
    end

    # Set up relaxation/recovery
    E1 = exp(-TE/2T1)
    E2 = exp(-TE/2T2)
    E = Diagonal([E2, E2, E1])

    # Set up magnetization configuration
    M = zeros(ComplexF64, 3, 2nechoes + 1)

    # Get rough estimate of initial steady-state magnetization prior to excitation
    M[3,1] = M0 * (1 - exp(-TR / T1))

    # Set up output
    signal = zeros(nechoes)

    # Simulate MESE

    # Initial excitation
    M[:,1] = Rex * M[:,1]

    # Loop through each echo
    for e = 1:nechoes

        # Relaxation prior to refocussing pulse
        M = E * M
        M[3,1] += 1 - E1

        # Gradient prior to refocussing pulse
        T!(M)

        # Refocussing pulse
        M = Rinv[e] * M

        # Gradient following refocussing pulse
        T!(M)

        # Relaxation prior to echo time
        M = E * M
        M[3,1] += 1 - E1

        # Store echo amplitude
        signal[e] = abs(M[1,1])

    end

    return signal

end

function meseepg(M0, T1, T2, κ, TR, αex, αinv::Real, TE, nechoes)

    meseepg(M0, T1, T2, κ, TR, αex, fill(αinv, nechoes), TE, nechoes)

end

function rotateepg(θ, α)

    return [                  cos(α/2)^2      exp(2im*θ)*sin(α/2)^2 -im*exp( im*θ)*sin(α);
                  exp(-2im*θ)*sin(α/2)^2                 cos(α/2)^2  im*exp(-im*θ)*sin(α);
            -0.5im*exp(-im*θ)*sin(α)     0.5im*exp(im*θ)*sin(α)                    cos(α)]

end
