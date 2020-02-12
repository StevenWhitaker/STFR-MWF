function crlb(
    M0::Real,
    ff::Real,
    T1f::Real,
    T1s::Real,
    T2f::Real,
    T2s::Real,
    Δωf::Real,
    Δω::Real,
    κ::Real,
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1},
    computebiasgrad = true;
    Δωknown::Bool = true,
    κknown::Bool = true,
    σ::Real = 4e-3,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    N::Integer = 1000, # Number of noise realizations
    T::Integer = 20000,
    H::Integer = 4000,
    updateM0DistTrain::Bool = true,
    M0DistTrain = Uniform(eps(), 1.2),
    ffDistTrain = Uniform(0.03, 0.31),
    T1fDistTrain = Uniform(320, 480),
    T1sDistTrain = Uniform(800, 1200),
    T2fDistTrain = Uniform(16, 24),
    T2sDistTrain = Uniform(64, 96),
    ΔωfDistTrain = Uniform(0 * 2π, 35 * 2π),
    ΔωDistTrain = Uniform(-50 * 2π, 50 * 2π),
    κDistTrain = Uniform(0.8, 1.2),
    resetRNG::Bool = true
)

    if resetRNG
        Random.seed!(20191111)
    end

    I = inv(fisher([M0], [ff], [T1f], [T1s], [T2f], [T2s], [Δωf], [Δω], [κ],
        scanDesign, Δωknown = Δωknown, κknown = κknown, σ = σ)[1,:,:]) # [L,L]

    if computebiasgrad

        D = length(scanDesign[1][1])

        kwargs = (ρ = ρ, λ = λ,
            noisetest  = complex.(σ .* randn(D, N), σ .* randn(D, N)),
            noisetrain = complex.(σ .* randn(D, T), σ .* randn(D, T)),
            M0Train  = rand(M0DistTrain,  T),
            ffTrain  = rand(ffDistTrain,  T),
            T1fTrain = rand(T1fDistTrain, T),
            T1sTrain = rand(T1sDistTrain, T),
            T2fTrain = rand(T2fDistTrain, T),
            T2sTrain = rand(T2sDistTrain, T),
            ΔωfTrain = rand(ΔωfDistTrain, T),
            ΔωTrain  = rand(ΔωDistTrain,  T),
            κTrain   = rand(κDistTrain,   T),
            M0params = (updateM0DistTrain ? randparams(N = 10000, [1],
                ffDistTrain, T1fDistTrain, T1sDistTrain, T2fDistTrain,
                T2sDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain)[2:end] :
                nothing),
            f = randn(H, D + Δωknown + κknown),
            ph = rand(H)
        )

        grad = biasgrad(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, scanDesign;
            Δωknown = Δωknown, κknown = κknown, kwargs...)
        return grad' * I * grad

    else

        return I

    end

end

function fisher(
    M0::AbstractArray{<:Real,1}, # [N]
    ff::AbstractArray{<:Real,1}, # [N]
    T1f::AbstractArray{<:Real,1}, # [N]
    T1s::AbstractArray{<:Real,1}, # [N]
    T2f::AbstractArray{<:Real,1}, # [N]
    T2s::AbstractArray{<:Real,1}, # [N]
    Δωf::AbstractArray{<:Real,1}, # [N]
    Δω::AbstractArray{<:Real,1}, # [N]
    κ::AbstractArray{<:Real,1}, # [N]
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    Δωknown::Bool = true,
    κknown::Bool = true,
    σ::Real = 4e-3
)

    # Get STFR gradient w.r.t. unknown parameters
    (gradx,) = getstfrgradfunctions_v2(Δωknown = Δωknown, κknown = κknown)

    # Collect unknown and known parameters
    (x, ν) = collectdists(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Δωknown,
                                                                         κknown)

    # Make sqrt-covariance matrix
    Σ = [fill(σ, length(scanDesign[1][1]))]

    # Compute the Fisher information matrix
    return ScanDesign.fisher([gradx], x, ν, scanDesign, Σ) # [N,L,L]

end

function biasgrad(
    M0::Real,
    ff::Real,
    T1f::Real,
    T1s::Real,
    T2f::Real,
    T2s::Real,
    Δωf::Real,
    Δω::Real,
    κ::Real,
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    Δωknown::Bool = true,
    κknown::Bool = true,
    kwargs...
)

    # Collect unknown and known parameters
    (x, ν) = collectdists(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Δωknown,
                                                                         κknown)

    # If ν is empty, then make it an array, not an array of arrays
    isempty(ν) && (ν = Array{Float64,1}())

    # Create function to differentiate
    f = x -> stfrperkexpectation(x, ν, scanDesign; Δωknown = Δωknown,
                                 κknown = κknown, kwargs...)

    # Return the gradient
    grad = ForwardDiff.gradient(f, x)

    return grad

end

function stfrperkexpectation(
    xtest::AbstractArray{<:Real,1}, # [L]
    νtest::AbstractArray{<:Real,1}, # [K]
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    noisetest::AbstractArray{<:Complex}, # [D,N]
    Δωknown::Bool = true,
    κknown::Bool = true,
    noisetrain::AbstractArray{<:Complex}, # [D,T]
    M0params = nothing, # Output of randparams (excluding first entry), or nothing
    M0Train::AbstractArray{<:Real,1}, # [T]
    ffTrain::AbstractArray{<:Real,1}, # [T]
    T1fTrain::AbstractArray{<:Real,1}, # [T]
    T1sTrain::AbstractArray{<:Real,1}, # [T]
    T2fTrain::AbstractArray{<:Real,1}, # [T]
    T2sTrain::AbstractArray{<:Real,1}, # [T]
    ΔωfTrain::AbstractArray{<:Real,1}, # [T]
    ΔωTrain::AbstractArray{<:Real,1}, # [T]
    κTrain::AbstractArray{<:Real,1}, # [T]
    f::AbstractArray{<:Real,2}, # RFF frequency (unscaled) [H,D+K]
    ph::AbstractArray{<:Real,1} # RFF phase [H]
)

    # Set up the signal model
    signalModels = getstfrsignalModels_bias(scanDesign, Δωknown, κknown)

    # Generate noisy test data
    y = vcat([hcat(signalModels[i](xtest..., νtest...))
        for i = 1:length(signalModels)]...) # [D]
    y = abs.(y .+ noisetest) # [D,N]

    N = size(y, 2)

    # Run PERK
    (xhat,) = stfrperk_bias(y, scanDesign, (Δωknown ? fill(νtest[1], N) : nothing),
        (κknown ? fill(νtest[end], N) : nothing), noisetrain, M0params, M0Train,
        ffTrain, T1fTrain, T1sTrain, T2fTrain, T2sTrain, ΔωfTrain, ΔωTrain,
        κTrain, f, ph, ρ = ρ, λ = λ)

    # Compute expected value of MWF
    return mean(xhat[2,:])

end

function stfrperkexpectation(
    M0::Real,
    ff::Real,
    T1f::Real,
    T1s::Real,
    T2f::Real,
    T2s::Real,
    Δωf::Real,
    Δω::Real,
    κ::Real,
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    Δωknown::Bool = true,
    κknown::Bool = true,
    σ::Real = 4e-3,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    N::Integer = 1000, # Number of noise realizations
    T::Integer = 20000,
    H::Integer = 4000,
    updateM0DistTrain::Bool = true,
    M0DistTrain = Uniform(eps(), 1.2),
    ffDistTrain = Uniform(0.03, 0.31),
    T1fDistTrain = Uniform(320, 480),
    T1sDistTrain = Uniform(800, 1200),
    T2fDistTrain = Uniform(16, 24),
    T2sDistTrain = Uniform(64, 96),
    ΔωfDistTrain = Uniform(0 * 2π, 35 * 2π),
    ΔωDistTrain = Uniform(-50 * 2π, 50 * 2π),
    κDistTrain = Uniform(0.8, 1.2),
    resetRNG::Bool = true
)

    if resetRNG
        Random.seed!(20191121)
    end

    # Collect unknown and known parameters
    (x, ν) = collectdists(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Δωknown,
                                                                         κknown)

    # If ν is empty, then make it an array, not an array of arrays
    isempty(ν) && (ν = Array{Float64,1}())

    # Create kwargs
    D = length(scanDesign[1][1])

    kwargs = (ρ = ρ, λ = λ,
        noisetest  = complex.(σ .* randn(D, N), σ .* randn(D, N)),
        noisetrain = complex.(σ .* randn(D, T), σ .* randn(D, T)),
        M0Train  = rand(M0DistTrain,  T),
        ffTrain  = rand(ffDistTrain,  T),
        T1fTrain = rand(T1fDistTrain, T),
        T1sTrain = rand(T1sDistTrain, T),
        T2fTrain = rand(T2fDistTrain, T),
        T2sTrain = rand(T2sDistTrain, T),
        ΔωfTrain = rand(ΔωfDistTrain, T),
        ΔωTrain  = rand(ΔωDistTrain,  T),
        κTrain   = rand(κDistTrain,   T),
        M0params = (updateM0DistTrain ? randparams(N = 10000, [1],
            ffDistTrain, T1fDistTrain, T1sDistTrain, T2fDistTrain,
            T2sDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain)[2:end] :
            nothing),
        f = randn(H, D + Δωknown + κknown),
        ph = rand(H)
    )

    # Return expected value of MWF
    return stfrperkexpectation(x, ν, scanDesign; Δωknown = Δωknown,
                               κknown = κknown, kwargs...)

end

function stfrblochsimperkexpectation(
    xtest::AbstractArray{<:Real,1}, # [L]
    νtest::AbstractArray{<:Real,1}, # [K]
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    noisetest::AbstractArray{<:Complex}, # [D,N]
    Δωknown::Bool = true,
    κknown::Bool = true,
    noisetrain::AbstractArray{<:Complex}, # [D,T]
    M0params = nothing, # Output of randparams (excluding first entry), or nothing
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
    ph::AbstractArray{<:Real,1} # RFF phase [H]
)

    # Set up the signal model
    signalModels = getstfrblochsimsignalModels_bias(scanDesign, Δωknown, κknown)

    # Generate noisy test data
    y = vcat([hcat(signalModels[i](xtest..., νtest...))
        for i = 1:length(signalModels)]...) # [D]
    y = abs.(y .+ noisetest) # [D,N]

    N = size(y, 2)

    # Run PERK
    (xhat,) = stfrblochsimperk_bias(y, scanDesign,
        (Δωknown ? fill(νtest[1], N) : nothing),
        (κknown ? fill(νtest[end], N) : nothing), noisetrain, M0params, M0Train,
        ffTrain, fmTrain, T1fTrain, T1sTrain, T1mTrain, T2fTrain, T2sTrain,
        T2mTrain, τfsTrain, τfmTrain, ΔωfTrain, ΔωTrain, κTrain, f, ph, ρ = ρ,
        λ = λ)

    # Compute expected value of MWF
    return mean(xhat[2,:])

end

function stfrblochsimperkexpectation(
    M0::Real,
    ff::Real,
    fm::Real,
    T1f::Real,
    T1s::Real,
    T1m::Real,
    T2f::Real,
    T2s::Real,
    T2m::Real,
    τfs::Real,
    τfm::Real,
    Δωf::Real,
    Δω::Real,
    κ::Real,
    scanDesign::AbstractArray{<:AbstractArray{<:AbstractArray{<:Real,1},1},1};
    Δωknown::Bool = true,
    κknown::Bool = true,
    σ::Real = 4e-3,
    ρ::Real = 2.0^-60,
    λ::Real = 2.0^3.5,
    N::Integer = 1000, # Number of noise realizations
    T::Integer = 20000,
    H::Integer = 4000,
    updateM0DistTrain::Bool = true,
    M0DistTrain = Uniform(eps(), 1.2),
    ffDistTrain = Uniform(0.03, 0.31),
    fmDistTrain = Uniform(0.03, 0.31),
    T1fDistTrain = Uniform(320, 480),
    T1sDistTrain = Uniform(800, 1200),
    T1mDistTrain = Uniform(800, 3000),
    T2fDistTrain = Uniform(16, 24),
    T2sDistTrain = Uniform(64, 96),
    T2mDistTrain = Uniform(0.01, 0.1),
    τfsDistTrain = Uniform(80, 150),
    τfmDistTrain = Uniform(40, 75),
    ΔωfDistTrain = Uniform(0 * 2π, 35 * 2π),
    ΔωDistTrain = Uniform(-50 * 2π, 50 * 2π),
    κDistTrain = Uniform(0.8, 1.2),
    resetRNG::Bool = true
)

    if resetRNG
        Random.seed!(20191204)
    end

    # Collect unknown and known parameters
    (x, ν) = collectdists(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm,
                          Δωf, Δω, κ, Δωknown, κknown)

    # If ν is empty, then make it an array, not an array of arrays
    isempty(ν) && (ν = Array{Float64,1}())

    # Create kwargs
    D = length(scanDesign[1][1])

    kwargs = (ρ = ρ, λ = λ,
        noisetest  = complex.(σ .* randn(D, N), σ .* randn(D, N)),
        noisetrain = complex.(σ .* randn(D, T), σ .* randn(D, T)),
        M0Train  = rand(M0DistTrain,  T),
        ffTrain  = rand(ffDistTrain,  T),
        fmTrain  = rand(fmDistTrain,  T),
        T1fTrain = rand(T1fDistTrain, T),
        T1sTrain = rand(T1sDistTrain, T),
        T1mTrain = rand(T1mDistTrain, T),
        T2fTrain = rand(T2fDistTrain, T),
        T2sTrain = rand(T2sDistTrain, T),
        T2mTrain = rand(T2mDistTrain, T),
        τfsTrain = rand(τfsDistTrain, T),
        τfmTrain = rand(τfmDistTrain, T),
        ΔωfTrain = rand(ΔωfDistTrain, T),
        ΔωTrain  = rand(ΔωDistTrain,  T),
        κTrain   = rand(κDistTrain,   T),
        M0params = (updateM0DistTrain ? randparams(N = 10000, [1],
            ffDistTrain, fmDistTrain, T1fDistTrain, T1sDistTrain, T1mDistTrain,
            T2fDistTrain, T2sDistTrain, T2mDistTrain, τfsDistTrain,
            τfmDistTrain, ΔωfDistTrain, ΔωDistTrain, κDistTrain)[2:end] :
            nothing),
        f = randn(H, D + Δωknown + κknown),
        ph = rand(H)
    )

    # Return expected value of MWF
    return stfrblochsimperkexpectation(x, ν, scanDesign; Δωknown = Δωknown,
                                       κknown = κknown, kwargs...)

end

"""
    getstfrgradfunctions_v2(; M0known, ffknown, T1fknown, T1sknown, T2fknown,
        T2sknown, Δωfknown, Δωknown, κknown)

Generate functions to be used to evaluate the gradient of the two-compartment
STFR signal model without exchange.

# Arguments
- `M0known::Bool = false`: Whether M0 is known
- `ffknown::Bool = false`: Whether ff is known
- `T1fknown::Bool = false`: Whether T1f is known
- `T1sknown::Bool = false`: Whether T1s is known
- `T2fknown::Bool = false`: Whether T2f is known
- `T2sknown::Bool = false`: Whether T2s is known
- `Δωfknown::Bool = false`: Whether Δωf is known
- `Δωknown::Bool = true`: Whether Δω is known
- `κknown::Bool = true`: Whether κ is known

# Return
- `gradx::Function`: Takes latent, known, and scan parameters and outputs the
  partial derivatives of the STFR signal model with respect to the latent
  parameters
- `gradxp::Function`: Takes latent, known, and scan parameters and outputs the
  mixed partial derivatives of the STFR signal model with respect to the latent
  parameters and the scan parameters
"""
function getstfrgradfunctions_v2(;
    M0known::Bool = false,
    ffknown::Bool = false,
    T1fknown::Bool = false,
    T1sknown::Bool = false,
    T2fknown::Bool = false,
    T2sknown::Bool = false,
    Δωfknown::Bool = false,
    Δωknown::Bool = true,
    κknown::Bool = true,
)

    # Set up symbolic variables
    # Latent/known parameters
    M0  = symbols("M0",    positive    = true)
    ff  = symbols("ff",    nonnegative = true)
    T1f = symbols("T1f",   positive    = true)
    T1s = symbols("T1s",   positive    = true)
    T2f = symbols("T2f",   positive    = true)
    T2s = symbols("T2s",   positive    = true)
    Δωf = symbols("dwf",   real        = true)
    Δω  = symbols("dw",    real        = true)
    κ   = symbols("kappa", positive    = true)
    # Scan parameters
    Tfree = symbols("Tfree", positive    = true)
    Tg    = symbols("Tg",    nonnegative = true)
    α     = symbols("alpha", positive    = true)
    β     = symbols("beta",  nonnegative = true)
    ϕ     = symbols("phi",   real        = true)
    TE    = symbols("TE",    nonnegative = true)

    # Make list of latent and known parameters
    latent = Array{SymPy.Sym,1}()
    known  = Array{SymPy.Sym,1}()
    if M0known  push!(known, M0)  else push!(latent, M0)  end
    if ffknown  push!(known, ff)  else push!(latent, ff)  end
    if T1fknown push!(known, T1f) else push!(latent, T1f) end
    if T1sknown push!(known, T1s) else push!(latent, T1s) end
    if T2fknown push!(known, T2f) else push!(latent, T2f) end
    if T2sknown push!(known, T2s) else push!(latent, T2s) end
    if Δωfknown push!(known, Δωf) else push!(latent, Δωf) end
    if Δωknown  push!(known, Δω)  else push!(latent, Δω)  end
    if κknown   push!(known, κ)   else push!(latent, κ)   end

    # Collect scan parameters
    scan = [Tfree, Tg, α, β, ϕ, TE]

    # Collect all parameters
    params = [latent; known; scan]

    # Create the signal equation
    n = exp(-Tg/T1f) * (1 - exp(-Tfree/T1f)) * cos(κ * β) + (1 - exp(-Tg/T1f))
    d = 1 - exp(-Tg/T1f - Tfree/T2f) * sin(κ * α) * sin(κ * β) *
      cos((Δωf + Δω) * Tfree/1000 - ϕ) - exp(-Tg/T1f - Tfree/T1f) *
      cos(κ * α) * cos(κ * β)

    Mf = M0 * sin(κ * α) * n / d
    Ms = subs(Mf, T1f => T1s, T2f => T2s, Δωf => 0)
    M = ff * Mf + (1 - ff) * Ms

    # Calculate the gradients
    gradx_sym  = [diff(M, x) for x in latent]
    gradxp_sym = [diff(g, x) for g in gradx_sym, x in scan]

    # Convert to Julia functions
    gradx_arr  = [lambdify(g, params) for g in gradx_sym]
    gradxp_arr = [lambdify(g, params) for g in gradxp_sym]
    gradx  = (args...) -> [gradx_arr[i].(args...)    for i = 1:length(latent)]
    gradxp = (args...) -> [gradxp_arr[i,j].(args...) for i = 1:length(latent),
                                                         j = 1:length(scan)]

    return (gradx, gradxp)

end
