function testdata_noise(;
    N = 1000,
    σ::Real = 2e-3
)

    noise = complex.(rand(Normal(0, σ), N), rand(Normal(0, σ), N))
    return abs.(noise)

end

function testdata_stfr(;
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1} = ones(N),
    ff::AbstractArray{<:Real,1} = rand(Uniform(0.03, 0.31), N),
    T1f::AbstractArray{<:Real,1} = rand(Uniform(320, 480), N),
    T1s::AbstractArray{<:Real,1} = rand(Uniform(800, 1200), N),
    T2f::AbstractArray{<:Real,1} = rand(Uniform(16, 24), N),
    T2s::AbstractArray{<:Real,1} = rand(Uniform(64, 96), N),
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0 * 2π, 35 * 2π), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3
)

    (Tfree, Tg, α, β, ϕ, TE) = load(scanDesignFile, "P")[1] # [D]
    y = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δω', κ', Tfree, Tg, α, β, ϕ, TE) # [D,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # # Inspect SNR
    # noise = testdata_noise(N = 1000, σ = σ)
    # SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    # @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ]

    return (y, x)

end

# Simulate test data using 3 compartment Bloch simulation
# f is for fast-relaxing compartment (myelin water)
# s is for slow-relaxing compartment (non-myelin water)
# m is for macromolecular pool
# fm, T1m, and T2m were inspired by http://mriquestions.com/how-to-predict-t1--t2.html
# τfs was selected to reflect the general impression I got when looking at the literature
# τfm was selected based on 12.5 / kBA ≈ 2 in Figure 4 in "Characterizing WM with MT and T2"
#     by Stanisz et al.
function testdata_stfrblochsim(;
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1} = ones(N),
    ff::AbstractArray{<:Real,1} = rand(Uniform(0.03, 0.31), N),
    fm::AbstractArray{<:Real,1} = rand(Uniform(0.03, 0.31), N), # Just use the same as ff
    T1f::AbstractArray{<:Real,1} = rand(Uniform(320, 480), N),
    T1s::AbstractArray{<:Real,1} = rand(Uniform(800, 1200), N),
    T1m::AbstractArray{<:Real,1} = rand(Uniform(800, 3000), N),
    T2f::AbstractArray{<:Real,1} = rand(Uniform(16, 24), N),
    T2s::AbstractArray{<:Real,1} = rand(Uniform(64, 96), N),
    T2m::AbstractArray{<:Real,1} = rand(Uniform(0.01, 0.1), N),
    τfs::AbstractArray{<:Real,1} = rand(Uniform(80, 150), N), # Ballpark value
    τfm::AbstractArray{<:Real,1} = rand(Uniform(40, 75), N), # About twice as fast as τfs
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0 * 2π, 35 * 2π), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3
)

    (Tfree, Tg, α, β, ϕ, TE) = load(scanDesignFile, "P")[1] # [D]
    fs = 1 .- ff .- fm # [N]
    τsf = τfs .* fs ./ ff # [N] Assume chemical equilibrium
    τsm = fill(Inf, N) # [N] No exchange
    τmf = fill(Inf, N) # [N] No exchange
    τms = fill(Inf, N) # [N] No exchange
    frac = permutedims([[ff[i], fs[i], fm[i]] for i = 1:N])
    T1 = permutedims([[T1f[i], T1s[i], T1m[i]] for i = 1:N])
    T2 = permutedims([[T2f[i], T2s[i], T2m[i]] for i = 1:N])
    Δωtmp = permutedims([[Δω[i] + Δωf[i], Δω[i], Δω[i]] for i = 1:N])
    τ = permutedims([[τfs[i], τfm[i], τsf[i], τsm[i], τmf[i], τms[i]] for i = 1:N])
    y = stfrblochsim.(M0', frac, T1, T2, Δωtmp, τ, κ', Tfree, Tg, α, β, ϕ, TE) # [D,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ]

    return (y, x)

end

# Simulate test data using 3 compartment Bloch simulation
# f is for fast-relaxing compartment (myelin water)
# s is for slow-relaxing compartment (non-myelin water)
# m is for macromolecular pool
# fm, T1m, and T2m were inspired by http://mriquestions.com/how-to-predict-t1--t2.html
# τfs was selected to reflect the general impression I got when looking at the literature
# τfm was selected based on 12.5 / kBA ≈ 2 in Figure 4 in "Characterizing WM with MT and T2"
#     by Stanisz et al.
function testdata_meseblochsim(;
    TR::Real = 1200,
    αex::Real = π/2,
    αinv::Real = π,
    TE::Real = 10,
    nechoes::Integer = 32,
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1} = rand(Delta(1), N),
    ff::AbstractArray{<:Real,1} = rand(Uniform(0.03, 0.31), N),
    fm::AbstractArray{<:Real,1} = rand(Uniform(0.03, 0.31), N), # Just use the same as ff
    T1f::AbstractArray{<:Real,1} = rand(Uniform(320, 480), N),
    T1s::AbstractArray{<:Real,1} = rand(Uniform(800, 1200), N),
    T1m::AbstractArray{<:Real,1} = rand(Uniform(800, 3000), N),
    T2f::AbstractArray{<:Real,1} = rand(Uniform(16, 24), N),
    T2s::AbstractArray{<:Real,1} = rand(Uniform(64, 96), N),
    T2m::AbstractArray{<:Real,1} = rand(Uniform(0.01, 0.1), N),
    τfs::AbstractArray{<:Real,1} = rand(Uniform(80, 150), N), # Ballpark value
    τfm::AbstractArray{<:Real,1} = rand(Uniform(40, 75), N), # About twice as fast as τfs
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0, 35), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50, 50), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3,
    nspins::Integer = 30, # Set to 1 to use ideal spoiling
    ncycles::Real = 1
)

    fs = 1 .- ff .- fm # [N]
    τsf = τfs .* fs ./ ff # [N] Assume chemical equilibrium
    τsm = fill(Inf, N) # [N] No exchange
    τmf = fill(Inf, N) # [N] No exchange
    τms = fill(Inf, N) # [N] No exchange
    frac = [[ff[i], fs[i], fm[i]] for i = 1:N]
    T1 = [[T1f[i], T1s[i], T1m[i]] for i = 1:N]
    T2 = [[T2f[i], T2s[i], T2m[i]] for i = 1:N]
    Δωtmp = [[Δω[i] + Δωf[i], Δω[i], Δω[i]] for i = 1:N]
    τ = [[τfs[i], τfm[i], τsf[i], τsm[i], τmf[i], τms[i]] for i = 1:N]
    if nspins == 1
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes) # [N][nechoes]
    else
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes,
                          nspins, ncycles) # [N][nechoes]
    end
    y = reduce(hcat, y) # [nechoes,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf, Δω, κ]

    return (y, x)

end

function testdata_brainweb(scan::Symbol, tissues = (:wm, :gm);
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    Δω::AbstractArray{<:Real,2} = (genmap(trues(181, 217), 50, 0, 100) .- 50) .* 2π,
    κ::AbstractArray{<:Real,2} = genmap(trues(181, 217), 1, 0.8, 1.2),
    σ::Real = 2e-3,
    meseScanDesign = (TR = 1200, αex = π/2, αinv = π, TE = 10, nechoes = 32),
    useexchange::Bool = true,
    D::Int = scan == :mese ? 32 : 11 # 11 STFR scans, 32 MESE scans
)

    # Read in the brain
    obj = readbrainweb(useexchange ? 3 : 300)

    # Preallocate the output
    data = zeros(size(obj)..., D)

    # Select the tissues of interest
    (mask, M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf) = getparams(obj, tissues...)
    N = count(mask)

    # Simulate the scan
    if scan == :stfr

        (y,) = testdata_stfr(scanDesignFile = scanDesignFile, N = N,
            M0 = M0, ff = ff, T1f = T1f, T1s = T1s, T2f = T2f, T2s = T2s,
            Δωf = Δωf, Δω = Δω[mask], κ = κ[mask], σ = σ)

    elseif scan == :stfrblochsim

        (y,) = testdata_stfrblochsim(scanDesignFile = scanDesignFile, N = N,
            M0 = M0, ff = ff, fm = fm, T1f = T1f, T1s = T1s, T1m = T1m,
            T2f = T2f, T2s = T2s, T2m = T2m, τfs = τfs, τfm = τfm,
            Δωf = Δωf, Δω = Δω[mask], κ = κ[mask], σ = σ)

    elseif scan == :mese

        (y,) = testdata_meseblochsim(N = N,
            M0 = M0, ff = ff, fm = fm, T1f = T1f, T1s = T1s, T1m = T1m,
            T2f = T2f, T2s = T2s, T2m = T2m, τfs = τfs, τfm = τfm,
            Δωf = Δωf, Δω = Δω[mask], κ = κ[mask], σ = σ; meseScanDesign...)

    else

        throw(ArgumentError("invalid scan"))

    end

    # Reshape data
    for d = 1:D

        @view(data[:,:,d])[mask] = y[d,:]

    end

    return (data, obj, mask)

end

function testdata_stfrblochsim_9comp(;
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1},
    ff1::AbstractArray{<:Real,1},
    ff2::AbstractArray{<:Real,1},
    ff3::AbstractArray{<:Real,1},
    fs1::AbstractArray{<:Real,1},
    fs2::AbstractArray{<:Real,1},
    fs3::AbstractArray{<:Real,1},
    fm1::AbstractArray{<:Real,1},
    fm2::AbstractArray{<:Real,1},
    fm3::AbstractArray{<:Real,1},
    T1f1::AbstractArray{<:Real,1},
    T1f2::AbstractArray{<:Real,1},
    T1f3::AbstractArray{<:Real,1},
    T1s1::AbstractArray{<:Real,1},
    T1s2::AbstractArray{<:Real,1},
    T1s3::AbstractArray{<:Real,1},
    T1m1::AbstractArray{<:Real,1},
    T1m2::AbstractArray{<:Real,1},
    T1m3::AbstractArray{<:Real,1},
    T2f1::AbstractArray{<:Real,1},
    T2f2::AbstractArray{<:Real,1},
    T2f3::AbstractArray{<:Real,1},
    T2s1::AbstractArray{<:Real,1},
    T2s2::AbstractArray{<:Real,1},
    T2s3::AbstractArray{<:Real,1},
    T2m1::AbstractArray{<:Real,1},
    T2m2::AbstractArray{<:Real,1},
    T2m3::AbstractArray{<:Real,1},
    τf1f2::AbstractArray{<:Real,1},
    τf1f3::AbstractArray{<:Real,1},
    τf1s1::AbstractArray{<:Real,1},
    τf1s2::AbstractArray{<:Real,1},
    τf1s3::AbstractArray{<:Real,1},
    τf1m1::AbstractArray{<:Real,1},
    τf1m2::AbstractArray{<:Real,1},
    τf1m3::AbstractArray{<:Real,1},
    τf2f3::AbstractArray{<:Real,1},
    τf2s1::AbstractArray{<:Real,1},
    τf2s2::AbstractArray{<:Real,1},
    τf2s3::AbstractArray{<:Real,1},
    τf2m1::AbstractArray{<:Real,1},
    τf2m2::AbstractArray{<:Real,1},
    τf2m3::AbstractArray{<:Real,1},
    τf3s1::AbstractArray{<:Real,1},
    τf3s2::AbstractArray{<:Real,1},
    τf3s3::AbstractArray{<:Real,1},
    τf3m1::AbstractArray{<:Real,1},
    τf3m2::AbstractArray{<:Real,1},
    τf3m3::AbstractArray{<:Real,1},
    τs1s2::AbstractArray{<:Real,1},
    τs1s3::AbstractArray{<:Real,1},
    τs2s3::AbstractArray{<:Real,1},
    τm1m2::AbstractArray{<:Real,1},
    τm1m3::AbstractArray{<:Real,1},
    τm2m3::AbstractArray{<:Real,1},
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0, 35), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3
)

    (Tfree, Tg, α, β, ϕ, TE) = load(scanDesignFile, "P")[1] # [D]
    if !all(+(ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3) .≈ 1)
        throw(ArgumentError("fractions add to values in the range $(extrema(+(
            ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3))), not 1"))
    end
    τf2f1 = τf1f2 .* ff2 ./ ff1 # [N] Assume chemical equilibrium
    τf3f1 = τf1f3 .* ff3 ./ ff1 # [N] Assume chemical equilibrium
    τf3f2 = τf2f3 .* ff3 ./ ff2 # [N] Assume chemical equilibrium
    τs2s1 = τs1s2 .* fs2 ./ fs1 # [N] Assume chemical equilibrium
    τs3s1 = τs1s3 .* fs3 ./ fs1 # [N] Assume chemical equilibrium
    τs3s2 = τs2s3 .* fs3 ./ fs2 # [N] Assume chemical equilibrium
    τm2m1 = τm1m2 .* fm2 ./ fm1 # [N] Assume chemical equilibrium
    τm3m1 = τm1m3 .* fm3 ./ fm1 # [N] Assume chemical equilibrium
    τm3m2 = τm2m3 .* fm3 ./ fm2 # [N] Assume chemical equilibrium
    τs1f1 = τf1s1 .* fs1 ./ ff1 # [N] Assume chemical equilibrium
    τs1f2 = τf2s1 .* fs1 ./ ff2 # [N] Assume chemical equilibrium
    τs1f3 = τf3s1 .* fs1 ./ ff3 # [N] Assume chemical equilibrium
    τs2f1 = τf1s2 .* fs2 ./ ff1 # [N] Assume chemical equilibrium
    τs2f2 = τf2s2 .* fs2 ./ ff2 # [N] Assume chemical equilibrium
    τs2f3 = τf3s2 .* fs2 ./ ff3 # [N] Assume chemical equilibrium
    τs3f1 = τf1s3 .* fs3 ./ ff1 # [N] Assume chemical equilibrium
    τs3f2 = τf2s3 .* fs3 ./ ff2 # [N] Assume chemical equilibrium
    τs3f3 = τf3s3 .* fs3 ./ ff3 # [N] Assume chemical equilibrium
    τs1m1 = fill(Inf, N) # [N] No exchange
    τs1m2 = fill(Inf, N) # [N] No exchange
    τs1m3 = fill(Inf, N) # [N] No exchange
    τs2m1 = fill(Inf, N) # [N] No exchange
    τs2m2 = fill(Inf, N) # [N] No exchange
    τs2m3 = fill(Inf, N) # [N] No exchange
    τs3m1 = fill(Inf, N) # [N] No exchange
    τs3m2 = fill(Inf, N) # [N] No exchange
    τs3m3 = fill(Inf, N) # [N] No exchange
    τm1f1 = fill(Inf, N) # [N] No exchange
    τm1f2 = fill(Inf, N) # [N] No exchange
    τm1f3 = fill(Inf, N) # [N] No exchange
    τm2f1 = fill(Inf, N) # [N] No exchange
    τm2f2 = fill(Inf, N) # [N] No exchange
    τm2f3 = fill(Inf, N) # [N] No exchange
    τm3f1 = fill(Inf, N) # [N] No exchange
    τm3f2 = fill(Inf, N) # [N] No exchange
    τm3f3 = fill(Inf, N) # [N] No exchange
    τm1s1 = fill(Inf, N) # [N] No exchange
    τm1s2 = fill(Inf, N) # [N] No exchange
    τm1s3 = fill(Inf, N) # [N] No exchange
    τm2s1 = fill(Inf, N) # [N] No exchange
    τm2s2 = fill(Inf, N) # [N] No exchange
    τm2s3 = fill(Inf, N) # [N] No exchange
    τm3s1 = fill(Inf, N) # [N] No exchange
    τm3s2 = fill(Inf, N) # [N] No exchange
    τm3s3 = fill(Inf, N) # [N] No exchange
    frac = permutedims([[ff1[i], ff2[i], ff3[i], fs1[i], fs2[i], fs3[i],
                         fm1[i], fm2[i], fm3[i]] for i = 1:N])
    T1 = permutedims([[T1f1[i], T1f2[i], T1f3[i], T1s1[i], T1s2[i], T1s3[i],
                       T1m1[i], T1m2[i], T1m3[i]] for i = 1:N])
    T2 = permutedims([[T2f1[i], T2f2[i], T2f3[i], T2s1[i], T2s2[i], T2s3[i],
                       T2m1[i], T2m2[i], T2m3[i]] for i = 1:N])
    Δωtmp = permutedims([[Δω[i] + Δωf[i], Δω[i] + Δωf[i], Δω[i] + Δωf[i],
                          Δω[i], Δω[i], Δω[i], Δω[i], Δω[i], Δω[i]] for i = 1:N])
    τ = permutedims([[τf1f2[i], τf1f3[i], τf1s1[i], τf1s2[i], τf1s3[i], τf1m1[i], τf1m2[i], τf1m3[i],
                      τf2f1[i], τf2f3[i], τf2s1[i], τf2s2[i], τf2s3[i], τf2m1[i], τf2m2[i], τf2m3[i],
                      τf3f1[i], τf3f2[i], τf3s1[i], τf3s2[i], τf3s3[i], τf3m1[i], τf3m2[i], τf3m3[i],
                      τs1f1[i], τs1f2[i], τs1f3[i], τs1s2[i], τs1s3[i], τs1m1[i], τs1m2[i], τs1m3[i],
                      τs2f1[i], τs2f2[i], τs2f3[i], τs2s1[i], τs2s3[i], τs2m1[i], τs2m2[i], τs2m3[i],
                      τs3f1[i], τs3f2[i], τs3f3[i], τs3s1[i], τs3s2[i], τs3m1[i], τs3m2[i], τs3m3[i],
                      τm1f1[i], τm1f2[i], τm1f3[i], τm1s1[i], τm1s2[i], τm1s3[i], τm1m2[i], τm1m3[i],
                      τm2f1[i], τm2f2[i], τm2f3[i], τm2s1[i], τm2s2[i], τm2s3[i], τm2m1[i], τm2m3[i],
                      τm3f1[i], τm3f2[i], τm3f3[i], τm3s1[i], τm3s2[i], τm3s3[i], τm3m1[i], τm3m2[i]] for i = 1:N])
    y = stfrblochsim.(M0', frac, T1, T2, Δωtmp, τ, κ', Tfree, Tg, α, β, ϕ, TE) # [D,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3, T1f1, T1f2, T1f3,
         T1s1, T1s2, T1s3, T1m1, T1m2, T1m3, T2f1, T2f2, T2f3, T2s1, T2s2, T2s3,
         T2m1, T2m2, T2m3, τf1f2, τf1f3, τf1s1, τf1s2, τf1s3, τf1m1, τf1m2,
         τf1m3, τf2f3, τf2s1, τf2s2, τf2s3, τf2m1, τf2m2, τf2m3, τf3s1, τf3s2,
         τf3s3, τf3m1, τf3m2, τf3m3, τs1s2, τs1s3, τs2s3, τm1m2, τm1m3, τm2m3,
         Δωf, Δω, κ]

    return (y, x)

end

function testdata_meseblochsim_9comp(;
    TR::Real = 1200,
    αex::Real = π/2,
    αinv::Real = π,
    TE::Real = 10,
    nechoes::Integer = 32,
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1},
    ff1::AbstractArray{<:Real,1},
    ff2::AbstractArray{<:Real,1},
    ff3::AbstractArray{<:Real,1},
    fs1::AbstractArray{<:Real,1},
    fs2::AbstractArray{<:Real,1},
    fs3::AbstractArray{<:Real,1},
    fm1::AbstractArray{<:Real,1},
    fm2::AbstractArray{<:Real,1},
    fm3::AbstractArray{<:Real,1},
    T1f1::AbstractArray{<:Real,1},
    T1f2::AbstractArray{<:Real,1},
    T1f3::AbstractArray{<:Real,1},
    T1s1::AbstractArray{<:Real,1},
    T1s2::AbstractArray{<:Real,1},
    T1s3::AbstractArray{<:Real,1},
    T1m1::AbstractArray{<:Real,1},
    T1m2::AbstractArray{<:Real,1},
    T1m3::AbstractArray{<:Real,1},
    T2f1::AbstractArray{<:Real,1},
    T2f2::AbstractArray{<:Real,1},
    T2f3::AbstractArray{<:Real,1},
    T2s1::AbstractArray{<:Real,1},
    T2s2::AbstractArray{<:Real,1},
    T2s3::AbstractArray{<:Real,1},
    T2m1::AbstractArray{<:Real,1},
    T2m2::AbstractArray{<:Real,1},
    T2m3::AbstractArray{<:Real,1},
    τf1f2::AbstractArray{<:Real,1},
    τf1f3::AbstractArray{<:Real,1},
    τf1s1::AbstractArray{<:Real,1},
    τf1s2::AbstractArray{<:Real,1},
    τf1s3::AbstractArray{<:Real,1},
    τf1m1::AbstractArray{<:Real,1},
    τf1m2::AbstractArray{<:Real,1},
    τf1m3::AbstractArray{<:Real,1},
    τf2f3::AbstractArray{<:Real,1},
    τf2s1::AbstractArray{<:Real,1},
    τf2s2::AbstractArray{<:Real,1},
    τf2s3::AbstractArray{<:Real,1},
    τf2m1::AbstractArray{<:Real,1},
    τf2m2::AbstractArray{<:Real,1},
    τf2m3::AbstractArray{<:Real,1},
    τf3s1::AbstractArray{<:Real,1},
    τf3s2::AbstractArray{<:Real,1},
    τf3s3::AbstractArray{<:Real,1},
    τf3m1::AbstractArray{<:Real,1},
    τf3m2::AbstractArray{<:Real,1},
    τf3m3::AbstractArray{<:Real,1},
    τs1s2::AbstractArray{<:Real,1},
    τs1s3::AbstractArray{<:Real,1},
    τs2s3::AbstractArray{<:Real,1},
    τm1m2::AbstractArray{<:Real,1},
    τm1m3::AbstractArray{<:Real,1},
    τm2m3::AbstractArray{<:Real,1},
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0, 35), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3,
    nspins::Integer = 30, # Set to 1 to use ideal spoiling
    ncycles::Real = 1
)

    if !all(+(ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3) .≈ 1)
        throw(ArgumentError("fractions add to values in the range $(extrema(+(
            ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3))), not 1"))
    end
    τf2f1 = τf1f2 .* ff2 ./ ff1 # [N] Assume chemical equilibrium
    τf3f1 = τf1f3 .* ff3 ./ ff1 # [N] Assume chemical equilibrium
    τf3f2 = τf2f3 .* ff3 ./ ff2 # [N] Assume chemical equilibrium
    τs2s1 = τs1s2 .* fs2 ./ fs1 # [N] Assume chemical equilibrium
    τs3s1 = τs1s3 .* fs3 ./ fs1 # [N] Assume chemical equilibrium
    τs3s2 = τs2s3 .* fs3 ./ fs2 # [N] Assume chemical equilibrium
    τm2m1 = τm1m2 .* fm2 ./ fm1 # [N] Assume chemical equilibrium
    τm3m1 = τm1m3 .* fm3 ./ fm1 # [N] Assume chemical equilibrium
    τm3m2 = τm2m3 .* fm3 ./ fm2 # [N] Assume chemical equilibrium
    τs1f1 = τf1s1 .* fs1 ./ ff1 # [N] Assume chemical equilibrium
    τs1f2 = τf2s1 .* fs1 ./ ff2 # [N] Assume chemical equilibrium
    τs1f3 = τf3s1 .* fs1 ./ ff3 # [N] Assume chemical equilibrium
    τs2f1 = τf1s2 .* fs2 ./ ff1 # [N] Assume chemical equilibrium
    τs2f2 = τf2s2 .* fs2 ./ ff2 # [N] Assume chemical equilibrium
    τs2f3 = τf3s2 .* fs2 ./ ff3 # [N] Assume chemical equilibrium
    τs3f1 = τf1s3 .* fs3 ./ ff1 # [N] Assume chemical equilibrium
    τs3f2 = τf2s3 .* fs3 ./ ff2 # [N] Assume chemical equilibrium
    τs3f3 = τf3s3 .* fs3 ./ ff3 # [N] Assume chemical equilibrium
    τs1m1 = fill(Inf, N) # [N] No exchange
    τs1m2 = fill(Inf, N) # [N] No exchange
    τs1m3 = fill(Inf, N) # [N] No exchange
    τs2m1 = fill(Inf, N) # [N] No exchange
    τs2m2 = fill(Inf, N) # [N] No exchange
    τs2m3 = fill(Inf, N) # [N] No exchange
    τs3m1 = fill(Inf, N) # [N] No exchange
    τs3m2 = fill(Inf, N) # [N] No exchange
    τs3m3 = fill(Inf, N) # [N] No exchange
    τm1f1 = fill(Inf, N) # [N] No exchange
    τm1f2 = fill(Inf, N) # [N] No exchange
    τm1f3 = fill(Inf, N) # [N] No exchange
    τm2f1 = fill(Inf, N) # [N] No exchange
    τm2f2 = fill(Inf, N) # [N] No exchange
    τm2f3 = fill(Inf, N) # [N] No exchange
    τm3f1 = fill(Inf, N) # [N] No exchange
    τm3f2 = fill(Inf, N) # [N] No exchange
    τm3f3 = fill(Inf, N) # [N] No exchange
    τm1s1 = fill(Inf, N) # [N] No exchange
    τm1s2 = fill(Inf, N) # [N] No exchange
    τm1s3 = fill(Inf, N) # [N] No exchange
    τm2s1 = fill(Inf, N) # [N] No exchange
    τm2s2 = fill(Inf, N) # [N] No exchange
    τm2s3 = fill(Inf, N) # [N] No exchange
    τm3s1 = fill(Inf, N) # [N] No exchange
    τm3s2 = fill(Inf, N) # [N] No exchange
    τm3s3 = fill(Inf, N) # [N] No exchange
    frac = [[ff1[i], ff2[i], ff3[i], fs1[i], fs2[i], fs3[i],
             fm1[i], fm2[i], fm3[i]] for i = 1:N]
    T1 = [[T1f1[i], T1f2[i], T1f3[i], T1s1[i], T1s2[i], T1s3[i],
           T1m1[i], T1m2[i], T1m3[i]] for i = 1:N]
    T2 = [[T2f1[i], T2f2[i], T2f3[i], T2s1[i], T2s2[i], T2s3[i],
           T2m1[i], T2m2[i], T2m3[i]] for i = 1:N]
    Δωtmp = [[Δω[i] + Δωf[i], Δω[i] + Δωf[i], Δω[i] + Δωf[i],
              Δω[i], Δω[i], Δω[i], Δω[i], Δω[i], Δω[i]] for i = 1:N]
    τ = [[τf1f2[i], τf1f3[i], τf1s1[i], τf1s2[i], τf1s3[i], τf1m1[i], τf1m2[i], τf1m3[i],
          τf2f1[i], τf2f3[i], τf2s1[i], τf2s2[i], τf2s3[i], τf2m1[i], τf2m2[i], τf2m3[i],
          τf3f1[i], τf3f2[i], τf3s1[i], τf3s2[i], τf3s3[i], τf3m1[i], τf3m2[i], τf3m3[i],
          τs1f1[i], τs1f2[i], τs1f3[i], τs1s2[i], τs1s3[i], τs1m1[i], τs1m2[i], τs1m3[i],
          τs2f1[i], τs2f2[i], τs2f3[i], τs2s1[i], τs2s3[i], τs2m1[i], τs2m2[i], τs2m3[i],
          τs3f1[i], τs3f2[i], τs3f3[i], τs3s1[i], τs3s2[i], τs3m1[i], τs3m2[i], τs3m3[i],
          τm1f1[i], τm1f2[i], τm1f3[i], τm1s1[i], τm1s2[i], τm1s3[i], τm1m2[i], τm1m3[i],
          τm2f1[i], τm2f2[i], τm2f3[i], τm2s1[i], τm2s2[i], τm2s3[i], τm2m1[i], τm2m3[i],
          τm3f1[i], τm3f2[i], τm3f3[i], τm3s1[i], τm3s2[i], τm3s3[i], τm3m1[i], τm3m2[i]] for i = 1:N]
    if nspins == 1
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes) # [N][nechoes]
    else
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes,
                          nspins, ncycles) # [N][nechoes]
    end
    y = reduce(hcat, y) # [nechoes,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3, T1f1, T1f2, T1f3,
         T1s1, T1s2, T1s3, T1m1, T1m2, T1m3, T2f1, T2f2, T2f3, T2s1, T2s2, T2s3,
         T2m1, T2m2, T2m3, τf1f2, τf1f3, τf1s1, τf1s2, τf1s3, τf1m1, τf1m2,
         τf1m3, τf2f3, τf2s1, τf2s2, τf2s3, τf2m1, τf2m2, τf2m3, τf3s1, τf3s2,
         τf3s3, τf3m1, τf3m2, τf3m3, τs1s2, τs1s3, τs2s3, τm1m2, τm1m3, τm2m3,
         Δωf, Δω, κ]

    return (y, x)

end

function testdata_brainweb_9comp(scan::Symbol, tissues = (:wm, :gm);
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    Δω::AbstractArray{<:Real,2} = (genmap(trues(181, 217), 50, 0, 100) .- 50) .* 2π,
    κ::AbstractArray{<:Real,2} = genmap(trues(181, 217), 1, 0.8, 1.2),
    σ::Real = 2e-3,
    meseScanDesign = (TR = 1200, αex = π/2, αinv = π, TE = 10, nechoes = 32),
    D::Int = scan == :mese ? 32 : 11 # 11 STFR scans, 32 MESE scans
)

    # Read in the brain
    obj = readbrainweb(9)

    # Preallocate the output
    data = zeros(size(obj)..., D)

    # Select the tissues of interest
    (mask, M0, ff1, ff2, ff3, fs1, fs2, fs3, fm1, fm2, fm3, T1f1, T1f2, T1f3,
     T1s1, T1s2, T1s3, T1m1, T1m2, T1m3, T2f1, T2f2, T2f3, T2s1, T2s2, T2s3,
     T2m1, T2m2, T2m3, τf1f2, τf1f3, τf1s1, τf1s2, τf1s3, τf1m1, τf1m2,
     τf1m3, τf2f3, τf2s1, τf2s2, τf2s3, τf2m1, τf2m2, τf2m3, τf3s1, τf3s2,
     τf3s3, τf3m1, τf3m2, τf3m3, τs1s2, τs1s3, τs2s3, τm1m2, τm1m3, τm2m3, Δωf) =
        getparams(obj, tissues...)
    N = count(mask)

    # Simulate the scan
    if scan == :stfrblochsim

        (y,) = testdata_stfrblochsim_9comp(scanDesignFile = scanDesignFile, N = N,
            M0 = M0, ff1 = ff1, ff2 = ff2, ff3 = ff3, fs1 = fs1, fs2 = fs2,
            fs3 = fs3, fm1 = fm1, fm2 = fm2, fm3 = fm3, T1f1 = T1f1,
            T1f2 = T1f2, T1f3 = T1f3, T1s1 = T1s1, T1s2 = T1s2, T1s3 = T1s3,
            T1m1 = T1m1, T1m2 = T1m2, T1m3 = T1m3, T2f1 = T2f1, T2f2 = T2f2,
            T2f3 = T2f3, T2s1 = T2s1, T2s2 = T2s2, T2s3 = T2s3, T2m1 = T2m1,
            T2m2 = T2m2, T2m3 = T2m3, τf1f2 = τf1f2, τf1f3 = τf1f3,
            τf1s1 = τf1s1, τf1s2 = τf1s2, τf1s3 = τf1s3, τf1m1 = τf1m1,
            τf1m2 = τf1m2, τf1m3 = τf1m3, τf2f3 = τf2f3, τf2s1 = τf2s1,
            τf2s2 = τf2s2, τf2s3 = τf2s3, τf2m1 = τf2m1, τf2m2 = τf2m2,
            τf2m3 = τf2m3, τf3s1 = τf3s1, τf3s2 = τf3s2, τf3s3 = τf3s3,
            τf3m1 = τf3m1, τf3m2 = τf3m2, τf3m3 = τf3m3, τs1s2 = τs1s2,
            τs1s3 = τs1s3, τs2s3 = τs2s3, τm1m2 = τm1m2, τm1m3 = τm1m3,
            τm2m3 = τm2m3, Δωf = Δωf, Δω = Δω[mask], κ = κ[mask], σ = σ)

    elseif scan == :mese

        (y,) = testdata_meseblochsim_9comp(N = N,
            M0 = M0, ff1 = ff1, ff2 = ff2, ff3 = ff3, fs1 = fs1, fs2 = fs2,
            fs3 = fs3, fm1 = fm1, fm2 = fm2, fm3 = fm3, T1f1 = T1f1,
            T1f2 = T1f2, T1f3 = T1f3, T1s1 = T1s1, T1s2 = T1s2, T1s3 = T1s3,
            T1m1 = T1m1, T1m2 = T1m2, T1m3 = T1m3, T2f1 = T2f1, T2f2 = T2f2,
            T2f3 = T2f3, T2s1 = T2s1, T2s2 = T2s2, T2s3 = T2s3, T2m1 = T2m1,
            T2m2 = T2m2, T2m3 = T2m3, τf1f2 = τf1f2, τf1f3 = τf1f3,
            τf1s1 = τf1s1, τf1s2 = τf1s2, τf1s3 = τf1s3, τf1m1 = τf1m1,
            τf1m2 = τf1m2, τf1m3 = τf1m3, τf2f3 = τf2f3, τf2s1 = τf2s1,
            τf2s2 = τf2s2, τf2s3 = τf2s3, τf2m1 = τf2m1, τf2m2 = τf2m2,
            τf2m3 = τf2m3, τf3s1 = τf3s1, τf3s2 = τf3s2, τf3s3 = τf3s3,
            τf3m1 = τf3m1, τf3m2 = τf3m2, τf3m3 = τf3m3, τs1s2 = τs1s2,
            τs1s3 = τs1s3, τs2s3 = τs2s3, τm1m2 = τm1m2, τm1m3 = τm1m3,
            τm2m3 = τm2m3, Δωf = Δωf, Δω = Δω[mask], κ = κ[mask], σ = σ;
            meseScanDesign...)

    else

        throw(ArgumentError("invalid scan"))

    end

    # Reshape data
    for d = 1:D

        @view(data[:,:,d])[mask] = y[d,:]

    end

    return (data, obj, mask)

end

function testdata_stfrblochsim_4comp(;
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1},
    ff::AbstractArray{<:Real,1},
    fm::AbstractArray{<:Real,1},
    fa::AbstractArray{<:Real,1},
    T1f::AbstractArray{<:Real,1},
    T1s::AbstractArray{<:Real,1},
    T1m::AbstractArray{<:Real,1},
    T1a::AbstractArray{<:Real,1},
    T2f::AbstractArray{<:Real,1},
    T2s::AbstractArray{<:Real,1},
    T2m::AbstractArray{<:Real,1},
    T2a::AbstractArray{<:Real,1},
    τfs::AbstractArray{<:Real,1},
    τfm::AbstractArray{<:Real,1},
    τfa::AbstractArray{<:Real,1},
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0, 35), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3
)

    (Tfree, Tg, α, β, ϕ, TE) = load(scanDesignFile, "P")[1] # [D]
    fs = 1 .- ff .- fm .- fa # [N]
    τsf = τfs .* fs ./ ff # [N] Assume chemical equilibrium
    τaf = τfa .* fa ./ ff # [N] Assume chemical equilibrium
    τsm = fill(Inf, N) # [N] No exchange
    τsa = fill(Inf, N) # [N] No exchange
    τmf = fill(Inf, N) # [N] No exchange
    τms = fill(Inf, N) # [N] No exchange
    τma = fill(Inf, N) # [N] No exchange
    τas = fill(Inf, N) # [N] No exchange
    τam = fill(Inf, N) # [N] No exchange
    frac = permutedims([[ff[i], fs[i], fm[i], fa[i]] for i = 1:N])
    T1 = permutedims([[T1f[i], T1s[i], T1m[i], T1a[i]] for i = 1:N])
    T2 = permutedims([[T2f[i], T2s[i], T2m[i], T2a[i]] for i = 1:N])
    Δωtmp = permutedims([[Δω[i] + Δωf[i], Δω[i], Δω[i], Δω[i]] for i = 1:N])
    τ = permutedims([[τfs[i], τfm[i], τfa[i], τsf[i], τsm[i], τsa[i], τmf[i],
                      τms[i], τma[i], τaf[i], τas[i], τam[i]] for i = 1:N])
    y = stfrblochsim.(M0', frac, T1, T2, Δωtmp, τ, κ', Tfree, Tg, α, β, ϕ, TE) # [D,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff, fm, fa, T1f, T1s, T1m, T1a, T2f, T2s, T2m, T2a, τfs, τfm, τfa,
         Δωf, Δω, κ]

    return (y, x)

end

function testdata_meseblochsim_4comp(;
    TR::Real = 1200,
    αex::Real = π/2,
    αinv::Real = π,
    TE::Real = 10,
    nechoes::Integer = 32,
    N::Integer = 1000, # Number test points
    M0::AbstractArray{<:Real,1},
    ff::AbstractArray{<:Real,1},
    fm::AbstractArray{<:Real,1},
    fa::AbstractArray{<:Real,1},
    T1f::AbstractArray{<:Real,1},
    T1s::AbstractArray{<:Real,1},
    T1m::AbstractArray{<:Real,1},
    T1a::AbstractArray{<:Real,1},
    T2f::AbstractArray{<:Real,1},
    T2s::AbstractArray{<:Real,1},
    T2m::AbstractArray{<:Real,1},
    T2a::AbstractArray{<:Real,1},
    τfs::AbstractArray{<:Real,1},
    τfm::AbstractArray{<:Real,1},
    τfa::AbstractArray{<:Real,1},
    Δωf::AbstractArray{<:Real,1} = rand(Uniform(0, 35), N),
    Δω::AbstractArray{<:Real,1} = rand(Uniform(-50 * 2π, 50 * 2π), N),
    κ::AbstractArray{<:Real,1} = rand(Uniform(0.8, 1.2), N),
    σ::Real = 2e-3,
    nspins::Integer = 30, # Set to 1 to use ideal spoiling
    ncycles::Real = 1
)

    fs = 1 .- ff .- fm .- fa # [N]
    τsf = τfs .* fs ./ ff # [N] Assume chemical equilibrium
    τaf = τfa .* fa ./ ff # [N] Assume chemical equilibrium
    τsm = fill(Inf, N) # [N] No exchange
    τsa = fill(Inf, N) # [N] No exchange
    τmf = fill(Inf, N) # [N] No exchange
    τms = fill(Inf, N) # [N] No exchange
    τma = fill(Inf, N) # [N] No exchange
    τas = fill(Inf, N) # [N] No exchange
    τam = fill(Inf, N) # [N] No exchange
    frac = [[ff[i], fs[i], fm[i], fa[i]] for i = 1:N]
    T1 = [[T1f[i], T1s[i], T1m[i], T1a[i]] for i = 1:N]
    T2 = [[T2f[i], T2s[i], T2m[i], T2a[i]] for i = 1:N]
    Δωtmp = [[Δω[i] + Δωf[i], Δω[i], Δω[i], Δω[i]] for i = 1:N]
    τ = [[τfs[i], τfm[i], τfa[i], τsf[i], τsm[i], τsa[i], τmf[i],
          τms[i], τma[i], τaf[i], τas[i], τam[i]] for i = 1:N]
    if nspins == 1
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes) # [N][nechoes]
    else
        y = meseblochsim.(M0, frac, T1, T2, Δωtmp, τ, κ, TR, αex, αinv, TE, nechoes,
                          nspins, ncycles) # [N][nechoes]
    end
    y = reduce(hcat, y) # [nechoes,N]
    y .+= complex.(rand(Normal(0, σ), size(y)...), rand(Normal(0, σ), size(y)...))
    y = abs.(y)

    # Inspect SNR
    noise = testdata_noise(N = 1000, σ = σ)
    SNR = [mean(y[d,:]) / (std(noise) / sqrt(2 - π/2)) for d = 1:size(y,1)]
    @show(σ, SNR)

    # Package up the true values for the unknown (and known) parameters to return
    x = [M0, ff, fm, fa, T1f, T1s, T1m, T1a, T2f, T2s, T2m, T2a, τfs, τfm, τfa,
         Δωf, Δω, κ]

    return (y, x)

end

function testdata_brainweb_4comp(scan::Symbol, tissues = (:wm, :gm);
    scanDesignFile::String = modulepath("estimation/data/designAsorted.jld"),
    Δω::AbstractArray{<:Real,2} = (genmap(trues(181, 217), 50, 0, 100) .- 50) .* 2π,
    κ::AbstractArray{<:Real,2} = genmap(trues(181, 217), 1, 0.8, 1.2),
    σ::Real = 2e-3,
    meseScanDesign = (TR = 1200, αex = π/2, αinv = π, TE = 10, nechoes = 32),
    D::Int = scan == :mese ? 32 : 11 # 11 STFR scans, 32 MESE scans
)

    # Read in the brain
    obj = readbrainweb(4)

    # Preallocate the output
    data = zeros(size(obj)..., D)

    # Select the tissues of interest
    (mask, M0, ff, fm, fa, T1f, T1s, T1m, T1a, T2f, T2s, T2m, T2a, τfs, τfm,
     τfa, Δωf) = getparams(obj, tissues...)
    N = count(mask)

    # Simulate the scan
    if scan == :stfrblochsim

        (y,) = testdata_stfrblochsim_4comp(scanDesignFile = scanDesignFile,
            N = N, M0 = M0, ff = ff, fm = fm, fa = fa, T1f = T1f, T1s = T1s,
            T1m = T1m, T1a = T1a, T2f = T2f, T2s = T2s, T2m = T2m, T2a = T2a,
            τfs = τfs, τfm = τfm, τfa = τfa, Δωf = Δωf, Δω = Δω[mask],
            κ = κ[mask], σ = σ)

    elseif scan == :mese

        (y,) = testdata_meseblochsim_4comp(
            N = N, M0 = M0, ff = ff, fm = fm, fa = fa, T1f = T1f, T1s = T1s,
            T1m = T1m, T1a = T1a, T2f = T2f, T2s = T2s, T2m = T2m, T2a = T2a,
            τfs = τfs, τfm = τfm, τfa = τfa, Δωf = Δωf, Δω = Δω[mask],
            κ = κ[mask], σ = σ; meseScanDesign...)

    else

        throw(ArgumentError("invalid scan"))

    end

    # Reshape data
    for d = 1:D

        @view(data[:,:,d])[mask] = y[d,:]

    end

    return (data, obj, mask)

end
