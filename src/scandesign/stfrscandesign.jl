"""
    scandesignA()

Create scan design A, the design that models the additional myelin water
off-resonance Δωf.
"""
function scandesignA()

    # Make sure results/ folder exists
    isdir(modulepath("scandesign/results")) || mkdir(modulepath("scandesign/results"))

    (P, cost, flag) = stfrscandesign(Δωdist = Uniform(-50 * 2π, 50 * 2π),
                    Δωnsamp = 10, Δωknown = true, maxtime = 20 * 60 * 60,
                    savefilename = modulepath("scandesign/results/designA.jld"),
                    description = "account for Δωf")

    displaydesignA()

    return (P, cost, flag)

end

"""
    scandesignB()

Create scan design B, the design that ignores the additional myelin water
off-resonance Δωf.
"""
function scandesignB()

    # Make sure results/ folder exists
    isdir(modulepath("scandesign/results")) || mkdir(modulepath("scandesign/results"))

    (P, cost, flag) = stfrscandesign(Δωdist = Uniform(-50 * 2π, 50 * 2π),
                    Δωnsamp = 10, Δωknown = true,
                    Δωfdist = Normal(0 * 2π, eps() * 2π),
                    Δωfnsamp = 1, Δωfknown = true, maxtime = 20 * 60 * 60,
                    savefilename = modulepath("scandesign/results/designB.jld"),
                    description = "ignore Δωf")

    displaydesignB()

    return (P, cost, flag)

end

function stfrscandesign(;
    D::Integer = 9, # Number of scans, excluding SPGR scans
    Tfreelb::AbstractArray{<:Real,1} = 6 * ones(D), # Only use lb and ub if param not fixed
    Tfreeub::AbstractArray{<:Real,1} = Inf * ones(D),
    Tglb::AbstractArray{<:Real,1} = 0 * ones(D),
    Tgub::AbstractArray{<:Real,1} = Inf * ones(D),
    αlb::AbstractArray{<:Real,1} = 1 * π / 180 * ones(D),
    αub::AbstractArray{<:Real,1} = 15 * π / 180 * ones(D),
    βlb::AbstractArray{<:Real,1} = 0 * π / 180 * ones(D),
    βub::AbstractArray{<:Real,1} = 15 * π / 180 * ones(D),
    ϕlb::AbstractArray{<:Real,1} = -π * ones(D),
    ϕub::AbstractArray{<:Real,1} = π * ones(D),
    Tfree::Union{<:Real,<:AbstractArray{<:Real,1}} = 8.0,
    Tg::Union{<:Real,<:AbstractArray{<:Real,1}} = 2.8,
    α::Union{<:Real,<:AbstractArray{<:Real,1}} = LinRange(maximum(αlb), minimum(αub), D),
    β::Union{<:Real,<:AbstractArray{<:Real,1}} = LinRange(maximum(βlb), minimum(βub), D),
    ϕ::Union{<:Real,<:AbstractArray{<:Real,1}} = LinRange(maximum(ϕlb), minimum(ϕub), D),
    Tfreespgrlb::AbstractArray{<:Real,1} = [Tfreelb[1]], # Only use lb and ub if param not fixed
    Tfreespgrub::AbstractArray{<:Real,1} = [Tfreeub[1]],
    Tgspgrlb::AbstractArray{<:Real,1} = [Tglb[1]],
    Tgspgrub::AbstractArray{<:Real,1} = [Tgub[1]],
    αspgrlb::AbstractArray{<:Real,1} = [5 * π / 180],
    αspgrub::AbstractArray{<:Real,1} = [5 * π / 180],
    Tfreespgr::Union{<:Real,<:AbstractArray{<:Real,1}} = 8.0, # If array should be length 1
    Tgspgr::Union{<:Real,<:AbstractArray{<:Real,1}} = 2.8, # If array should be length 1
    αspgr::Union{<:Real,<:AbstractArray{<:Real,1}} = [5 * π / 180], # Set to 0 if not using SPGR, if array should be length 1
    ΔTE::Real = 2.3,
    σ::Real = 3.8607e-4,
    M0dist::Distribution{<:Univariate,<:ValueSupport} = Normal(1, 0.2),
    M0nsamp::Integer = 1,
    M0cfvar::Bool = false,
    M0weight::Real = 0, # Weight only matters if param not known
    M0known::Bool = false,
    ffdist::Distribution{<:Univariate,<:ValueSupport} = Uniform(0.03, 0.31),
    ffnsamp::Integer = 5,
    ffcfvar::Bool = true,
    ffweight::Real = 1,
    ffknown::Bool = false,
    T1fdist::Distribution{<:Univariate,<:ValueSupport} = Normal(400, 80),
    T1fnsamp::Integer = 3,
    T1fcfvar::Bool = false,
    T1fweight::Real = 0,
    T1fknown::Bool = false,
    T1sdist::Distribution{<:Univariate,<:ValueSupport} = Normal(1000, 200),
    T1snsamp::Integer = 3,
    T1scfvar::Bool = false,
    T1sweight::Real = 0,
    T1sknown::Bool = false,
    T2fdist::Distribution{<:Univariate,<:ValueSupport} = Normal(20, 4),
    T2fnsamp::Integer = 3,
    T2fcfvar::Bool = false,
    T2fweight::Real = 0,
    T2fknown::Bool = false,
    T2sdist::Distribution{<:Univariate,<:ValueSupport} = Normal(80, 16),
    T2snsamp::Integer = 3,
    T2scfvar::Bool = false,
    T2sweight::Real = 0,
    T2sknown::Bool = false,
    Δωfdist::Distribution{<:Univariate,<:ValueSupport} = Uniform(5 * 2π, 35 * 2π),
    Δωfnsamp::Integer = 5,
    Δωfcfvar::Bool = false,
    Δωfweight::Real = 0,
    Δωfknown::Bool = false,
    Δωdist::Distribution{<:Univariate,<:ValueSupport} = Normal(0 * 2π, eps() * 2π),
    Δωnsamp::Integer = 1,
    Δωcfvar::Bool = false,
    Δωweight::Real = 0,
    Δωknown::Bool = true,
    κdist::Distribution{<:Univariate,<:ValueSupport} = Uniform(0.8, 1.2),
    κnsamp::Integer = 3,
    κcfvar::Bool = false,
    κweight::Real = 0,
    κknown::Bool = true,
    maxtime::Real = 16 * 60 * 60, # 16 hours
    algorithm::Symbol = :G_MLSL,
    localalgorithm::Symbol = :LD_SLSQP,
    verbosecount::Integer = 10,
    savefilename::String = modulepath("scandesign/results/designtmp.jld"),
    description::String = "", # Explain what went into the scan design
    resetRNG::Bool = true
)

    if resetRNG
        Random.seed!(0)
    end

    usespgr = αspgr != 0

    # Helper functions
    stfrhelp = (a, b, c, d, e) -> begin
        tmp = hcat((Tfree isa Real ? zeros(D, 0) : a),
                   (Tg    isa Real ? zeros(D, 0) : b),
                   (α     isa Real ? zeros(D, 0) : c),
                   (β     isa Real ? zeros(D, 0) : d),
                   (ϕ     isa Real ? zeros(D, 0) : e))
        [[tmp[:,i] for i = 1:size(tmp,2)]]
    end
    if usespgr
        spgrhelp = (a, b, c) -> begin
            tmp = hcat((Tfreespgr isa Real ? zeros(1, 0) : a),
                       (Tgspgr    isa Real ? zeros(1, 0) : b),
                       (αspgr     isa Real ? zeros(1, 0) : c))
            [[tmp[:,i] for i = 1:size(tmp,2)]]
        end
    end

    # Set up initial scan design
    P0stfr = stfrhelp(Tfree, Tg, α, β, ϕ)
    if usespgr
        P0spgr = spgrhelp(Tfreespgr, Tgspgr, αspgr)
        P0 = [P0spgr; P0spgr; P0stfr]
    else
        P0 = P0stfr
    end

    # Set up the cost function
    # Set up the scan standard deviation (assumed to be constant across all scans)
    if usespgr
        Σ = [[σ], [σ], σ * ones(D)]
    else
        Σ = [σ * ones(D)]
    end

    # Set up prior information
    M0  = Dict("dist" => M0dist,  "nsamp" => M0nsamp,  "cfvar" => M0cfvar,  "weight" => M0weight)
    ff  = Dict("dist" => ffdist,  "nsamp" => ffnsamp,  "cfvar" => ffcfvar,  "weight" => ffweight)
    T1f = Dict("dist" => T1fdist, "nsamp" => T1fnsamp, "cfvar" => T1fcfvar, "weight" => T1fweight)
    T1s = Dict("dist" => T1sdist, "nsamp" => T1snsamp, "cfvar" => T1scfvar, "weight" => T1sweight)
    T2f = Dict("dist" => T2fdist, "nsamp" => T2fnsamp, "cfvar" => T2fcfvar, "weight" => T2fweight)
    T2s = Dict("dist" => T2sdist, "nsamp" => T2snsamp, "cfvar" => T2scfvar, "weight" => T2sweight)
    Δωf = Dict("dist" => Δωfdist, "nsamp" => Δωfnsamp, "cfvar" => Δωfcfvar, "weight" => Δωfweight)
    Δω  = Dict("dist" => Δωdist,  "nsamp" => Δωnsamp,  "cfvar" => Δωcfvar,  "weight" => Δωweight)
    κ   = Dict("dist" => κdist,   "nsamp" => κnsamp,   "cfvar" => κcfvar,   "weight" => κweight)
    xPriors = Array{Dict{String,Any},1}()
    νPriors = Array{Dict{String,Any},1}()
    if M0known  push!(νPriors, M0)  else push!(xPriors, M0)  end
    if ffknown  push!(νPriors, ff)  else push!(xPriors, ff)  end
    if T1fknown push!(νPriors, T1f) else push!(xPriors, T1f) end
    if T1sknown push!(νPriors, T1s) else push!(xPriors, T1s) end
    if T2fknown push!(νPriors, T2f) else push!(xPriors, T2f) end
    if T2sknown push!(νPriors, T2s) else push!(xPriors, T2s) end
    if Δωfknown push!(νPriors, Δωf) else push!(xPriors, Δωf) end
    if Δωknown  push!(νPriors, Δω)  else push!(xPriors, Δω)  end
    if κknown   push!(νPriors, κ)   else push!(xPriors, κ)   end

    # Get the gradient functions
    (gradx, gradxp) = getstfrgradfunctions(M0known = M0known, ffknown = ffknown,
                        T1fknown = T1fknown, T1sknown = T1sknown,
                        T2fknown = T2fknown, T2sknown = T2sknown,
                        Δωfknown = Δωfknown, Δωknown = Δωknown, κknown = κknown,
                        Tfreeval = Tfree isa Real ? Tfree : nothing,
                        Tgval    = Tg    isa Real ? Tg    : nothing,
                        αval     = α     isa Real ? α     : nothing,
                        βval     = β     isa Real ? β     : nothing,
                        ϕval     = ϕ     isa Real ? ϕ     : nothing,
                        Tfreespgrval = Tfreespgr isa Real ? Tfreespgr : nothing,
                        Tgspgrval    = Tgspgr    isa Real ? Tgspgr    : nothing,
                        αspgrval     = αspgr     isa Real ? αspgr     : nothing,
                        ΔTE = ΔTE)

    costfun = (P, computegrad) -> expectedcost(P, Σ, xPriors, νPriors, gradx,
                                            gradxp, computegrad = computegrad)

    # Set up lower and upper bounds on the scan parameters
    lbstfr = stfrhelp(Tfreelb, Tglb, αlb, βlb, ϕlb)
    if usespgr
        lbspgr = spgrhelp(Tfreespgrlb, Tgspgrlb, αspgrlb)
        lb = [lbspgr; lbspgr; lbstfr]
    else
        lb = lbstfr
    end
    ubstfr = stfrhelp(Tfreeub, Tgub, αub, βub, ϕub)
    if usespgr
        ubspgr = spgrhelp(Tfreespgrub, Tgspgrub, αspgrub)
        ub = [ubspgr; ubspgr; ubstfr]
    else
        ub = ubstfr
    end

    # Run the optimization
    (P, cost, flag) = scandesign(P0, costfun, lb = lb, ub = ub,
                                 maxtime = maxtime, algorithm = algorithm,
                                 localalgorithm = localalgorithm,
                                 verbosecount = verbosecount)

    # Save the results
    @save(savefilename, P, cost, flag, description)

    return (P, cost, flag)

end

"""
    getstfrgradfunctions(; M0known, ffknown, T1fknown, T1sknown, T2fknown,
        T2sknown, Δωfknown, Δωknown, κknown, Tfreeval, Tgval, αval, βval, ϕval,
        Tfreespgrval, Tgspgrval, αspgrval, ΔTE)

Generate functions to be used to evaluate the gradient of the STFR signal
model.

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
- `Tfreeval::Union{<:Real,Nothing} = nothing`: Value to use for Tfree if fixed;
  set to `nothing` if Tfree is not fixed
- `Tgval::Union{<:Real,Nothing} = nothing`: Value to use for Tg if fixed; set to
  `nothing` if Tg is not fixed
- `αval::Union{<:Real,Nothing} = nothing`: Value to use for α if fixed; set to
  `nothing` if α is not fixed
- `βval::Union{<:Real,Nothing} = nothing`: Value to use for β if fixed; set to
  `nothing` if β is not fixed
- `ϕval::Union{<:Real,Nothing} = nothing`: Value to use for ϕ if fixed; set to
  `nothing` if ϕ is not fixed
- `Tfreespgrval::Union{<:Real,Nothing} = 8.0`: Tfree to use for SPGR scans; does
  not include ΔTE; set to `nothing` if Tfreespgrval is not fixed
- `Tgspgrval::Union{<:Real,Nothing} = 2.8`: Tg to use for SPGR scans; set to
  `nothing` if Tgspgrval is not fixed
- `αspgrval::Union{<:Real,Nothing} = 0`: Flip angle to use for SPGR scans; set
  to `0` if not using SPGR scans; set to `nothing` if αspgrval is not fixed
- `ΔTE::Real = 2.3`: Difference in echo time between the two SPGR scans

# Return
- `gradx::Array{Function,1}`: Each entry of `gradx` takes latent, known, and
  scan parameters and outputs the partial derivatives of the STFR signal model
  with respect to the latent parameters; `length(gradx) == 3` if SPGR scans
  are used, in which case `gradx[1]` and `gradx[2]` refer to the SPGR scans
- `gradxp::Array{Function,1}`: Each entry of `gradxp` takes latent, known, and
  scan parameters and outputs the mixed partial derivatives of the STFR signal
  model with respect to the latent parameters and the scan parameters;
  `length(gradxp) == 3` if SPGR scans are used, in which case `gradxp[1]` and
  `gradxp[2]` refer to the SPGR scans
"""
function getstfrgradfunctions(;
    M0known::Bool = false,
    ffknown::Bool = false,
    T1fknown::Bool = false,
    T1sknown::Bool = false,
    T2fknown::Bool = false,
    T2sknown::Bool = false,
    Δωfknown::Bool = false,
    Δωknown::Bool = true,
    κknown::Bool = true,
    Tfreeval::Union{<:Real,Nothing} = nothing,
    Tgval::Union{<:Real,Nothing} = nothing,
    αval::Union{<:Real,Nothing} = nothing,
    βval::Union{<:Real,Nothing} = nothing,
    ϕval::Union{<:Real,Nothing} = nothing,
    Tfreespgrval::Union{<:Real,Nothing} = 8.0,
    Tgspgrval::Union{<:Real,Nothing} = 2.8,
    αspgrval::Union{<:Real,Nothing} = 0,
    ΔTE::Real = 2.3
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
    Tfree = isnothing(Tfreeval) ? symbols("Tfree", positive    = true) : Tfreeval
    Tg    = isnothing(Tgval)    ? symbols("Tg",    positive    = true) : Tgval
    α     = isnothing(αval)     ? symbols("alpha", positive    = true) : αval
    β     = isnothing(βval)     ? symbols("beta",  nonnegative = true) : βval
    ϕ     = isnothing(ϕval)     ? symbols("phi",   positive    = true) : ϕval
    # SPGR scan parameters
    usespgr = αspgrval != 0
    if usespgr
        Tfreespgr = isnothing(Tfreespgrval) ? symbols("Tfreespgr", positive = true) : Tfreespgrval
        Tgspgr    = isnothing(Tgspgrval)    ? symbols("Tgspgr",    positive = true) : Tgspgrval
        αspgr     = isnothing(αspgrval)     ? symbols("alphaspgr", positive = true) : αspgrval
    end

    # Make list of latent parameters
    latent = Array{SymPy.Sym,1}()
    if !M0known  push!(latent, M0)  end
    if !ffknown  push!(latent, ff)  end
    if !T1fknown push!(latent, T1f) end
    if !T1sknown push!(latent, T1s) end
    if !T2fknown push!(latent, T2f) end
    if !T2sknown push!(latent, T2s) end
    if !Δωfknown push!(latent, Δωf) end
    if !Δωknown  push!(latent, Δω)  end
    if !κknown   push!(latent, κ)   end

    # Make list of non-fixed scan parameters
    scan = Array{SymPy.Sym,1}()
    if isnothing(Tfreeval) push!(scan, Tfree) end
    if isnothing(Tgval)    push!(scan, Tg)    end
    if isnothing(αval)     push!(scan, α)     end
    if isnothing(βval)     push!(scan, β)     end
    if isnothing(ϕval)     push!(scan, ϕ)     end

    # Make list of non-fixed SPGR scan parameters
    if usespgr
        scanspgr = Array{SymPy.Sym,1}()
        if isnothing(Tfreespgrval) push!(scanspgr, Tfreespgr) end
        if isnothing(Tgspgrval)    push!(scanspgr, Tgspgr)    end
        if isnothing(αspgrval)     push!(scanspgr, αspgr)     end
    end

    # Create the signal equation
    n = exp(-Tg/T1f) * (1 - exp(-Tfree/T1f)) * cos(κ * β) + (1 - exp(-Tg/T1f))
    d = 1 - exp(-Tg/T1f - Tfree/T2f) * sin(κ * α) * sin(κ * β) *
        cos((Δωf + Δω) * Tfree/1000 - ϕ) - exp(-Tg/T1f - Tfree/T1f) *
        cos(κ * α) * cos(κ * β)

    Mf = M0 * sin(κ * α) * n / d
    Ms = subs(Mf, T1f => T1s, T2f => T2s, Δωf => 0)
    M = ff * Mf + (1 - ff) * Ms

    # Create the SPGR signal equation
    if usespgr
        TfΔTE = Tfreespgr + ΔTE
        nspgr1 = exp(-Tgspgr/T1f) * (1 - exp(-TfΔTE/T1f)) + (1 - exp(-Tgspgr/T1f))
        dspgr1 = 1 - exp(-Tgspgr/T1f - TfΔTE/T1f) * cos(κ * αspgr)

        Mfspgr1 = M0 * sin(κ * αspgr) * nspgr1 / dspgr1
        Msspgr1 = subs(Mfspgr1, T1f => T1s, T2f => T2s, Δωf => 0)
        Mspgr1 = ff * Mfspgr1 + (1 - ff) * Msspgr1

        nspgr2 = exp(-Tgspgr/T1f) * (1 - exp(-TfΔTE/T1f)) + (1 - exp(-Tgspgr/T1f))
        dspgr2 = 1 - exp(-Tgspgr/T1f - TfΔTE/T1f) * cos(κ * αspgr)

        Mfspgr2 = M0 * sin(κ * αspgr) * nspgr2 / dspgr2
        Msspgr2 = subs(Mfspgr2, T1f => T1s, T2f => T2s, Δωf => 0)
        Mspgr2 = ff * Mfspgr2 + (1 - ff) * Msspgr2
    end

    # Calculate the gradients
    gradx_sym  = [diff(M, x) for x in latent]
    gradxp_sym = [diff(g, x) for g in gradx_sym, x in scan]
    if usespgr
        gradxspgr1_sym  = [diff(Mspgr1, x) for x in latent]
        gradxpspgr1_sym = [diff(g, x) for g in gradxspgr1_sym, x in scanspgr]
        gradxspgr2_sym  = [diff(Mspgr2, x) for x in latent]
        gradxpspgr2_sym = [diff(g, x) for g in gradxspgr2_sym, x in scanspgr]
    end

    # Convert to Julia functions
    gradx_arr  = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, scan...])
                  for g in gradx_sym]
    gradxp_arr = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, scan...])
                  for g in gradxp_sym]
    gradx  = (args...) -> [gradx_arr[i].(args...)    for i in 1:length(latent)]
    gradxp = (args...) -> [gradxp_arr[i,j].(args...) for i in 1:length(latent),
                                                         j in 1:length(scan)]
    if usespgr
        gradxspgr1_arr  = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ,
                                        scanspgr...]) for g in gradxspgr1_sym]
        gradxpspgr1_arr = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ,
                                        scanspgr...]) for g in gradxpspgr1_sym]
        gradxspgr1  = (args...) -> [gradxspgr1_arr[i].(args...)
                                    for i in 1:length(latent)]
        gradxpspgr1 = (args...) -> [gradxpspgr1_arr[i,j].(args...)
                                    for i in 1:length(latent),
                                        j in 1:length(scanspgr)]

        gradxspgr2_arr  = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ,
                                        scanspgr...]) for g in gradxspgr2_sym]
        gradxpspgr2_arr = [lambdify(g, [M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ,
                                        scanspgr...]) for g in gradxpspgr2_sym]
        gradxspgr2  = (args...) -> [gradxspgr2_arr[i].(args...)
                                    for i in 1:length(latent)]
        gradxpspgr2 = (args...) -> [gradxpspgr2_arr[i,j].(args...)
                                    for i in 1:length(latent),
                                        j in 1:length(scanspgr)]
    end

    if usespgr
        return ([gradxspgr1, gradxspgr2, gradx], [gradxpspgr1, gradxpspgr2, gradxp])
    else
        return ([gradx], [gradxp])
    end

end
