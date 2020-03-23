"""
    figure2()

Create Figure 2. Figure 2 plots the expected unbiased Crámer-Rao Lower Bound for
both scan design A and B as a function of the additional myelin water
off-resonance frequency Δωf.
"""
function figure2()

    # Run figure2data() if it hasn't run yet
    isfile(modulepath("scandesign/results/compare_designs_range.jld")) ||
        figure2data()

    (Δff, cost) = load(modulepath("scandesign/results/compare_designs_range.jld"), "Δff", "cost")
    p = plot(title = "Comparison of CRLBs of Scan Designs",
             xlabel = L"Δω_{\mathrm{f}} \; \mathrm{ (Hz)}",
             ylabel = "CRLB of Standard Deviation of MWF",
             foreground_color_grid = :lightgrey, gridalpha = 1.0,
             size = (600, 400), dpi = 300, reuse = false)
    plot!(p, Δff, cost[:,2], line = (:red),
          label = L"\mathrm{Design \; B: Ignore \;} Δω_{\mathrm{f}} \mathrm{\; during \; scan \; design}",
          marker = (:rect, :red))
    plot!(p, Δff, cost[:,1], line = (:blue),
          label = L"\mathrm{Design \; A: Optimize \; over \; range \; of \;} Δω_{\mathrm{f}}",
          marker = (:circle, :blue))
    display(p)

end

"""
    figure2data()

Create data for Figure 2.
"""
function figure2data()

    # Make sure scan design files exist
    createdesignA_scandesign()
    createdesignB_scandesign()

    PA = load(modulepath("scandesign/results/designA.jld"), "P")
    PB = load(modulepath("scandesign/results/designB.jld"), "P")
    (Σ, xPriors, νPriors, gradx) = getargs(Δωfdist = Normal(0 * 2π, eps()),
                                           Δωfnsamp = 1,
                                           Δωdist = Uniform(-50 * 2π, 50 * 2π),
                                           Δωnsamp = 10,
                                           σ = 4e-3)
    Δff = 0:1:40
    cost = zeros(length(Δff),2)
    prog = Progress(2length(Δff), dt = 0.1, barglyphs = BarGlyphs("[=> ]"))
    for i = 1:length(Δff)
        xPriors[7]["dist"] = Normal(Δff[i] * 2π, eps())
        cost[i,1] = expectedcost(PA, Σ, xPriors, νPriors, gradx, [x->0])
        next!(prog)
        cost[i,2] = expectedcost(PB, Σ, xPriors, νPriors, gradx, [x->0])
        next!(prog)
    end
    cost = sqrt.(cost)

    @save(modulepath("scandesign/results/compare_designs_range.jld"), Δff, cost)

end

function getargs(;
    D::Integer = 9,
    Tfreeval::Union{<:Real,Nothing} = 8.0,
    Tgval::Union{<:Real,Nothing} = 2.8,
    αval::Union{<:Real,Nothing} = nothing,
    βval::Union{<:Real,Nothing} = nothing,
    ϕval::Union{<:Real,Nothing} = nothing,
    Tfreespgrval::Union{<:Real,Nothing} = 8.0,
    Tgspgrval::Union{<:Real,Nothing} = 2.8,
    αspgrval::Union{<:Real,Nothing} = nothing,
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
)

    if αspgrval != 0
        Σ = [[σ], [σ], σ * ones(D)]
    else
        Σ = [σ * ones(D)]
    end

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

    (gradx, _) = getstfrgradfunctions(M0known = M0known, ffknown = ffknown,
                      T1fknown = T1fknown, T1sknown = T1sknown,
                      T2fknown = T2fknown, T2sknown = T2sknown,
                      Δωfknown = Δωfknown, Δωknown = Δωknown, κknown = κknown,
                      Tfreeval = Tfreeval, Tgval = Tgval, αval = αval,
                      βval = βval, ϕval = ϕval, Tfreespgrval = Tfreespgrval,
                      Tgspgrval = Tgspgrval, αspgrval = αspgrval, ΔTE = ΔTE)

    return (Σ, xPriors, νPriors, gradx)

end
