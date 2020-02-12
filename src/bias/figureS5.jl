"""
    figureS5()

Create Supporting Information Figure S5. Figure S5 compares the unbiased and
biased Crámer-Rao Lower Bounds for scan design A for white matter tissue
parameter values as a function of the additional myelin water off-resonance
frequency Δωf.
"""
function figureS5()

    # Run figureS5data() if it hasn't run yet
    isfile(modulepath("bias/results/biased_crlb_vs_Dwf.jld")) ||
        figureS5data()

    (Δff, biased, unbiased) = load(modulepath("bias/results/biased_crlb_vs_Dwf.jld"),
                                   "Δff", "biased", "unbiased")
    pyplot()
    p = plot(title = "Unbiased vs Biased CRLB for WM",
             xlabel = L"Δω_{\mathrm{f}} \; \mathrm{ (Hz)}",
             ylabel = "CRLB of Standard Deviation of MWF",
             ylims = (-0.05, 1), yticks = 0:0.1:1, xticks = 0:5:35,
             foreground_color_grid = :lightgrey, gridalpha = 1.0,
             size = (600, 400), dpi = 300)
    plot!(p, Δff, unbiased, line = (:red),
          label = "Unbiased CRLB",
          marker = (:rect, :red))
    plot!(p, Δff, biased, line = (:blue),
          label = "Biased CRLB",
          marker = (:circle, :blue))
    display(p)
    gr()

end

"""
    figureS5data()

Create data for Supporting Information Figure S5.
"""
function figureS5data()

    # Make sure scan design file exists
    createdesignA_bias()

    P = load(modulepath("bias/data/designA.jld"), "P")
    Δff = 0:1:35
    biased = zeros(length(Δff))
    unbiased = zeros(length(Δff))
    f = (Δff, biased) -> crlb(0.77, 0.15, 400, 832, 20, 80, Δff * 2π, 0, 1,
                              P, biased)
    prog = Progress(length(Δff), dt = 0.1, barglyphs = BarGlyphs("[=> ]"))
    for i = 1:length(Δff)
        biased[i] = f(Δff[i], true)
        unbiased[i] = f(Δff[i], false)[2,2]
        next!(prog)
    end
    biased .= sqrt.(biased)
    unbiased .= sqrt.(unbiased)

    # Make sure results/ folder exists
    isdir(modulepath("bias/results")) || mkdir(modulepath("bias/results"))

    @save(modulepath("bias/results/biased_crlb_vs_Dwf.jld"), Δff, biased, unbiased)

end
