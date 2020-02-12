"""
    figureS6()

Create Supporting Information Figure S6. Figure S6 plots the bias of the
proposed STFR3-PERK estimator versus true MWF value for low and high SNR.
"""
function figureS6()

    # Run figureS6data() if it hasn't run yet
    isfile(modulepath("bias/results/expectation_blochsim_vs_MWF.jld")) ||
        figureS6data()

    (ff, expectation_lowSNR, expectation_highSNR) =
        load(modulepath("bias/results/expectation_blochsim_vs_MWF.jld"),
                        "ff", "expectation_lowSNR", "expectation_highSNR")
    pyplot()
    p = plot(ff[[1,end]], ff[[1,end]], line = (:green, :dash), label = "",
             title = "Expected MWF Estimate for WM",
             xlabel = "MWF",
             ylabel = "Expected MWF Estimate",
             ylims = (-0.015, 0.3), yticks = 0:0.05:0.3,
             xticks = 0:0.05:0.3,
             aspect_ratio = :equal,
             foreground_color_grid = :lightgrey, gridalpha = 1.0,
             size = (440, 400), dpi = 300)
    plot!(p, ff, expectation_lowSNR, line = (:blue),
          label = "Low SNR",
          marker = (:circle, :blue))
    plot!(p, ff, expectation_highSNR, line = (:orange),
          label = "High SNR",
          marker = (:square, :orange))
    display(p)
    gr()

end

"""
    figureS6data()

Create data for Supporting Information Figure S6.
"""
function figureS6data()

    # Make sure scan design file exists
    createdesignA_bias()

    P = load(modulepath("bias/data/designA.jld"), "P")
    ff = 0:0.01:0.3

    expectation_lowSNR = map(ff -> stfrblochsimperkexpectation(0.77, ff, 0.1,
        400, 832, 1000, 20, 80, 0.02, 100, 50, 15 * 2π, 0, 1, P), ff)

    expectation_highSNR = map(ff -> stfrblochsimperkexpectation(0.77, ff, 0.1,
        400, 832, 1000, 20, 80, 0.02, 100, 50, 15 * 2π, 0, 1, P, σ = 3.8607e-4),
        ff)

    # Make sure results/ folder exists
    isdir(modulepath("bias/results")) || mkdir(modulepath("bias/results"))

    @save(modulepath("bias/results/expectation_blochsim_vs_MWF.jld"), ff,
        expectation_lowSNR, expectation_highSNR)

end
