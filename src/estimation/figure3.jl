"""
    figure3()

Create Figure 3. Figure 3 compares different methods of accounting for the
additional myelin water off-resonance frequency Δωf. In one case, Δωf is
ignored in both the scan design and when training PERK. In another case, Δωf is
ignored in the scan design but modeled in the PERK training data. In the third
case, Δωf is modeled in both the scan design and in the PERK training data.
This simulation uses the two-compartment non-exchanging model.
"""
function figure3()

    # Run figure3data() if it hasn't run yet
    isfile(modulepath("estimation/results/compare_designs_range.jld")) || figure3data()

    (Δff, rmse) = load(modulepath("estimation/results/compare_designs_range.jld"), "Δff", "rmse")
    pyplot()
    p = plot(title = L"\mathrm{Comparison \; of \; How \; to \; Account \; for \;} Δω_{\mathrm{f}}",
             xlabel = L"Δω_{\mathrm{f}} \; \mathrm{ (Hz)}",
             ylabel = "RMSE of MWF Estimates",
             ylims = (0, 0.07), yticks = 0:0.01:0.07, xticks = 0:5:35,
             foreground_color_grid = :lightgrey, gridalpha = 1.0,
             size = (600, 400), dpi = 300)
    plot!(p, Δff, rmse[:,1], line = (:green),
          label = L"\mathrm{Ignore \;} Δω_{\mathrm{f}} \mathrm{\; in \; scan \; design \; and \; in \; training}",
          marker = (:diamond, :green))
    plot!(p, Δff, rmse[:,2], line = (:red),
          label = L"\mathrm{Account \; for \;} Δω_{\mathrm{f}} \mathrm{\; in \; training \; but \; ignore \; in \; scan \; design}",
          marker = (:rect, :red))
    plot!(p, Δff, rmse[:,3], line = (:blue),
          label = L"\mathrm{Account \; for \;} Δω_{\mathrm{f}} \mathrm{\; in \; scan \; design \; and \; in \; training}",
          marker = (:circle, :blue))
    display(p)
    gr()

end

"""
    figure3data()

Create data for Figure 3.
"""
function figure3data()

    # Make sure scan design files exist
    createdesignA_estimation()
    createdesignB_estimation()

    N = 1000
    Nnoise = 1000
    σ = 4e-3
    Δff = 0:1:35 # Hz
    rmse = zeros(length(Δff),3)
    prog = Progress(3length(Δff), dt = 0.1, barglyphs = BarGlyphs("[=> ]"))

    getdata = (scanDesignFile, Δωf) -> begin
        (y, x) = testdata_stfr(scanDesignFile = modulepath(scanDesignFile),
            N = N, T1f = rand(truncated(Normal(400, 80), 320, 480), N),
            T1s = rand(truncated(Normal(1000, 200), 800, 1200), N),
            T2f = rand(truncated(Normal(20, 4), 16, 24), N),
            T2s = rand(truncated(Normal(80, 16), 64, 96), N),
            Δωf = fill(Δωf, N), σ = σ)
        noise = testdata_noise(N = Nnoise, σ = σ)
        (ff, Δω, κ) = x[[2,8,9]]
        return (y, noise, Δω, κ, ff)
    end

    for i = 1:length(Δff)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata("estimation/data/designB.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ, ignoreΔωf = true,
            scanDesignFile = modulepath("estimation/data/designB.jld"))[1][2,:]
        rmse[i,1] = norm(ffhat - ff) / sqrt(N)
        next!(prog)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata("estimation/data/designB.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designB.jld"))[1][2,:]
        rmse[i,2] = norm(ffhat - ff) / sqrt(N)
        next!(prog)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata("estimation/data/designA.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designA.jld"))[1][2,:]
        rmse[i,3] = norm(ffhat - ff) / sqrt(N)
        next!(prog)

    end

    # Make sure results/ folder exists
    isdir(modulepath("estimation/results")) || mkdir(modulepath("estimation/results"))

    @save(modulepath("estimation/results/compare_designs_range.jld"), Δff, rmse)

end
