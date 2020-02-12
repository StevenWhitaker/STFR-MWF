"""
    figureS1()

Create Supporting Information Figure S1. Figure S1 compares the MWF estimates
from the two scan designs for both white matter and gray matter tissue
parameters. This simulation uses the two-compartment non-exchanging model.
"""
function figureS1()

    # Run figureS1data() if it hasn't run yet
    isfile(modulepath("estimation/results/compare_designs_wmgm.jld")) || figureS1data()

    (Δff, rmse) = load(modulepath("estimation/results/compare_designs_wmgm.jld"), "Δff", "rmse")
    pyplot()
    p = plot(title = "RMSE in White Matter and Gray Matter",
             xlabel = L"Δω_{\mathrm{f}} \; \mathrm{ (Hz)}",
             ylabel = "RMSE of MWF Estimates",
             ylims = (0, 0.05), xticks = 0:5:35, legend = :bottomleft,
             foreground_color_grid = :lightgrey, gridalpha = 1.0,
             size = (600, 400), dpi = 300)
    plot!(p, Δff, rmse[:,3], line = (:blue, :dot),
          label = "Design A: Gray Matter",
          marker = (:circle, :blue))
    plot!(p, Δff, rmse[:,4], line = (:red, :dot),
          label = "Design B: Gray Matter",
          marker = (:rect, :red))
    plot!(p, Δff, rmse[:,2], line = (:red),
          label = "Design B: White Matter",
          marker = (:rect, :red))
    plot!(p, Δff, rmse[:,1], line = (:blue),
          label = "Design A: White Matter",
          marker = (:circle, :blue))
    display(p)

end

"""
    figureS1data()

Create data for Supporting Information Figure S1.
"""
function figureS1data()

    # Make sure scan design files exist
    createdesignA_estimation()
    createdesignB_estimation()

    N = 1000
    Nnoise = 1000
    σ = 4e-3
    Δff = 0:1:35
    rmse = zeros(length(Δff),4)
    avg = zeros(length(Δff),4)
    stdev = zeros(length(Δff),4)
    prog = Progress(4length(Δff), dt = 0.1, barglyphs = BarGlyphs("[=> ]"))

    getdata_wm = (scanDesignFile, Δωf) -> begin
        (y, x) = testdata_stfr(scanDesignFile = modulepath(scanDesignFile),
            N = N, M0 = fill(0.77, N), ff = fill(0.17, N),
            T1f = fill(400, N), T1s = fill(832, N), T2f = fill(20, N),
            T2s = fill(80, N), Δωf = fill(Δωf, N), σ = σ)
        noise = testdata_noise(N = Nnoise, σ = σ)
        (ff, Δω, κ) = x[[2,8,9]]
        return (y, noise, Δω, κ, ff)
    end

    getdata_gm = (scanDesignFile, Δωf) -> begin
        (y, x) = testdata_stfr(scanDesignFile = modulepath(scanDesignFile),
            N = N, M0 = fill(0.86, N), ff = fill(0.03, N),
            T1f = fill(500, N), T1s = fill(1331, N), T2f = fill(20, N),
            T2s = fill(80, N), Δωf = fill(Δωf, N), σ = σ)
        noise = testdata_noise(N = Nnoise, σ = σ)
        (ff, Δω, κ) = x[[2,8,9]]
        return (y, noise, Δω, κ, ff)
    end

    for i = 1:length(Δff)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata_wm("estimation/data/designA.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designA.jld"))[1][2,:]
        rmse[i,1] = norm(ffhat - ff) / sqrt(N)
        avg[i,1] = mean(ffhat)
        stdev[i,1] = std(ffhat)
        next!(prog)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata_wm("estimation/data/designB.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designB.jld"))[1][2,:]
        rmse[i,2] = norm(ffhat - ff) / sqrt(N)
        avg[i,2] = mean(ffhat)
        stdev[i,2] = std(ffhat)
        next!(prog)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata_gm("estimation/data/designA.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designA.jld"))[1][2,:]
        rmse[i,3] = norm(ffhat - ff) / sqrt(N)
        avg[i,3] = mean(ffhat)
        stdev[i,3] = std(ffhat)
        next!(prog)

        Random.seed!(20190824)
        (y, noise, Δω, κ, ff) = getdata_gm("estimation/data/designB.jld", Δff[i] * 2π)
        ffhat = stfrperk(y, noise, Δω = Δω, κ = κ,
            scanDesignFile = modulepath("estimation/data/designB.jld"))[1][2,:]
        rmse[i,4] = norm(ffhat - ff) / sqrt(N)
        avg[i,4] = mean(ffhat)
        stdev[i,4] = std(ffhat)
        next!(prog)

    end

    # Make sure results/ folder exists
    isdir(modulepath("estimation/results")) || mkdir(modulepath("estimation/results"))

    @save(modulepath("estimation/results/compare_designs_wmgm.jld"), Δff, rmse)

end
