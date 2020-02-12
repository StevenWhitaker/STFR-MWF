"""
    figure5()

Create Figure 5. Figure 5 compares various MWF estimation techniques using white
matter and gray matter tissues with a nine-compartment model with exchange.
Also, display values analogous to those found in Table 3.
"""
function figure5()

    # Run figure5data() if it hasn't run yet
    (isfile(modulepath("estimation/results/brainweb_9comp_STFR2-PERK.jld")) &&
     isfile(modulepath("estimation/results/brainweb_9comp_STFR3-PERK.jld")) &&
     isfile(modulepath("estimation/results/brainweb_9comp_MESE-NNLS.jld")) &&
     isfile(modulepath("estimation/results/brainweb_9comp_MESE-PERK.jld")) &&
     isfile(modulepath("estimation/results/brainweb_9comp_STFR3-PERK-JE.jld"))) || figure4table3data()

    (ff1, ffmap1, rmse_wm1, rmse_gm1, mean_wm1, mean_gm1, std_wm1, std_gm1, t1) =
        load(modulepath("estimation/results/brainweb_9comp_STFR2-PERK.jld"), "ff", "ffmap",
        "rmse_wm", "rmse_gm", "mean_wm", "mean_gm", "std_wm", "std_gm", "t")
    (ff2, ffmap2, rmse_wm2, rmse_gm2, mean_wm2, mean_gm2, std_wm2, std_gm2, t2) =
        load(modulepath("estimation/results/brainweb_9comp_STFR3-PERK.jld"), "ff", "ffmap",
        "rmse_wm", "rmse_gm", "mean_wm", "mean_gm", "std_wm", "std_gm", "t")
    (ff3, ffmap3, rmse_wm3, rmse_gm3, mean_wm3, mean_gm3, std_wm3, std_gm3, t3) =
        load(modulepath("estimation/results/brainweb_9comp_MESE-NNLS.jld"), "ff", "ffmap",
        "rmse_wm", "rmse_gm", "mean_wm", "mean_gm", "std_wm", "std_gm", "t")
    (ff4, ffmap4, rmse_wm4, rmse_gm4, mean_wm4, mean_gm4, std_wm4, std_gm4, t4) =
        load(modulepath("estimation/results/brainweb_9comp_MESE-PERK.jld"), "ff", "ffmap",
        "rmse_wm", "rmse_gm", "mean_wm", "mean_gm", "std_wm", "std_gm", "t")
    (ff5, ffmap5, rmse_wm5, rmse_gm5, mean_wm5, mean_gm5, std_wm5, std_gm5, t5) =
        load(modulepath("estimation/results/brainweb_9comp_STFR3-PERK-JE.jld"), "ff",
        "ffmap", "rmse_wm", "rmse_gm", "mean_wm", "mean_gm", "std_wm",
        "std_gm", "t")

    # Ground truth MWF
    @assert ff1 == ff2 == ff3 == ff4 == ff5
    ff = ff1

    # Display MWF maps
    ff = [ffmap1 ff; ffmap5 ffmap2; ffmap4 ffmap3]
    p = heatmap(ff', color = :fire, aspect_ratio = :equal,
        xticks = [], yticks = [], clims = (0, 0.3),
        dpi = 300, size = (900, 600))
    display(p)

    # Helper functions for displaying statistics
    roundall = (rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t) -> begin
        (round(rmse_wm, digits = 3),
         round(rmse_gm, digits = 3),
         round(mean_wm, digits = 3),
         round(mean_gm, digits = 3),
         round(std_wm,  digits = 3),
         round(std_gm,  digits = 3),
         round(t,       digits = 1))
    end

    # Display stats
    (rmse_wm1, rmse_gm1, mean_wm1, mean_gm1, std_wm1, std_gm1, t1) =
        roundall(rmse_wm1, rmse_gm1, mean_wm1, mean_gm1, std_wm1, std_gm1, t1)

    (rmse_wm2, rmse_gm2, mean_wm2, mean_gm2, std_wm2, std_gm2, t2) =
        roundall(rmse_wm2, rmse_gm2, mean_wm2, mean_gm2, std_wm2, std_gm2, t2)

    (rmse_wm3, rmse_gm3, mean_wm3, mean_gm3, std_wm3, std_gm3, t3) =
        roundall(rmse_wm3, rmse_gm3, mean_wm3, mean_gm3, std_wm3, std_gm3, t3)

    (rmse_wm4, rmse_gm4, mean_wm4, mean_gm4, std_wm4, std_gm4, t4) =
        roundall(rmse_wm4, rmse_gm4, mean_wm4, mean_gm4, std_wm4, std_gm4, t4)

    (rmse_wm5, rmse_gm5, mean_wm5, mean_gm5, std_wm5, std_gm5, t5) =
        roundall(rmse_wm5, rmse_gm5, mean_wm5, mean_gm5, std_wm5, std_gm5, t5)

    println("              RMSE (WM)  Mean (WM)  Std (WM)  RMSE (GM)  Mean (GM)  Std (GM)  Time (s)")
    println("STFR2-PERK:    $rmse_wm1      $mean_wm1      $std_wm1     $rmse_gm1      $mean_gm1      $std_gm1     $t1")
    println("STFR3-PERK:    $rmse_wm2      $mean_wm2      $std_wm2     $rmse_gm2      $mean_gm2      $std_gm2     $t2")
    println("STFR3-PERK-JE: $rmse_wm5      $mean_wm5      $std_wm5     $rmse_gm5      $mean_gm5      $std_gm5     $t5")
    println("MESE-NNLS:     $rmse_wm3      $mean_wm3      $std_wm3     $rmse_gm3      $mean_gm3      $std_gm3     $t3")
    println("MESE-PERK:     $rmse_wm4      $mean_wm4      $std_wm4     $rmse_gm4      $mean_gm4      $std_gm4     $t4")

end

"""
    figure5data()

Create data for Figure 5.
"""
function figure5data()

    # Noise standard deviation
    σ = 4e-3

    # Make sure scan design file exists and create it if not
    createdesignAsorted_estimation()

    # Make sure data exists and download it if not
    isdir(modulepath("estimation/data/invivo")) || mkpath(modulepath("estimation/data/invivo"))
    isfile(modulepath("estimation/data/invivo/recon.jld")) || error("Download estimation/data/invivo/recon.jld")
    isfile(modulepath("estimation/data/invivo/b1slice5.mat")) || error("Download estimation/data/invivo/b1slice5.mat")

    # Make sure results/ folder exists
    isdir(modulepath("estimation/results")) || mkdir(modulepath("estimation/results"))

    # B0 and B1+ maps
    datapath = modulepath("estimation/data/invivo")
    b0 = load(joinpath(datapath, "recon.jld"), "b0map")[41:161,180:-1:21,1]
    Δω = (clamp.(imresize(b0, 181, 217), -40, 40) ./ (4/3)) .* 2π
    b1 = matread(joinpath(datapath, "b1slice5.mat"))["map"][341:461,160:-1:41]
    κ = clamp.(imresize(b1, 181, 217), 0.8, 1.2)

    # Helper functions for computing statistics
    rmse = (tissue, obj, ffmap, ff) -> begin
        m = getmask(obj, tissue)
        n = count(m)
        return norm(ffmap[m] - ff[m]) / sqrt(n)
    end

    avg = (tissue, obj, ffmap) -> begin
        m = getmask(obj, tissue)
        return mean(ffmap[m])
    end

    stdev = (tissue, obj, ffmap) -> begin
        m = getmask(obj, tissue)
        return std(ffmap[m])
    end

    # STFR2-PERK
    Random.seed!(20191202)
    (data, obj, mask) = testdata_brainweb_9comp(:stfrblochsim, Δω = Δω, κ = κ, σ = σ)
    Nnoise = count(getmask(obj, :background))
    noise = testdata_noise(N = Nnoise, σ = σ)
    y = reduce(vcat, (transpose(data[:,:,d][mask]) for d = 1:11))
    t = @elapsed ffhat = stfrperk(y, noise, Δω = Δω[mask], κ = κ[mask])[1][2,:]
    (ffmap = zeros(size(obj)))[mask] = ffhat
    (ff = obj.ff)[.!mask] .= 0

    rmse_wm = rmse(:wm, obj, ffmap, ff)
    rmse_gm = rmse(:gm, obj, ffmap, ff)

    mean_wm = avg(:wm, obj, ffmap)
    mean_gm = avg(:gm, obj, ffmap)

    std_wm = stdev(:wm, obj, ffmap)
    std_gm = stdev(:gm, obj, ffmap)

    @save(modulepath("estimation/results/brainweb_9comp_STFR2-PERK.jld"),
        ff, ffmap, rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t)

    # STFR3-PERK
    Random.seed!(20191202)
    (data, obj, mask) = testdata_brainweb_9comp(:stfrblochsim, Δω = Δω, κ = κ, σ = σ)
    Nnoise = count(getmask(obj, :background))
    noise = testdata_noise(N = Nnoise, σ = σ)
    y = reduce(vcat, (transpose(data[:,:,d][mask]) for d = 1:11))
    t = @elapsed ffhat = stfrblochsimperk(y, noise, Δω = Δω[mask], κ = κ[mask])[1][2,:]
    (ffmap = zeros(size(obj)))[mask] = ffhat
    (ff = obj.ff)[.!mask] .= 0

    rmse_wm = rmse(:wm, obj, ffmap, ff)
    rmse_gm = rmse(:gm, obj, ffmap, ff)

    mean_wm = avg(:wm, obj, ffmap)
    mean_gm = avg(:gm, obj, ffmap)

    std_wm = stdev(:wm, obj, ffmap)
    std_gm = stdev(:gm, obj, ffmap)

    @save(modulepath("estimation/results/brainweb_9comp_STFR3-PERK.jld"),
        ff, ffmap, rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t)

    # STFR3-PERK-JE
    Random.seed!(20191202)
    (data, obj, mask) = testdata_brainweb_9comp(:stfrblochsim, Δω = Δω, κ = κ, σ = σ)
    Nnoise = count(getmask(obj, :background))
    noise = testdata_noise(N = Nnoise, σ = σ)
    y = reduce(vcat, (transpose(data[:,:,d][mask]) for d = 1:11))
    t = @elapsed ffhat = stfrblochsimperk(y, noise)[1][2,:]
    (ffmap = zeros(size(obj)))[mask] = ffhat
    (ff = obj.ff)[.!mask] .= 0

    rmse_wm = rmse(:wm, obj, ffmap, ff)
    rmse_gm = rmse(:gm, obj, ffmap, ff)

    mean_wm = avg(:wm, obj, ffmap)
    mean_gm = avg(:gm, obj, ffmap)

    std_wm = stdev(:wm, obj, ffmap)
    std_gm = stdev(:gm, obj, ffmap)

    @save(modulepath("estimation/results/brainweb_9comp_STFR3-PERK-JE.jld"),
        ff, ffmap, rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t)

    # MESE-PERK
    Random.seed!(20191202)
    (data, obj, mask) = testdata_brainweb_9comp(:mese, Δω = Δω, κ = κ, σ = σ)
    Nnoise = count(getmask(obj, :background))
    noise = testdata_noise(N = Nnoise, σ = σ)
    y = reduce(vcat, (transpose(data[:,:,d][mask]) for d = 1:32))
    t = @elapsed ffhat = meseblochsimperk(y, noise)[1][2,:]
    (ffmap = zeros(size(obj)))[mask] = ffhat
    (ff = obj.ff)[.!mask] .= 0

    rmse_wm = rmse(:wm, obj, ffmap, ff)
    rmse_gm = rmse(:gm, obj, ffmap, ff)

    mean_wm = avg(:wm, obj, ffmap)
    mean_gm = avg(:gm, obj, ffmap)

    std_wm = stdev(:wm, obj, ffmap)
    std_gm = stdev(:gm, obj, ffmap)

    @save(modulepath("estimation/results/brainweb_9comp_MESE-PERK.jld"),
        ff, ffmap, rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t)

    # MESE-NNLS
    Random.seed!(20191202)
    (data, obj, mask) = testdata_brainweb_9comp(:mese, Δω = Δω, κ = κ, σ = σ)
    y = reduce(vcat, (transpose(data[:,:,d][mask]) for d = 1:32))
    β = 0.0151 * mean(y[1,:])
    y = [y[:,n] for n = 1:count(mask)]
    t = @elapsed (ffhat,) = mesennls(y, β = β)
    (ffmap = zeros(size(obj)))[mask] = ffhat
    (ff = obj.ff)[.!mask] .= 0

    rmse_wm = rmse(:wm, obj, ffmap, ff)
    rmse_gm = rmse(:gm, obj, ffmap, ff)

    mean_wm = avg(:wm, obj, ffmap)
    mean_gm = avg(:gm, obj, ffmap)

    std_wm = stdev(:wm, obj, ffmap)
    std_gm = stdev(:gm, obj, ffmap)

    @save(modulepath("estimation/results/brainweb_9comp_MESE-NNLS.jld"),
        ff, ffmap, rmse_wm, rmse_gm, mean_wm, mean_gm, std_wm, std_gm, t)

end
