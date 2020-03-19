"""
    figure6table4()

Create Figure 6 and Table 4. Figure 6 and Table 4 compare various MWF estimation
techniques using in vivo data from a healthy volunteer.
"""
function figure6table4()

    # Run figure6table4data() if it hasn't run yet
    isfile(modulepath("estimation/results/invivo.jld")) || figure6table4data()

    (ffmapstfr1, ffmapstfr2, ffmapmese1, ffmapmese2, tstfr1, tstfr2, tmese1, tmese2) =
        load(modulepath("estimation/results/invivo.jld"), "ffmapstfr1", "ffmapstfr2",
        "ffmapmese1", "ffmapmese2", "tstfr1", "tstfr2", "tmese1", "tmese2")
    ff = [ffmapstfr1 ffmapmese1; ffmapstfr2 ffmapmese2]
    p = heatmap(ff', color = :fire, aspect_ratio = :equal, yflip = true,
        xticks = [], yticks = [], clims = (0, 0.3),
        dpi = 300, size = (600, 600))
    display(p)

    # Make sure data exists and download it if not
    isdir(modulepath("estimation/data/invivo")) || mkpath(modulepath("estimation/data/invivo"))
    isfile(modulepath("estimation/data/invivo/rois.mat")) ||
        getdata("https://doi.org/10.7302/9t0w-2970", modulepath("estimation/data/invivo/rois.mat"))

    masks = matread(modulepath("estimation/data/invivo/rois.mat"))
    stats = x -> begin
        avg = [round(mean(x[m]), digits = 3) for (key, m) in masks if key != "mic"]
        stdev = [round(std(x[m]), digits = 3) for (key, m) in masks if key != "mic"]
        return (avg, stdev)
    end

    (meanstfr1, stdstfr1) = stats(ffmapstfr1)
    (meanstfr2, stdstfr2) = stats(ffmapstfr2)
    (meanmese1, stdmese1) = stats(ffmapmese1)
    (meanmese2, stdmese2) = stats(ffmapmese2)

    tstfr1 = round(tstfr1, digits = 1)
    tstfr2 = round(tstfr2, digits = 1)
    tmese1 = round(tmese1, digits = 1)
    tmese2 = round(tmese2, digits = 1)

    todisplay = ["ROI" "STFR2-PERK Mean" "STFR2-PERK St. Dev." "STFR3-PERK Mean" "STFR3-PERK St. Dev." "MESE-NNLS Mean" "MESE-NNLS St. Dev." "MESE-PERK Mean" "MESE-PERK St. Dev.";
        [k for k in keys(masks) if k != "mic"] meanstfr1 stdstfr1 meanstfr2 stdstfr2 meanmese1 stdmese1 meanmese2 stdmese2]

    display(todisplay[:,[1,2,3,4,5]])
    display(todisplay[:,[1,6,7,8,9]])

    println("STFR2-PERK: $tstfr1 s")
    println("STFR3-PERK: $tstfr2 s")
    println("MESE-NNLS:  $tmese1 s")
    println("MESE-PERK:  $tmese2 s")

end

"""
    figure6table4data()

Create data for Figure 6 and Table 4.
"""
function figure6table4data()

    # Make sure scan design file exists and create it if not
    createdesignAsorted_estimation()

    # Make sure data exists and download it if not
    isdir(modulepath("estimation/data/invivo")) || mkpath(modulepath("estimation/data/invivo"))
    isfile(modulepath("estimation/data/invivo/recon.jld")) ||
        getdata("https://doi.org/10.7302/y5ye-2706", modulepath("estimation/data/invivo/recon.jld"))
    isfile(modulepath("estimation/data/invivo/b1slice5.mat")) ||
        getdata("https://doi.org/10.7302/epwx-vb86", modulepath("estimation/data/invivo/b1slice5.mat"))
    isfile(modulepath("estimation/data/invivo/meseslice5.mat")) ||
        getdata("https://doi.org/10.7302/vat5-gw24", modulepath("estimation/data/invivo/meseslice5.mat"))
    isfile(modulepath("estimation/data/invivo/rois.mat")) ||
        getdata("https://doi.org/10.7302/9t0w-2970", modulepath("estimation/data/invivo/rois.mat"))

    datapath = modulepath("estimation/data/invivo")

    (img, mask) = load(joinpath(datapath, "recon.jld"), "img", "mask")
    b1 = matread(joinpath(datapath, "b1slice5.mat"))["map"][301:500,:]
    mese = matread(joinpath(datapath, "meseslice5.mat"))["mese"][301:500,:,:]

    # Noise masks (found through trial and error)
    nmaskx = [1:25 176:200]
    nmasky = 10:190

    # Convert data into inputs to stfrperk() and mesennls() and meseperk()
    ystfr = transpose(reduce(hcat, [img[:,:,d][mask] for d = 1:size(img,3)]))
    κ = b1[mask]
    nstfr = vec(reduce(vcat, [img[:,:,d][nmaskx,nmasky] for d = 1:size(img,3)]))
    ymese = transpose(reduce(hcat, (mese[:,:,e][mask] for e = 1:32))) # [D,N]
    nmese = vec(reduce(vcat, (mese[:,:,e][nmaskx,nmasky] for e = 1:32)))
    β = 0.0151 * mean(ymese[1,:])

    # Compute SNR
    masks = matread(joinpath(datapath, "rois.mat"))
    maskroi = reduce((x, y) -> x .| y, (m for (key, m) in masks if key != "mgm"))
    nstfrstd = [std(img[:,:,d][nmaskx,nmasky]) for d = 1:11] ./ sqrt(2 - π/2)
    nmesestd = [std(mese[:,:,d][nmaskx,nmasky]) for d = 1:32] ./ sqrt(2 - π/2)
    snrstfr = [mean(img[:,:,d][maskroi]) / nstfrstd[d] for d = 1:11]
    snrmese = [mean(mese[:,:,d][maskroi]) / nmesestd[d] for d = 1:32]

    # Run PERK and NNLS
    @time tstfr1 = @elapsed ffstfr1 = stfrperk(ystfr, nstfr, κ = κ)[1][2,:]
    @time tstfr2 = @elapsed ffstfr2 = stfrblochsimperk(ystfr, nstfr, κ = κ)[1][2,:]
    @time tmese2 = @elapsed ffmese2 = meseblochsimperk(ymese, nmese)[1][2,:]
    @time tmese1 = @elapsed (ffmese1,) = mesennls([ymese[:,n] for n = 1:count(mask)], β = β)

    (ffmapstfr1 = zeros(size(mask)))[mask] = ffstfr1
    (ffmapstfr2 = zeros(size(mask)))[mask] = ffstfr2
    (ffmapmese1 = zeros(size(mask)))[mask] = ffmese1
    (ffmapmese2 = zeros(size(mask)))[mask] = ffmese2

    # Make sure results/ folder exists
    isdir(modulepath("estimation/results")) || mkdir(modulepath("estimation/results"))

    @save(modulepath("estimation/results/invivo.jld"), ffmapstfr1, ffmapstfr2,
        ffmapmese1, ffmapmese2, tstfr1, tstfr2, tmese1, tmese2, snrstfr, snrmese)

end
