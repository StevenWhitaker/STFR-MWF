"""
    figureS4()

Create Supporting Information Figure S4. Figure S4 shows the 11 STFR scans and
the separately estimated B0 map.
"""
function figureS4()

    # Make sure data exists and download it if not
    isdir(modulepath("estimation/data/invivo")) || mkpath(modulepath("estimation/data/invivo"))
    isfile(modulepath("estimation/data/invivo/recon.jld")) ||
        getdata("https://doi.org/10.7302/y5ye-2706", modulepath("estimation/data/invivo/recon.jld"))
    isfile(modulepath("estimation/data/invivo/headmask.mat")) ||
        getdata("https://doi.org/10.7302/xdrn-gh31", modulepath("estimation/data/invivo/headmask.mat"))

    datapath = modulepath("estimation/data/invivo")

    img = load(joinpath(datapath, "recon.jld"), "img")

    img = [img[:,:,1] img[:,:,5] img[:,:,9];
        img[:,:,2] img[:,:,6] img[:,:,10];
        img[:,:,3] img[:,:,7] img[:,:,11];
        img[:,:,4] img[:,:,8] zeros(size(img[:,:,1]))]

    p = heatmap(img', color = :grays, aspect_ratio = :equal, yflip = true,
        xticks = [], yticks = [], clims = (0, 0.8), cbar = :none,
        dpi = 300, size = (1200, 900))
    display(p)

    b0 = load(joinpath(datapath, "recon.jld"), "b0map")
    mask = matread(joinpath(datapath, "headmask.mat"))["mbrain"]
    b0[.!mask] .= 0
    pb0 = heatmap(b0', color = :coolwarm, aspect_ratio = :equal,
        xticks = [], yticks = [], clims = (-80, 80), yflip = true,
        dpi = 300, size = (200, 200))
    display(pb0)

end
