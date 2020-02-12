import Pkg

# Add packages that are not registered in Julia's general registry.
# This is necessary to avoid errors when adding the STFRMWF package.
Pkg.add(Pkg.PackageSpec(url = "https://github.com/rdeits/NNLS.jl", rev = "5833899d31c840a318de6395b0e400da507e00c1"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/StevenWhitaker/STFR.jl", rev = "v0.0.2"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/StevenWhitaker/ScanDesign.jl", rev = "v0.0.2"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/StevenWhitaker/PERK.jl", rev = "v0.0.4"))

# Add PyPlot
Pkg.add("PyPlot")
