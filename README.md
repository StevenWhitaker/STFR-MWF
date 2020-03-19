# STFR-MWF
This repository contains code
for reproducing the results in the paper:

Steven T. Whitaker, Gopal Nataraj, Jon-Fredrik Nielsen, and Jeffrey A. Fessler \
Myelin Water Fraction Estimation Using Small-Tip Fast Recovery MRI \
Magnetic Resonance in Medicine, 2020, Accepted

## Getting Started
The code in this repository is structured as a Julia module.
It was coded using Julia 1.1.1,
though other versions may work as well.

1. Download or clone this repository.
2. Download the BrainWeb dataset [here](https://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?alias=phantom_1.0mm_normal_crisp&download=1).
   - Choose raw byte (unsigned) for the file format.
   - Choose none for compression.
   - Place the downloaded file in `<path_to_this_repo>/STFR-MWF/src/estimation/data/BrainWeb/`,
     where `<path_to_this_repo>` is where you downloaded this repo on your computer.
     (You will have to create the `data` and `BrainWeb` directories.)
3. Download Julia 1.1.1 [here](https://julialang.org/downloads/oldreleases/).
4. Run Julia.
5. Change directories to this repo with
   ```julia
   julia> cd("<path_to_this_repo>")
   ```
6. Install some necessary packages and load the code by running `setup.jl` via
   ```julia
   julia> include("setup.jl")
   ```
7. Run any function with
    ```julia
    julia> STFRMWF.func() # Replace func with actual function name
    ```

Steps 1 through 3 only need to be done once,
and steps 4 through 6 only need to be done once each time you start Julia.
(Note that running `setup.jl` again will not reinstall packages,
so will run faster after the first time.)

## Reproducing Results
To reproduce the results in the paper, call the corresponding function.
For example, to reproduce Figure 3 run
```julia
julia> STFRMWF.figure3()
```
Some results are used for both a figure and a table.
For example, Figure 4 and Table 3 use the same data,
so to reproduce both run
```julia
julia> STFRMWF.figure4table3()
```

The first time calling these functions may take a while
(especially the biased CRLB and NNLS results),
but after running once the results will be saved,
so future calls will simply load the results and display them.

## Data
The data used in the paper and in this code is available at
[Deep Blue Data](https://doi.org/10.7302/nw6e-1d66).
