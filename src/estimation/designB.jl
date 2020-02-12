"""
    createdesignB_estimation()

Create file from which to load scan design B.
"""
function createdesignB_estimation()
    
    # Check to make sure the data/ file exists and create it if not
    isdir(modulepath("estimation/data")) || mkdir(modulepath("estimation/data"))

    # Create design A if it does not already exist
    if !isfile(modulepath("estimation/data/designB.jld"))

        Tfree = [10.3, 10.3, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]
        Tg = [2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8]
        α = [0.0872665, 0.0872665, 0.261799, 0.261799, 0.261799, 0.261799,
             0.261799, 0.261799, 0.261799, 0.261799, 0.261799]
        β = [0.0, 0.0, 0.261799, 0.253355, 0.259702, 0.258688, 0.258688,
             0.259702, 0.0, 0.253355, 0.261799]
        ϕ = [0.0, 0.0, -2.43067, -1.97728, -1.11159, -0.249261, 0.249262,
             1.11159, 1.45159, 1.97728, 2.43067]
        TE = [4.0, 6.3, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
        
        P = [[Tfree, Tg, α, β, ϕ, TE]]
        @save(modulepath("estimation/data/designB.jld"), P)

    end

end

