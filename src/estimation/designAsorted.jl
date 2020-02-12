"""
    createdesignAsorted_estimation()

Create file from which to load scan design A that was implemented on the
scanner.
"""
function createdesignAsorted_estimation()
    
    # Check to make sure the data/ file exists and create it if not
    isdir(modulepath("estimation/data")) || mkdir(modulepath("estimation/data"))

    # Create design A if it does not already exist
    if !isfile(modulepath("estimation/data/designAsorted.jld"))

        Tfree = [10.156, 10.156, 7.856, 7.856, 7.856, 7.856, 7.856, 7.856,
                 7.856, 7.856, 7.856] 
        Tg = [3.144, 3.144, 3.144, 3.144, 3.144, 3.144, 3.144, 3.144, 3.144,
              3.144, 3.144]
        α = [0.0698132, 0.0698132, 0.261799, 0.261799, 0.261637, 0.261799,
             0.261799, 0.261799, 0.261799, 0.199619, 0.261799]
        β = [0.0, 0.0, 0.261799, 0.261799, 0.202236, 0.261799, 0.232563,
             0.261799, 0.260596, 0.00439381, 0.25111]
        ϕ = [0.0, 0.0, -2.43209, -1.88717, -1.15143, -0.488278, 0.452156,
             1.12332, 1.81616, 2.55347, 3.02002]
        TE = [3.928, 6.228, 3.928, 3.928, 3.928, 3.928, 3.928, 3.928, 3.928,
              3.928, 3.928]
        
        P = [[Tfree, Tg, α, β, ϕ, TE]]
        @save(modulepath("estimation/data/designAsorted.jld"), P)

    end

end

