"""
    createdesignA_bias()

Create file from which to load scan design A.
"""
function createdesignA_bias()
    
    # Check to make sure the data/ file exists and create it if not
    isdir(modulepath("bias/data")) || mkdir(modulepath("bias/data"))

    # Create design A if it does not already exist
    if !isfile(modulepath("bias/data/designA.jld"))

        Tfree = [10.3, 10.3, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]
        Tg = [2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8]
        α = [0.0872665, 0.0872665, 0.261799, 0.261799, 0.261637, 0.261799,
             0.261799, 0.261799, 0.261799, 0.199619, 0.261799]
        β = [0.0, 0.0, 0.261799, 0.261799, 0.202236, 0.261799, 0.232563,
             0.261799, 0.260596, 0.00439381, 0.25111]
        ϕ = [0.0, 0.0, -2.43209, -1.88717, -1.15143, -0.488278, 0.452156,
             1.12332, 1.81616, 2.55347, 3.02002]
        TE = [4.0, 6.3, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
        
        P = [[Tfree, Tg, α, β, ϕ, TE]]
        @save(modulepath("bias/data/designA.jld"), P)

    end

end

