"""
    createdesignA_scandesign()

Create file from which to load scan design A.
"""
function createdesignA_scandesign()

    # Check to make sure the results/ file exists and create it if not
    isdir(modulepath("scandesign/results")) || mkdir(modulepath("scandesign/results"))

    # Create design A if it does not already exist
    if !isfile(modulepath("scandesign/results/designA.jld"))

        P = [
             [[0.0872665]],
             [[0.0872665]],
             [
              [0.261799, 0.261799, 0.261799, 0.261637, 0.261799, 0.199619,
               0.261799, 0.261799, 0.261799],
              [0.261799, 0.261799, 0.232563, 0.202236, 0.261799, 0.00439381,
               0.25111, 0.261799, 0.260596],
              [1.12332, -0.488278, 0.452156, -1.15143, -2.43209, 2.55347,
               3.02002, -1.88717, 1.81616]
             ]
            ]
        cost = 0.00177008
        flag = :MAXTIME_REACHED
        description = "ΔΔω known from -50 to 50 Hz, nsamp = 10, 20 hours"

        @save(modulepath("scandesign/results/designA.jld"), P, cost, flag,
              description)

    end

end

