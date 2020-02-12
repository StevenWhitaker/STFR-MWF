"""
    createdesignB_scandesign()

Create file from which to load scan design B.
"""
function createdesignB_scandesign()

    # Check to make sure the results/ file exists and create it if not
    isdir(modulepath("scandesign/results")) || mkdir(modulepath("scandesign/results"))

    # Create design B if it does not already exist
    if !isfile(modulepath("scandesign/results/designB.jld"))

        P = [
             [[0.0872665]],
             [[0.0872665]],
             [
              [0.261799, 0.261799, 0.261799, 0.261799, 0.261799, 0.261799,
               0.261799, 0.261799, 0.261799],
              [0.258688, 0.261799, 0.258688, 0.259702, 0.253355, 0.259702,
               0.261799, 0.253355, 0.0],
              [-0.249261, -2.43067, 0.249262, -1.11159, -1.97728, 1.11159,
               2.43067, 1.97728, 1.45159]
             ]
            ]
        cost = 0.0366846
        flag = :MAXTIME_REACHED
        description = "ΔΔω known from -50 to 50 Hz, nsamp = 10, 20 hours, " *
                      "Δωf = 0 known"

        @save(modulepath("scandesign/results/designB.jld"), P, cost, flag,
              description)

    end

end

