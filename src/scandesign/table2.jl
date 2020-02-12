"""
    table2()

Display the data used to create Table 2. Table 2 reports the optimized scan
parameters for two scan designs: design A that accounts for the additional
myelin water off-resonance frequency Δωf, and design B that ignores Δωf.
"""
function table2()

    println("Design A")
    displaydesignA()
    println("Design B")
    displaydesignB()

end

"""
    displaydesignA()

Display the scan parameters for design A.
"""
function displaydesignA()

    # Make sure scan design file exists
    createdesignA_scandesign()

    displaydesign(modulepath("scandesign/results/designA.jld"))

end

"""
    displaydesignB()

Display the scan parameters for design B.
"""
function displaydesignB()

    # Make sure scan design file exists
    createdesignB_scandesign()

    displaydesign(modulepath("scandesign/results/designB.jld"))

end

"""
    displaydesign(filename)

Display the scan design found with the given file name.
"""
function displaydesign(filename)

    P = load(filename, "P")
    Psort = sortdesign(P)
    Pannotated = ["Scan #" transpose(1:size(Psort, 2)); ["α", "β", "ϕ"] Psort]
    display(Pannotated)

end

"""
    sortdesign(P)

Sort the given scan design by putting the fixed SPGR scans first and then
sorting the STFR scans by tip-up phase ϕ.
"""
function sortdesign(P)

    if length(P) == 3
        Pspgr = [rad2deg(P[1][1][1]) rad2deg(P[2][1][1]); 0 0; 0 0]
    else
        Pspgr = []
    end

    Pstfr = reduce(hcat, sort([(rad2deg.([P[end][1] P[end][2] P[end][3]]))[r,:]
        for r = 1:9], by = x -> x[3]))

    return round.([Pspgr Pstfr], digits = 1)

end
