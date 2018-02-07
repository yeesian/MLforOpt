# Results for ICML

using JLD, OPFRecourse

include("WindowSize.jl")


# Set parameters
alpha = 0.05
epsilon = 0.04
delta = 0.01
gamma = 2

# Out of sample test size
testsize = 20_000

M = Vector{Float64}(15)
K_M = Vector{Float64}(15)
W = Vector{Float64}(15)
R_MW = Vector{Float64}(15)
R_OS = Vector{Float64}(15)
W_test = Vector{Float64}(15)
#K_M = Dict{String,Any}()
#W = Dict{String,Any}()
#R_WM = Dict{String,Any}()


# Running over OPF Benchmark Cases
for (k,f) in enumerate([
        "case3_lmbd",
        "case5_pjm",
        "case14_ieee",
        "case24_ieee_rts",
        #"case30_as",
        #"case30_fsr",
        "case30_ieee",
        "case39_epri",
        "case57_ieee",
        "case73_ieee_rts", # will not terminate for uniform
        # "case89_pegase", # (infeasible)
        "case118_ieee",
        "case162_ieee_dtc",
        "case200_pserc", # uniform = 6000-7000 samples
        "case240_pserc", # will not terminate for uniform
        "case300_ieee", # uniform = 9000-10000 smples
        # "case1354_pegase", # (infeasible)
         "case1888_rte", # less than 100 samples
         "case1951_rte" # less than 100 samples
])
         # "case2383wp_k" # (infeasible)
         #"case2736sp_k",
#         "case2737sop_k", #seems fine
         #"case2746wop_k", # (infeasible)
         #"case2746wp_k", # (infeasible)
#         "case2848_rte", #seems fine
#         "case2868_rte", #seems fine
#         "case2869_pegase", #seems fine (but seems to take longer)
#         "case2869_pegase", #seems fine (but also longer)
         #"case3012wp_k", # (infeasible)
         #"case3120sp_k", # (infeasible)
         #"case3375wp_k", # (singluar)
#         "case6468_rte", # seems fine (abou 1 min 30s for 10 samples)
#         "case6495_rte", # seems fine (abou 1 min 37s for 10 samples)
#         "case6515_rte", # seems fine (abou 1 min 40s for 10 samples)
         #"case9241_pegase", # infeasible
         #"case13659_pegase"
#])

    # Clumsy indexing....

    print("Working on $f: ")
    #data_file = string("darwin/results_gaussian/", f, "3.jld")
    data_file = string("darwin/results_uniform/", f, "3.jld")

    scenarios = JLD.load(data_file, "scenarios")

    Mmin = 1

    if f=="case73_ieee_rts"
        Mmin = 23000
    end
    if f=="case200_pserc"
        Mmin = 6000
    end
    if f=="case240_pserc"
        Mmin = 23000
    end
    if f=="case300_ieee"
        Mmin = 9000
    end

    # Run the learning algorithm
    (M[k], W[k], R_MW[k], K_M[k], R_OS[k], W_test[k]) = RunLearningAlgorithm(alpha, delta, epsilon, gamma, scenarios, testsize, Mmin)

    # Compute probability of each basis
    #basisM.whichscenario = scenarios.whichscenario[1+offset:m+offset,:]
    #prob_M = [length(basisM.whichscenario[(uniqueM[b,1],uniquesM[b,2])]) for b in 1:size(uniquescenarios,1)]/M
    #basisW.whichscenario = scenarios.whichscenario[m+1+offset:m+w+offset,:]
    #prob_W = [length(basisW.whichscenario[(uniqueM[b,1],uniquesM[b,2])]) for b in 1:size(uniquescenarios,1)]/W




    println()
    print(" M: ", M[k], " W: ", W[k], " R_MW: ", R_MW[k], " K_M: ", K_M[k], " R_OS: ", R_OS[k], " Test size: ", W_test[k])
    println()

end

JLD.save("summary_uniform.jld", "M", M, "K_M", K_M, "R_MW", R_MW, "W", W, "R_OS", R_OS, "W_test", W_test)
#JLD.save("results_uniform/summary_gaussian.jld", "M", M, "K_M", K_M, "R_MW", R_MW, "W", W)
