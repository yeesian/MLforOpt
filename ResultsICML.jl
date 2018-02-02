# Results for ICML

using JLD, OPFRecourse

include("WindowSize.jl")


# Set parameters
alpha = 0.1
epsilon = 0.05
delta = 0.01
gamma = 2


#M = Dict{String,Any}()
#K_M = Dict{String,Any}()
#W = Dict{String,Any}()
#R_WM = Dict{String,Any}()


k = 0 # clumsy indexing...

# Running over OPF Benchmark Cases
for f in [
        "case3_lmbd",
        "case5_pjm",
        "case14_ieee",
        "case24_ieee_rts",
        #"case30_as",
        #"case30_fsr",
        "case30_ieee",
        "case39_epri",
        "case57_ieee",
        "case73_ieee_rts",
        # "case89_pegase", # (infeasible)
        "case118_ieee",
        "case162_ieee_dtc",
        "case200_pserc",
        "case240_pserc",
        "case300_ieee",
        # "case1354_pegase", # (infeasible)
         "case1888_rte",
         "case1951_rte",
         # "case2383wp_k" # (infeasible)
         #"case2736sp_k",
         "case2737sop_k", #seems fine
         #"case2746wop_k", # (infeasible)
         #"case2746wp_k", # (infeasible)
         "case2848_rte", #seems fine
         "case2868_rte", #seems fine
         "case2869_pegase", #seems fine (but seems to take longer)
         "case2869_pegase", #seems fine (but also longer)
         #"case3012wp_k", # (infeasible)
         #"case3120sp_k", # (infeasible)
         #"case3375wp_k", # (singluar)
         "case6468_rte", # seems fine (abou 1 min 30s for 10 samples)
         "case6495_rte", # seems fine (abou 1 min 37s for 10 samples)
         "case6515_rte", # seems fine (abou 1 min 40s for 10 samples)
         #"case9241_pegase", # infeasible
         "case13659_pegase"
    ]

    # Clumsy indexing....
    k=k+1

    print("Working on $f: ")
    data_file = string("darwin/results/", f, "3.jld")

    (M, W, R_MW, K_M) = RunLearningAlgorithm(alpha, delta, epsilon, gamma, data_file)

    println()
    print(["M: ", M, "W: ", W, "R_MW: ", R_MW, "K_M: ", K_M])
    println()

end
