using JLD, OPFRecourse

include("WindowSize.jl")
include("LearningAlgorithms_DCOPF.jl")


f = "case24_ieee_rts"


Mmin = 1

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
# Evaluate stopping criterion

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

data_file = string("darwin/results_uniform/", f, "3.jld")
scenarios = JLD.load(data_file, "scenarios")

(M[k], W[k], R_MW[k], K_M[k], R_OS[k], W_test[k]) = RunStreamingAlgorithm(alpha, delta, epsilon, gamma, scenarios, testsize, Mmin)
println(" M: ", M[k], " W: ", W[k], " R_MW: ", R_MW[k], " K_M: ", K_M[k], " R_OS: ", R_OS[k], " Test size: ", W_test[k])


end
