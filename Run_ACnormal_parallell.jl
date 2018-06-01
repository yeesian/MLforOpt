
using Ipopt
using Distributions
using DataStructures
using JuMP
using PowerModels
using ProgressMeter
using JLD



# Running over OPF Benchmark Cases
f = ARGS[1]

print("Working on $f: ")
data_file = string("pglib-opf/pglib_opf_", f, ".m")
#data_file = string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/pglib_opf_", f, ".m")

include("find_active_set.jl")
include("basis_policy.jl")
include("WindowSize.jl")
include("LearningAlgorithms_ACOPF.jl")

NLsolver = IpoptSolver(print_level=0)

alpha = 0.1
delta = 0.01
epsilon = 0.05
gamma = 2
Minitial = 1

M, W, RoD, K_M, results = RunStreamingAlgorithmAC(alpha, delta, epsilon, gamma, Minitial, data_file, NLsolver, tol = 1e-5)

JLD.save("results_ACnormal_180531/$(f).jld", "results", results, "M", M, "W", W, "RoD", RoD, "K_M", K_M)

println()
println("!!!!!!!!!!!!!!! FINISHED WORKING ON $f !!!!!!!!!!!!!!!!!!!!")
println()

#for i=3 #in 1:5
#     @time ref = OPFRecourse.NetworkReference(data_file, σscaling=0.01*i);
#     m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver());
#
#     # Generate new uncertainties based on uniform distribution
#     srand(1234)
#     nsamples = 50_000
#     ωsamples = zeros(ref.nbus, nsamples)
#     nonzeroindices = (1:length(ref.stdω))[ref.stdω .> 1e-5]
#     for b in nonzeroindices
#         ωsamples[b,:] = rand(Uniform(-3*ref.stdω[b], 3*ref.stdω[b]), nsamples)
#     end
#     # solve OPF for uniform distribution realizations
#     scenarios = OPFRecourse.OPFScenarios(ref, m, ωsamples);
#     JLD.save("results_uniform/$(f)$(i).jld", "scenarios", scenarios)
#     print("$i ")
# #end
