using Ipopt
using Distributions
using DataStructures
using JuMP
using PowerModels
using ProgressMeter
using JLD


include("find_active_set.jl")
include("basis_policy.jl")
include("WindowSize.jl")
include("LearningActiveSet_ACOPF.jl")
include("LearningActiveConstraints_ACOPF.jl")
include("Learning_AllRoD_ACOPF.jl")

NLsolver = IpoptSolver(print_level=0)

tol = 1e-5

#filename = "data/nesta_case1397sp_eir.m"
#filename = "pglib-opf/pglib_opf_case5_pjm.m"
#filename = "pglib-opf/pglib_opf_case24_ieee_rts.m"
#filename = "pglib-opf/pglib_opf_case118_ieee.m"
#filename = "pglib-opf/pglib_opf_case57_ieee.m"
#filename = "pglib-opf/pglib_opf_case240_pserc.m"
filename = "pglib-opf/pglib_opf_case300_ieee.m"

network_data = PowerModels.parse_file(filename)
ref = PowerModels.build_ref(network_data)[:nw][0]

#result = run_ac_opf(filename,NLsolver)

m = Model(solver = NLsolver)

#jm, const_refs, var_refs = post_ac_opf_withref(network_data,m)
#jm, const_refs, var_refs, nl_refs = post_ac_opf_withref_uncertainty(network_data,m)

#status = solve(jm)


alpha = 0.9
delta = 0.1
epsilon = 0.4
gamma = 2
Minitial = 1

M, W, RoD, K_M, results = RunStreamingAlgorithmAC_AllRoD(alpha, delta, epsilon, gamma, Minitial, filename, NLsolver)
#M, W, RoD, K_M, results = RunStreamingAlgorithmAC(alpha, delta, epsilon, gamma, Minitial, filename)

#JLD.save("$(filename).jld", "results", results, "M", M, "W", W, "RoD", RoD, "K_M", K_M)

JLD.save("testing.jld", "results", results, "M", M, "W", W, "RoD", RoD, "K_M", K_M)





#
# network_data = PowerModels.parse_file(filename)
# ref = PowerModels.build_ref(network_data)[:nw][0]
# nonzeroindices = [i for (i,loads) in ref[:bus_loads] if length(loads) > 0]
#
# m = Model(solver = NLsolver)
#
# #jm, const_refs, var_refs = post_ac_opf_withref(network_data,m)
# jm, const_refs, var_refs, nl_refs = post_ac_opf_withref_uncertainty(network_data,m)
#
# status = solve(jm)
#
# # Constructing distribution for samples
# sigma = 0.1
# load = [sum(ref[:load][l]["pd"] for l in ref[:bus_loads][i]) for i in nonzeroindices]
# w = Distributions.MvNormal(
#     zeros(length(nonzeroindices)),
#     diagm((sigma*load).^2)
# )
#
# # Run OPF to get active sets
# active_rows = Set{Vector{Int}}()
# active_cols_upper = Set{Vector{Int}}()
# active_cols_lower = Set{Vector{Int}}()
#
# dict_active_rows = Dict{Vector{Int}, Int}()
# dict_active_cols_upper = Dict{Vector{Int}, Int}()
# dict_active_cols_lower = Dict{Vector{Int}, Int}()
#
# observed_active_sets = Set{Int}()
# all_active_sets = Set{Tuple{Int,Int,Int}}()
# active_sets = Tuple{Int,Int,Int}[]
# active_index = Dict{Tuple{Int,Int,Int}, Int}()
#
# observedsets = Set{Int}()
# windowcount = Dict{Int,Int}() # active set index => # of times it appears in the window
# q = DataStructures.Queue(Int)
#
# nsamples = 100
# Minitial = 1
# W = 20
# threshold = 0.001
#
# # Initialization
#
# for i = 1:(Minitial+W)#:m+W
#     # 1. Generate sample
#     w_sample = rand(w,1)
#     # 2. Fix NL parameters
#     for (k,j) in enumerate(nonzeroindices)
#         setvalue(nl_refs["u"][j], w_sample[k])
#     end
#     # 3. Get active set
#     active_set = find_active_set(jm, const_refs, var_refs, tol)
#
#     # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)
#
#     # If a new active row or column is observed, add them to the set and create a new index
#     if !in(active_set["active_rows"], active_rows)
#         push!(active_rows, active_set["active_rows"])
#         dict_active_rows[active_set["active_rows"]] = length(active_rows)
#     end
#     if !in(active_set["active_cols_upper"], active_cols_upper)
#         push!(active_cols_upper, active_set["active_cols_upper"])
#         dict_active_cols_upper[active_set["active_cols_upper"]] = length(active_cols_upper)
#     end
#     if !in(active_set["active_cols_lower"], active_cols_lower)
#         push!(active_cols_lower, active_set["active_cols_lower"])
#         dict_active_cols_lower[active_set["active_cols_lower"]] = length(active_cols_lower)
#     end
#
#     # Put together the marker for the active set discovered now
#     new_active_set = (
#         dict_active_rows[active_set["active_rows"]],
#         dict_active_cols_upper[active_set["active_cols_upper"]],
#         dict_active_cols_lower[active_set["active_cols_lower"]]
#     )
#
#     # If the current active set has not been discovered
#     if !in(new_active_set, all_active_sets)
#         push!(all_active_sets, new_active_set)  # add the new active set to the (unordered) SET of discovered active sets
#         push!(active_sets, new_active_set)      # add the new active set to the (ordered) VECTORS of active sets
#         active_index[new_active_set] = length(active_sets)  # make a new index for the new active set
#     end
#
#     # Add the active set to the queue (the list of all active sets)
#     DataStructures.enqueue!(q, active_index[new_active_set])
#
#     # Update the count of how frequently the active sets appear in the queue
#     windowcount[active_index[new_active_set]] = get(windowcount, active_index[new_active_set], 0) + 1
# end
#
# println("before: $q")
# for i in 1:Minitial
#     curr_active_set = dequeue!(q)                   # remove the sample
#     push!(observed_active_sets, curr_active_set)    # put active set into observed active set
#     windowcount[curr_active_set] -= 1               # update windowcount
#     @assert windowcount[curr_active_set] >= 0
# end
# println("after: $q")
#
# # computing rate of discovery
# numerator = 0
# denominator = sum(values(windowcount))
# @assert denominator > 0
# for (as,v) in windowcount
#     if !in(as, observed_active_sets)
#         numerator += v
#     end
# end
# println("rate of discovery: $(numerator / denominator)")
#
# # Streaming
#
# M = Minitial;
# for j = 1:nsamples
#     M += 1;
#     Wold = W;
#     W += 1; # update to use WindowSize function which depends on M
#
#     println(W-Wold+1)
#
#     # add new samples to W
#     for i = 1:W-Wold+1 #add at least one sample:
#                        #replace the which will be moved in to M,
#                        #plus add new samples to cover an increase in the size of W
#         # 1. Generate sample
#         w_sample = rand(w,1)
#         # 2. Fix NL parameters
#         for (k,j) in enumerate(nonzeroindices)
#             setvalue(nl_refs["u"][j], w_sample[k])
#         end
#         # 3. Get active set
#         active_set = find_active_set(jm, const_refs, var_refs, tol)
#
#         # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)
#
#         # If a new active row or column is observed, add them to the set and create a new index
#         if !in(active_set["active_rows"], active_rows)
#             push!(active_rows, active_set["active_rows"])
#             dict_active_rows[active_set["active_rows"]] = length(active_rows)
#         end
#         if !in(active_set["active_cols_upper"], active_cols_upper)
#             push!(active_cols_upper, active_set["active_cols_upper"])
#             dict_active_cols_upper[active_set["active_cols_upper"]] = length(active_cols_upper)
#         end
#         if !in(active_set["active_cols_lower"], active_cols_lower)
#             push!(active_cols_lower, active_set["active_cols_lower"])
#             dict_active_cols_lower[active_set["active_cols_lower"]] = length(active_cols_lower)
#         end
#
#         # Put together the marker for the active set discovered now
#         new_active_set = (
#             dict_active_rows[active_set["active_rows"]],
#             dict_active_cols_upper[active_set["active_cols_upper"]],
#             dict_active_cols_lower[active_set["active_cols_lower"]]
#         )
#
#         # If the current active set has not been discovered
#         if !in(new_active_set, all_active_sets)
#             push!(all_active_sets, new_active_set)  # add the new active set to the (unordered) SET of discovered active sets
#             push!(active_sets, new_active_set)      # add the new active set to the (ordered) VECTORS of active sets
#             active_index[new_active_set] = length(active_sets)  # make a new index for the new active set
#         end
#
#         # Add the active set to the queue (the list of all active sets)
#         DataStructures.enqueue!(q, active_index[new_active_set])
#
#         # Update the count of how frequently the active sets appear in the queue
#         windowcount[active_index[new_active_set]] = get(windowcount, active_index[new_active_set], 0) + 1
#     end
#
#     # pushing one sample from W into M
#     for i in 1:Minitial
#         curr_active_set = dequeue!(q)                   # remove the sample
#         push!(observed_active_sets, curr_active_set)    # put active set into observed active set
#         windowcount[curr_active_set] -= 1               # update windowcount
#         @assert windowcount[curr_active_set] >= 0
#     end
#     println("Window $q at sample number $j")
#
#     # computing rate of discovery
#     numerator = 0
#     denominator = sum(values(windowcount))
#     @assert denominator > 0
#     for (as,v) in windowcount
#         if !in(as, observed_active_sets)
#             numerator += v
#         end
#     end
#     RoD = (numerator / denominator)
#     println("rate of discovery: $(numerator / denominator)")
#
#     if RoD <= threshold
#         break
#     end
#
# end
#
#
#
#
#
#
#
# # status = solve(jm)
#
# # active_set = find_active_set(jm, const_refs, var_refs, tol)
#
# # km = build_basis_policy_model(filename, NLsolver, active_set)
