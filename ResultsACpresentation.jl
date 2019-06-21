# Summarizing and plotting results for the AC OPF simulations

using JLD
using ProgressMeter

K_M = Vector{Int}(19)
M = Vector{Int}(19)
W = Vector{Int}(19)
RoD = Vector{Float64}(19)
RoD_Constraint = Vector{Float64}(19)
Rows = Vector{Int}(19)
Cols_upper = Vector{Int}(19)
Cols_lower = Vector{Int}(19)
Row_constraints = Vector{Int}(19)
Cols_upper_constraints = Vector{Int}(19)
Cols_lower_constraints = Vector{Int}(19)
All_Sets = Vector{Int}(19)


@showprogress 1 for (k,f) in enumerate([
        "results_AC_allRoD_180825/case3_lmbd",
        "results_AC_allRoD_180825/case5_pjm",
        "results_AC_allRoD_180825/case14_ieee",
        "results_AC_allRoD_180825/case24_ieee_rts",
        #"case30_as",
        #"case30_fsr",
        "results_AC_allRoD_180825/case30_ieee",
        "results_AC_allRoD_180825/case39_epri",
        "results_AC_allRoD_180825/case57_ieee",
        "results_AC_allRoD_180825/case73_ieee_rts", # will not terminate for uniform
        # "case89_pegase", # (infeasible)
        "results_AC_allRoD_180825/case118_ieee",
        "results_AC_allRoD_180825/case162_ieee_dtc",
        "results_AC_allRoD_180825/case200_pserc", # uniform = 6000-7000 samples
        "intermediate_results/pglib_opf_case240_pserc_iteration700", # will not terminate for uniform
        "intermediate_results/pglib_opf_case300_ieee_iteration1000", # uniform = 9000-10000 smples
        # "case1354_pegase", # (infeasible)
         "intermediate_results/pglib_opf_case1888_rte_iteration1100", # less than 100 samples
         "intermediate_results/pglib_opf_case1951_rte_iteration800", # less than 100 samples
         "intermediate_results/pglib_opf_case2737sop_k_iteration3600", # less than 100 samples
         "intermediate_results/pglib_opf_case2848_rte_iteration1000", # less than 100 samples
         "intermediate_results/pglib_opf_case2869_pegase_iteration700", # less than 100 samples
         "intermediate_results/pglib_opf_case6468_rte_iteration500" # less than 100 samples
])

    unique_rows = Vector{Int}()
    unique_cols_upper = Vector{Int}()
    unique_cols_lower = Vector{Int}()

    data_file = string(f, ".jld")
    data = JLD.load(data_file)

    K_M[k]  = data["K_M"]
    M[k]    = data["M"]
    W[k]    = data["W"]
    RoD[k]  = data["RoD"]
    Rows[k] = length(data["results"]["active_rows"])
    Cols_upper[k] = length(data["results"]["active_cols_upper"])
    Cols_lower[k] = length(data["results"]["active_cols_lower"])
    RoD_Constraint[k] = data["results"]["RoD_Constraint"]
    All_Sets[k] = length(data["results"]["active_sets"])

    unique_rows = unique(vcat(data["results"]["active_rows"]...))
    Row_constraints[k] = length(unique_rows)

    unique_cols_upper = unique(vcat(data["results"]["active_cols_upper"]...))
    Cols_upper_constraints[k] = length(unique_cols_upper)

    unique_cols_lower = unique(vcat(data["results"]["active_cols_lower"]...))
    Cols_lower_constraints[k] = length(unique_cols_lower)

end


# data_file = "intermediate_results/pglib_opf_case240_pserc_iteration700.jld"
#
# data = JLD.load(data_file)
#
# active_rows = data["results"]["active_rows"]
# active_cols_upper = data["results"]["active_cols_lower"]
# active_cols_lower = data["results"]["active_cols_upper"]
#
# num_rows = Vector{Int}(length(active_rows))
# num_cols_upper = Vector{Int}(length(active_cols_upper))
# num_cols_lower = Vector{Int}(length(active_cols_lower))
# unique_rows = Vector{Int}()
# unique_cols_upper = Vector{Int}()
# unique_cols_lower = Vector{Int}()
#
# for i=1:length(active_rows)
#   num_rows[i] = length(active_rows[i])
#   unique_rows = unique([unique_rows; active_rows[i]])
# end
#
# for i=1:length(active_cols_upper)
#   num_cols_upper[i] = length(active_cols_upper[i])
#   unique_cols_upper = unique([unique_cols_upper; active_cols_upper[i]])
# end
#
# for i=1:length(active_cols_lower)
#   num_cols_lower[i] = length(active_cols_lower[i])
#   unique_cols_lower = unique([unique_cols_lower; active_cols_lower[i]])
# end



# f = "results_AC_allRoD_180825/case200_pserc"
# #f = "intermediate_results/pglib_opf_case1951_rte_iteration800"
# data_file = string(f, ".jld")
# data = JLD.load(data_file)
#
# count = collect(1:length(data["results"]["iteration_RoD_sets"]))
# plot(count,[data["results"]["iteration_RoD_sets"], data["results"]["iteration_RoD_constraints"]], layout = (2,1))
