# Summarizing and plotting results for the AC OPF simulations

using JLD
using ProgressMeter
using PowerModels


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

    # Defining necessary data structures
    unique_rows = Vector{Int}()
    unique_cols_upper = Vector{Int}()
    unique_cols_lower = Vector{Int}()

    # Loading results
    data_file = string(f, ".jld")
    data = JLD.load(data_file)

    # Loading data about the system
    filename = data["results"]["filename"]
    network_data = PowerModels.parse_file(filename)
    ref = PowerModels.build_ref(network_data)[:nw][0]
    nL = length(ref[:branch])
    nG = length(ref[:gen])
    nB = length(ref[:bus])

    # Variable definition:
    # va = 1, ... , a            a = nB
    # vm = a+1, ... , b          b = 2*nB
    # pg = b+1, ... , c          c = 2*nB+nG
    # qg = c+1, ... , d          d = 2*nB+2*nG
    # p  = d+1, ... , e          e = 2*nB+2*nG+nL
    # q  = e+1, ... , f          f = 2*nB+2*nG+2*nL
    # pdc  = f+1, ... , g        g = 2*nB+2*nG+2*nL+nDC
    # qdc  = g+1, ... , h        h = 2*nB+2*nG+2*nL+2*nDC
    a = nB
    b = 2*nB
    c = 2*nB+nG
    d = 2*nB+2*nG
    e = 2*nB+2*nG+nL
    f = 2*nB+2*nG+2*nL
    g = 2*nB+2*nG+2*nL+nDC
    h = 2*nB+2*nG+2*nL+2*nDC


    # Summary data
    K_M[k]  = data["K_M"]
    M[k]    = data["M"]
    W[k]    = data["W"]

    RoD[k]  = data["RoD"]
    RoD_Constraint[k] = data["results"]["RoD_Constraint"]

    # Figuring out constraint combinations
    All_Sets[k] = length(data["results"]["active_sets"])
    Rows[k] = length(data["results"]["active_rows"])
    Cols_upper[k] = length(data["results"]["active_cols_upper"])
    Cols_lower[k] = length(data["results"]["active_cols_lower"])


    unique_rows = unique(vcat(data["results"]["active_rows"]...))
    Row_constraints[k] = length(unique_rows)

    unique_cols_upper = unique(vcat(data["results"]["active_cols_upper"]...))
    Cols_upper_constraints[k] = length(unique_cols_upper)

    unique_cols_lower = unique(vcat(data["results"]["active_cols_lower"]...))
    Cols_lower_constraints[k] = length(unique_cols_lower)

end




# Variable definition:

va = 1, ... , a            a = nB
vm = a+1, ... , b          b = 2*nB
pg = b+1, ... , c          c = 2*nB+nG
qg = c+1, ... , d          d = 2*nB+2*nG
p  = d+1, ... , e          e = 2*nB+2*nG+nL
q  = e+1, ... , f          f = 2*nB+2*nG+2*nL
pdc  = f+1, ... , g        g = 2*nB+2*nG+2*nL+nDC
qdc  = g+1, ... , h        h = 2*nB+2*nG+2*nL+2*nDC



# @variable(model, va[i in keys(ref[:bus])])
# @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
#
# @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
# @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
#
# @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
# @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
#
# @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
# @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])






# Printing results to tables
results_combinations=[M+W All_Sets Cols_lower Cols_lower_constraints Cols_upper Cols_upper_constraints Rows Row_constraints]

#results_numbers=[M+W All_Sets Cols_lower Cols_lower_constraints Cols_upper Cols_upper_constraints Rows Row_constraints]

data_file = "intermediate_results/pglib_opf_case240_pserc_iteration700.jld"

data = JLD.load(data_file)
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
