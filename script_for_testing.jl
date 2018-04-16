using Ipopt

include("find_active_set.jl")
include("basis_policy.jl")

NLsolver = IpoptSolver()

tol = 1e-5

filename = "data/nesta_case1397sp_eir.m"

network_data = PowerModels.parse_file(filename)
m = Model(solver = NLsolver)

jm, const_refs, var_refs = post_ac_opf_withref(network_data,m)

active_set = find_active_set(jm, const_refs, var_refs, tol)

m = build_basis_policy_model(filename, NLsolver, active_set)

status = solve(m)
