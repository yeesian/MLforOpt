using Ipopt

include("find_active_set.jl")
include("basis_policy.jl")

NLsolver = IpoptSolver()

tol = 1e-5

filename = "data/nesta_case1397sp_eir.m"

active_set = find_active_set(filename, tol, NLsolver)

m = build_basis_policy_model(filename, NLsolver, active_set)

status = solve(m)
