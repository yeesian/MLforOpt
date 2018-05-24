using Ipopt
using Distributions

include("find_active_set.jl")
include("basis_policy.jl")

NLsolver = IpoptSolver()

tol = 1e-5

#filename = "data/nesta_case1397sp_eir.m"
filename = "pglib-opf/pglib_opf_case5_pjm.m"

network_data = PowerModels.parse_file(filename)
ref = PowerModels.build_ref(network_data)[:nw][0]
nonzeroindices = [i for (i,bus) in ref[:bus] if bus["pd"]>1e-5]

m = Model(solver = NLsolver)

#jm, const_refs, var_refs = post_ac_opf_withref(network_data,m)
jm, const_refs, var_refs, nl_refs = post_ac_opf_withref_uncertainty(network_data,m)

status = solve(jm)

# Constructing distribution for samples
sigma = 0.1
load = [ref[:bus][i]["pd"] for i in nonzeroindices]
w = Distributions.MvNormal(
    zeros(length(nonzeroindices)),
    diagm((sigma*load).^2)
)

# Run OPF to get active sets
results=[]
active_rows = Set()
active_cols_upper = Set()
active_cols_lower = Set()

for i = 1:4#:m+W
    # 1. Generate sample
    w_sample = rand(w,1)
    # 2. Fix NL parameters
    for (k,j) in enumerate(nonzeroindices)
        setvalue(nl_refs["u"][j], w_sample[k])
    end
    # 3. Get active set
    active_set = find_active_set(jm, const_refs, var_refs, tol)
    push!(active_rows, active_set["active_rows"])
    push!(active_cols_upper, active_set["active_cols_upper"])
    push!(active_cols_lower, active_set["active_cols_lower"])
    # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)

end


status = solve(jm)

#active_set = find_active_set(jm, const_refs, var_refs, tol)

km = build_basis_policy_model(filename, NLsolver, active_set)
