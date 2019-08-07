using JLD, ProgressMeter, PowerModels, JuMP, Ipopt, Distributions, Base.Test

include("find_active_set.jl") # for post_ac_opf_withref_uncertainty
include("basis_policy_withref_uncertainty.jl") # for post_ac_opf_active_set_withref_uncertainty

nlsolver = Ipopt.IpoptSolver(print_level=0)
tol = 1e-5

f = "results_AC_allRoD_180825/case162_ieee_dtc" # 0.51 versus 0.65 sec
f = "results_AC_allRoD_180825/case200_pserc" # 0.74 versus 0.52 sec
f = "intermediate_results/pglib_opf_case240_pserc_iteration700" # 3.59 v.s. 1.34 sec
f = "intermediate_results/pglib_opf_case300_ieee_iteration1000" # 1.11 v.s. 0.85 sec
f = "intermediate_results/pglib_opf_case1888_rte_iteration1100" # 10.19 v.s. 7.05 sec
f = "intermediate_results/pglib_opf_case1951_rte_iteration800" # 12.18 v.s. 7.88 sec (~4 violated constraint for both models)
f = "intermediate_results/pglib_opf_case2737sop_k_iteration3600" # 12.91 v.s. 10.03 sec
f = "intermediate_results/pglib_opf_case2848_rte_iteration1000" # 16.03 v.s. 11.14 sec
f = "intermediate_results/pglib_opf_case2869_pegase_iteration700" # 27.92 v.s. 14.05 sec
f = "intermediate_results/pglib_opf_case6468_rte_iteration500" # 45.42 v.s. 30.82 seconds

data_file = string(f, ".jld")
data = JLD.load(data_file)
results = data["results"]

# Loading data about the system
filename = results["filename"]
network_data = PowerModels.parse_file(filename)
ref = PowerModels.build_ref(network_data)[:nw][0]
nonzeroindices = [i for (i,l) in ref[:load] if abs(l["pd"]) > 0.0]

# create uncertainty distribution
sigma = 0.1
load = [ref[:load][i]["pd"] for i in nonzeroindices]
w = Distributions.MvNormal(
    zeros(length(nonzeroindices)),
    diagm((sigma*load).^2)
)

# full ACOPF model
jm, const_refs, var_refs, nl_refs = post_ac_opf_withref_uncertainty(
    network_data,
    JuMP.Model(solver = nlsolver)
)
JuMP.solve(jm)

# for i in 1:(# of iterations)

w_sample = rand(w,1) # 1. Generate sample
for j in eachindex(nonzeroindices) # 2. Fix NL parameters
    setvalue(nl_refs["u"][nonzeroindices[j]], w_sample[j])
end
@time status = JuMP.solve(jm)
@test status == :Optimal
active_set = find_active_set(jm, const_refs, var_refs, tol) # 3. Get active set

# set up active set model
acopf_model = post_ac_opf_active_set_withref_uncertainty(
    network_data,
    Dict("active_rows" => active_set["active_rows"]),
    JuMP.Model(solver = nlsolver)
)
JuMP.solve(acopf_model.model)
for j in eachindex(nonzeroindices) # 2. Fix NL parameters
    setvalue(acopf_model.nl_refs["u"][nonzeroindices[j]], w_sample[j])
end
@time status = JuMP.solve(acopf_model.model)
@test status == :Optimal
new_active_set = find_active_set(
    acopf_model.model,
    acopf_model.const_refs_powerlimit,
    acopf_model.var_refs,
    tol
)

# following https://discourse.julialang.org/t/getting-constraint-values-using-jump/1996/7
d = JuMP.NLPEvaluator(jm)
MathProgBase.initialize(d, [:Grad])
g = zeros(MathProgBase.numconstr(d.m))
b = JuMP.constraintbounds(jm)
function find_violated_constraints(x)
    MathProgBase.eval_g(d, g, x)
    return find([g-b[1];b[2]-g] .< -1e-4)
end
@test length(find_violated_constraints(JuMP.internalmodel(jm).inner.x)) == 0
@test length(find_violated_constraints(JuMP.internalmodel(acopf_model.model).inner.x)) == 0

# @test active_set == new_active_set
