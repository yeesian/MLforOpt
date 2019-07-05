
# returns a dictionary with the active active_set
# inputs are the filename of the case file, tolerance for lagrange multipliers and
# variables to be considered active, and the non-linear IpoptSolver

function find_active_set(jm, const_refs, var_refs, tol = 1e-5)

    all_var_refs = [var_refs["pg"][:];var_refs["qg"][:]]

    row_duals = JuMP.getdual(const_refs)    # find dual of the constraints

    active_rows = find(abs.(row_duals).>tol)     # non-zero Lagrange multiplier
    active_cols_lower = find(abs.(jm.colVal - jm.colLower).< tol)    # variables that are at lower bound
    active_cols_upper = find(abs.(jm.colVal - jm.colUpper).< tol)    # variables that are at upper bound

    active_set = Dict{String,Vector{Int}}(
        "active_rows" => active_rows,
        "active_cols_lower" => active_cols_lower,
        "active_cols_upper" => active_cols_upper
    )

    return active_set
end



# """
# Given a JuMP model and a PowerModels network data structure,
# Builds an AC-OPF formulation of the given data and returns the JuMP model
# The withref modification also returns the variable and constraint references
# """
#
# function post_ac_opf_withref(data::Dict{String,Any}, model=Model())
#     ref = PowerModels.build_ref(data)[:nw][0]
#
#     @variable(model, va[i in keys(ref[:bus])])
#     @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
#
#     @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
#     @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
#
#     @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
#     @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
#
#     @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
#     @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])
#
#     from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
#     @objective(model, Min,
#         sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
#         sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
#     )
#
#     # @constraintref slack_ref
#     for (i,bus) in ref[:ref_buses]
#         # Reference Bus
#         slack_ref = @constraint(model, va[i] == 0)
#     end
#
#     # constraint reference arrays for bus constraints
#     # @constraintref kcl_p[1:length(ref[:bus])]
#     # @constraintref kcl_q[1:length(ref[:bus])]
#     for (i,bus) in ref[:bus]
#         # Bus KCL
#         @NLconstraint(model,
#             sum(p[a] for a in ref[:bus_arcs][i]) +
#             sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
#             sum(pg[g] for g in ref[:bus_gens][i]) -
#             bus["pd"] - bus["gs"]*vm[i]^2
#         )
#         @NLconstraint(model,
#             sum(q[a] for a in ref[:bus_arcs][i]) +
#             sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
#             sum(qg[g] for g in ref[:bus_gens][i]) -
#             bus["qd"] + bus["bs"]*vm[i]^2
#         )
#     end
#
#     # constraint reference arrays for branch constraints
#     # @constraintref ang_diff_max_ref[1:length(ref[:branch])
#     # @constraintref ang_diff_min_ref[1:length(ref[:branch])
#     # @constraintref sapp_fr_ref[1:length(ref[:branch])]
#     # @constraintref sapp_to_ref[1:length(ref[:branch])]
#     # branch_ctr = 0
#     @constraintref const_refs[1:2*length(ref[:branch])]
#     row_ctr = 0
#     for (i,branch) in ref[:branch]
#         f_idx = (i, branch["f_bus"], branch["t_bus"])
#         t_idx = (i, branch["t_bus"], branch["f_bus"])
#
#         p_fr = p[f_idx]
#         q_fr = q[f_idx]
#         p_to = p[t_idx]
#         q_to = q[t_idx]
#
#         vm_fr = vm[branch["f_bus"]]
#         vm_to = vm[branch["t_bus"]]
#         va_fr = va[branch["f_bus"]]
#         va_to = va[branch["t_bus"]]
#
#         # Line Flow
#         g, b = PowerModels.calc_branch_y(branch)
#         tr, ti = PowerModels.calc_branch_t(branch)
#         g_fr = branch["g_fr"]
#         b_fr = branch["b_fr"]
#         g_to = branch["g_to"]
#         b_to = branch["b_to"]
#         tm = branch["tap"]^2
#
#         # AC Line Flow Constraints
#         @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
#         @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
#
#         @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
#         @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
#
#         # Phase Angle Difference Limit
#
#         @constraint(model, va_fr - va_to <= branch["angmax"])
#         @constraint(model, va_fr - va_to >= branch["angmin"])
#
#         # Apparent Power Limit, From and To
#         row_ctr += 1
#         const_refs[row_ctr] = @NLconstraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
#         row_ctr += 1
#         const_refs[row_ctr] = @NLconstraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
#     end
#
#     for (i,dcline) in ref[:dcline]
#         # DC Line Flow Constraint
#         f_idx = (i, dcline["f_bus"], dcline["t_bus"])
#         t_idx = (i, dcline["t_bus"], dcline["f_bus"])
#
#         @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
#     end
#
#
#     # const_refs = Dict{String,Any}()
#     # const_refs["sapp_fr_ref"] = sapp_fr_ref
#     # const_refs["sapp_to_ref"] = sapp_to_ref
#
#     var_refs = Dict{String,Any}()
#     var_refs["pg"] = pg
#     var_refs["qg"] = qg
#     var_refs["p"] = p
#     var_refs["q"] = q
#     var_refs["va"] = va
#     var_refs["vm"] = vm
#
#     return model, const_refs, var_refs
# end



"""
Given a JuMP model and a PowerModels network data structure,
Builds an AC-OPF formulation of the given data and returns the JuMP model
The withref modification also returns the variable and constraint references
"""

function post_ac_opf_withref_uncertainty(data::Dict{String,Any}, model=Model())
    ref = PowerModels.build_ref(data)[:nw][0]
    nonzeroload = [i for (i,l) in ref[:load] if abs(l["pd"]) > 0.0]
    #randcost = [rand(Uniform(),1) for (i,gen) in ref[:gen]]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    @NLparameter(model, u[i in nonzeroload] == 0)

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
        + sum(qg[i]^2 for (i,gen) in ref[:gen])
    )

    # @constraintref slack_ref
    for (i,bus) in ref[:ref_buses]
        # Reference Bus
        slack_ref = @constraint(model, va[i] == 0)
    end

    # constraint reference arrays for bus constraints
    # @constraintref kcl_p[1:length(ref[:bus])]
    # @constraintref kcl_q[1:length(ref[:bus])]
    # function gamma(ref, load)
    #     abs(ref[:load][load]["pd"]) > 0 ? ref[:load][load]["qd"]/ref[:load][load]["pd"] : 0.0
    # end

    for i in keys(ref[:bus])
        qd = pd = gs = bs = 0.0
        for load in ref[:bus_loads][i]
            pd += ref[:load][load]["pd"]
            qd += ref[:load][load]["qd"]
        end
        for shunt in ref[:bus_shunts][i]
            gs += ref[:shunt][shunt]["gs"]
            bs += ref[:shunt][shunt]["bs"]
        end

        # Bus KCL
        @NLconstraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            pd - gs*vm[i]^2 + sum(u[load] for load in ref[:bus_loads][i] if load in nonzeroload)
        )

        @NLconstraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            qd + bs*vm[i]^2 +
            sum(u[load]*ref[:load][load]["qd"]/ref[:load][load]["pd"] for load in ref[:bus_loads][i] if abs(ref[:load][load]["pd"]) > 0)
        )

    end

    # constraint reference arrays for branch constraints
    # @constraintref ang_diff_max_ref[1:length(ref[:branch])
    # @constraintref ang_diff_min_ref[1:length(ref[:branch])
    # @constraintref sapp_fr_ref[1:length(ref[:branch])]
    # @constraintref sapp_to_ref[1:length(ref[:branch])]
    # branch_ctr = 0
    @constraintref const_refs[1:2*length(ref[:branch])]
    row_ctr = 0
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit

        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        row_ctr += 1
        const_refs[row_ctr] = @NLconstraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        row_ctr += 1
        const_refs[row_ctr] = @NLconstraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end


    # const_refs = Dict{String,Any}()
    # const_refs["sapp_fr_ref"] = sapp_fr_ref
    # const_refs["sapp_to_ref"] = sapp_to_ref

    var_refs = Dict{String,Any}()
    var_refs["pg"] = pg
    var_refs["qg"] = qg
    var_refs["p"] = p
    var_refs["q"] = q
    var_refs["va"] = va
    var_refs["vm"] = vm

    nl_refs = Dict{String,Any}()
    nl_refs["u"] = u

    return model, const_refs, var_refs, nl_refs
end
