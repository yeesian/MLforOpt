"""
Given a (i) JuMP model, (ii) PowerModels network data, and (iii) active
set (to be provided by the callee),
Builds an AC-OPF formulation of the given data and returns the JuMP model
The withref modification also returns the variable and constraint references

We do not provide a default value for active_set for now. As the method
    * fixing the variables to their limits (based on the active set)
    * adding line limits (if they are in the active set)
a sensible default for active set is to have none of the columns, and all of the
rows so that it behaves like post_ac_opf_withref_uncertainty().
"""
function post_ac_opf_active_set_withref_uncertainty(
        data::Dict{String,Any},
        active_set,
        model = Model(),
    )
    active_rows = active_set["active_rows"]
    active_cols_lower = active_set["active_cols_lower"]
    active_cols_upper = active_set["active_cols_upper"]

    col_ctr = 0
    row_ctr = 0

    ref = PowerModels.build_ref(data)[:nw][0]
    # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#L163 at ebc4c9e
        nonzeroload = [i for (i,l) in ref[:load] if abs(l["pd"]) > 0.0]

    # VARIABLES

    # voltage angles
    @variable(model, va[i in keys(ref[:bus])])
    col_ctr += length(ref[:bus])

    # voltage magnitudes
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
        # marked out as corresponding code for setting voltage magnitudes
        # based on active set. (It may be ignored in the future if we only
        # want to deal with active sets for the line limits.)
        ind_lower = [i for i in active_cols_lower if col_ctr+1 <= i <= col_ctr+length(ref[:bus])]
        ind_upper = [i for i in active_cols_upper if col_ctr+1 <= i <= col_ctr+length(ref[:bus])]
        idx = [i for i in keys(ref[:bus])]
        idx_lower = idx[ind_lower-col_ctr]
        idx_upper = idx[ind_upper-col_ctr]
        for i in idx_lower
            JuMP.fix(vm[i],ref[:bus][i]["vmin"])
        end
        for i in idx_upper
            JuMP.fix(vm[i],ref[:bus][i]["vmax"])
        end
    col_ctr += length(ref[:bus])

    # active power
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
        # marked out as corresponding code for setting active power
        # based on active set. (It may be ignored in the future if we only
        # want to deal with active sets for the line limits.)
        ind_lower = [i for i in active_cols_lower if col_ctr+1 <= i <= col_ctr+length(ref[:gen])]
        ind_upper = [i for i in active_cols_upper if col_ctr+1 <= i <= col_ctr+length(ref[:gen])]
        idx = [i for i in keys(ref[:gen])]
        idx_lower = idx[ind_lower-col_ctr]
        idx_upper = idx[ind_upper-col_ctr]
        for i in idx_lower
            JuMP.fix(pg[i],ref[:gen][i]["pmin"])
        end
        for i in idx_upper
            JuMP.fix(pg[i],ref[:gen][i]["pmax"])
        end
    col_ctr += length(ref[:gen])

    # reactive power
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
        # marked out as corresponding code for setting reactive power
        # based on active set. (It may be ignored in the future if we only
        # want to deal with active sets for the line limits.)
        ind_lower = [i for i in active_cols_lower if col_ctr+1 <= i <= col_ctr+length(ref[:gen])]
        ind_upper = [i for i in active_cols_upper if col_ctr+1 <= i <= col_ctr+length(ref[:gen])]
        idx = [i for i in keys(ref[:gen])]
        idx_lower = idx[ind_lower-col_ctr]
        idx_upper = idx[ind_upper-col_ctr]
        for i in idx_lower
            JuMP.fix(qg[i],ref[:gen][i]["qmin"])
        end
        for i in idx_upper
            JuMP.fix(qg[i],ref[:gen][i]["qmax"])
        end
    col_ctr += length(ref[:gen])

    # ac power flows
    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    col_ctr += length(ref[:branch])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    col_ctr += length(ref[:branch])


    # dc power flows (on dc lines)
    @variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    col_ctr += length(ref[:arcs_dc_param])
    @variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])
    col_ctr += length(ref[:arcs_dc_param])

    # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#L178 at ebc4c9e
    # TODO(yeesian): double-check whether we need to increase col_ctr here.
    #   I think we should be fine so long as it is introduced only after all the
    #   variables have been defined (which is assumed to be the case right now).
        @NLparameter(model, u[i in nonzeroload] == 0)

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#L181-185 at ebc4c9e
        @objective(model, Min,
            sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
            sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
            + sum(qg[i]^2 for (i,gen) in ref[:gen])
        )

    for (i,bus) in ref[:ref_buses]
        # Reference Bus
        @constraint(model, va[i] == 0)
    end

    for i in keys(ref[:bus])
        # Bus KCL
        qd = pd = gs = bs = 0.0
        for load in ref[:bus_loads][i]
            pd += ref[:load][load]["pd"]
            qd += ref[:load][load]["qd"]
        end
        for shunt in ref[:bus_shunts][i]
            gs += ref[:shunt][shunt]["gs"]
            bs += ref[:shunt][shunt]["bs"]
        end

        # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#213-226 at ebc4c9e
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

    # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#230-237 at ebc4c9e
        # constraint reference arrays for branch constraints
        # @constraintref ang_diff_max_ref[1:length(ref[:branch])
        # @constraintref ang_diff_min_ref[1:length(ref[:branch])
        # @constraintref sapp_fr_ref[1:length(ref[:branch])]
        # @constraintref sapp_to_ref[1:length(ref[:branch])]
        # branch_ctr = 0
        @constraintref const_refs_powerlimit[1:2*length(ref[:branch])]
        @constraintref const_refs_phaseangle[1:2*length(ref[:branch])]
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

        # AC Line Flow Constraints
        # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#255-266 at ebc4c9e
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]
            tm = branch["tap"]^2

            @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
            @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

            @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
            @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )


        # Phase Angle Difference Limit
        row_ctr += 1
        const_refs_phaseangle[row_ctr] = @constraint(model, va_fr - va_to <= branch["angmax"])
        row_ctr += 1
        const_refs_phaseangle[row_ctr] = @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        # differs from post_ac_opf_withref_uncertainty() in find_active_set.jl#274-277 at ebc4c9e
        # * here is where the active sets stuff happens.
        # * row_ctr should increase regardless of whether the if condition holds,
        #   so that it tracks the corresponding row in the original model.
        row_ctr += 1
        if row_ctr in active_rows
            const_refs_powerlimit[row_ctr] = @NLconstraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        end
        row_ctr += 1
        if row_ctr in active_rows
            const_refs_powerlimit[row_ctr] = @NLconstraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
        end
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    # mirrors post_ac_opf_withref_uncertainty() in find_active_set.jl#289-304 at ebc4c9e
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

        return model, const_refs_phaseangle, const_refs_powerlimit, var_refs, nl_refs
end
