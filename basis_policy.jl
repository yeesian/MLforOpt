using JuMP
using Ipopt
using PowerModels

"""
Given a JuMP model and a PowerModels network data structure,
Builds an AC-OPF formulation of the given data and returns the JuMP model
The active_set modification only adds constraints corresponding to the active
constraints and sets active variables to their upper and lower bounds.
"""

function post_ac_opf_active_set(data::Dict{String,Any}, NLsolver, active_set)
    active_rows = active_set["active_rows"]
    active_cols_lower = active_set["active_cols_lower"]
    active_cols_upper = active_set["active_cols_upper"]

    col_ctr = 0
    row_ctr = 0

    model = Model(solver = NLsolver)

    ref = PowerModels.build_ref(data)[:nw][0]

    # VARIABLES

    # voltage angles
    @variable(model, va[i in keys(ref[:bus])])
    col_ctr += length(ref[:bus])

    # voltage magnitudes
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
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


    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
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
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) - pd - gs*vm[i]^2
        )

        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) - qd + bs*vm[i]^2
        )
    end


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
        c = branch["b_fr"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @NLconstraint(model, p_fr == g/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+c)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        @NLconstraint(model, p_to == g*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+c)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )


        # Phase Angle Difference Limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        row_ctr += 1
        if row_ctr in active_rows
            @NLconstraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        end
        row_ctr += 1
        if row_ctr in active_rows
            @NLconstraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
        end
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end


    return model
end


function build_basis_policy_model(filename, NLsolver, active_set)
    data = PowerModels.parse_file(filename)
    return post_ac_opf_active_set(data::Dict{String,Any}, NLsolver, active_set)
end
