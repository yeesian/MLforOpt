function post_ac_opf_active_set_withref_uncertainty(
        data::Dict{String,Any},
        NLsolver,
        active_set
    )
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

# ======================================================================

"""
Run learning algorithm as a streaming algorithm for a given system and pre-computed set of samples,
using rate of discovery of active sets as the stopping criterion (as it is the most conservative one),
but also recording rate of discovery for all the other quantities:
New active sets, new active constraints, new active rows, new acitve columns (upper and lower)

Inputs:\\
alpha   - maximum unobserved mass \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1\\

Outputs:\\
"""

function RunStreamingAlgorithmAC_AllRoD(alpha, delta, epsilon, gamma, Minitial, filename, NLsolver; maxsamples = 50000, tol = 1e-5, sigma=0.03)

   # Keeping track of infeasible samples
   infeasible_samples = 0

   # Evaluate stopping criterion
   threshold = StoppingCriterion(alpha, epsilon)

   #Find maximum window length (since the number of samples is limited)
   #W_max = WindowSize(delta, epsilon, gamma, length(scenarios.whichbasis[:,1]))

   # Initialize
   W = WindowSize(delta, epsilon, gamma, Minitial)

   # Parse data
   network_data = PowerModels.parse_file(filename)
   ref = PowerModels.build_ref(network_data)[:nw][0]
   nonzeroload = [i for (i,l) in ref[:load] if abs(l["pd"]) > 0.0]
   #nonzeroindices = [i for (i,loads) in ref[:bus_loads] if length(loads) > 0]


   # Posting model with NL parameters for omega
   m_init = Model(solver = NLsolver)
   jm, const_refs, var_refs, nl_refs = post_ac_opf_withref_uncertainty(network_data,m_init)

   # Constructing distribution for samples
   #sigma = 0.03
   load = [l["pd"] for (i,l) in ref[:load] if abs(l["pd"]) > 0.0]
   w = Distributions.MvNormal(
       zeros(length(load)),
       diagm((sigma*abs.(load)).^2)
   )

   # Run OPF to get active sets
   active_rows = Set{Vector{Int}}()
   active_cols_upper = Set{Vector{Int}}()
   active_cols_lower = Set{Vector{Int}}()

   dict_active_rows = Dict{Vector{Int}, Int}()
   dict_active_cols_upper = Dict{Vector{Int}, Int}()
   dict_active_cols_lower = Dict{Vector{Int}, Int}()

   observed_active_sets = Set{Int}()
   all_active_sets = Set{Tuple{Int,Int,Int}}()
   active_sets = Tuple{Int,Int,Int}[]
   active_index = Dict{Tuple{Int,Int,Int}, Int}()

   observedsets = Set{Int}()
   windowcount = Dict{Int,Int}() # active set index => # of times it appears in the window
   q = DataStructures.Queue(Int)

   #scenario = Dict{String,Any}[]
   scenario_realization = Matrix{Float64}[]
   scenario_active_set = Int[]

   iteration_RoD_sets = Float64[]
   iteration_RoD_constraints = Float64[]

   # Initialization

   ProgressMeter.@showprogress 1 for i = 1:(Minitial+W)#:m+W

       w_sample = Array{Float64}

       for nthtry=1:100
           # 1. Generate sample
           w_sample = rand(w,1)
           # 2. Fix NL parameters
           for (k,j) in enumerate(nonzeroload)
               setvalue(nl_refs["u"][j], w_sample[k])
           end
           # 3. Solve problem
           status = solve(jm)
           if status == :Optimal
               if nthtry > 1
                   println("It took $nthtry tries to solve the problem")
               end
               break
           else
               infeasible_samples = infeasible_samples + 1
           end
       end
       # 3. Get active set
       active_set = find_active_set(jm, const_refs, var_refs, tol)

       # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)

       # If a new active row or column is observed, add them to the set and create a new index
       if !in(active_set["active_rows"], active_rows)
           push!(active_rows, active_set["active_rows"])
           dict_active_rows[active_set["active_rows"]] = length(active_rows)
       end
       if !in(active_set["active_cols_upper"], active_cols_upper)
           push!(active_cols_upper, active_set["active_cols_upper"])
           dict_active_cols_upper[active_set["active_cols_upper"]] = length(active_cols_upper)
       end
       if !in(active_set["active_cols_lower"], active_cols_lower)
           push!(active_cols_lower, active_set["active_cols_lower"])
           dict_active_cols_lower[active_set["active_cols_lower"]] = length(active_cols_lower)
       end

       # Put together the marker for the active set discovered now
       new_active_set = (
           dict_active_rows[active_set["active_rows"]],
           dict_active_cols_upper[active_set["active_cols_upper"]],
           dict_active_cols_lower[active_set["active_cols_lower"]]
       )

       # If the current active set has not been discovered
       if !in(new_active_set, all_active_sets)
           push!(all_active_sets, new_active_set)  # add the new active set to the (unordered) SET of discovered active sets
           push!(active_sets, new_active_set)      # add the new active set to the (ordered) VECTORS of active sets
           active_index[new_active_set] = length(active_sets)  # make a new index for the new active set
       end

       # Add the active set to the queue (the list of all active sets)
       DataStructures.enqueue!(q, active_index[new_active_set])

       push!(scenario_active_set, active_index[new_active_set])
       push!(scenario_realization, w_sample)

       #push!(scenario, Dict{String,Any}("realization"=>w_sample, "active_set"=>active_index[new_active_set])))

       # Update the count of how frequently the active sets appear in the queue
       windowcount[active_index[new_active_set]] = get(windowcount, active_index[new_active_set], 0) + 1
   end

   for i in 1:Minitial
       curr_active_set = dequeue!(q)                   # remove the sample
       push!(observed_active_sets, curr_active_set)    # put active set into observed active set
       windowcount[curr_active_set] -= 1               # update windowcount
       @assert windowcount[curr_active_set] >= 0
   end

   # computing rate of discovery
   numerator = 0
   denominator = sum(l for (i,l) in windowcount)
   @assert denominator > 0
   for (as,v) in windowcount
       if !in(as, observed_active_sets)
           numerator += v
       end
   end

   RoD = (numerator / denominator)
   println("rate of discovery: $(RoD)")
   push!(iteration_RoD_sets, RoD)



   # NEW PART ABOUT DISCOVERING CONSTRAINTS

   # General idea: Record the discovery of NEW ACTIVE CONSTRAINTS
   # Consider all active sets where the active constraints are
   # subsets of the union of all active constraints in the discovered active
   # sets to also be discovered

   # Computing the set of all active constraints that have been discovered
   # - compute the union of all active constraints inside the observed_active_sets
   # - for each active set in the window which has not yet been discovered:
   # -> Check if the active constraints are a subset of the active constraints within the observed active sets
   # -> If it is a subset: Mark as "discovered"
   # -> If it is not a subset: Leave as "undiscovered"
   # - Compute rate of discovery of active sets that have new active constraints



   # Initializing the constraint learning
   list_active_rows = Array{Vector{Int}}(length(active_rows))
   list_active_cols_upper = Array{Vector{Int}}(length(active_cols_upper))
   list_active_cols_lower = Array{Vector{Int}}(length(active_cols_lower))

   for (v,i) in dict_active_rows; list_active_rows[i] = v end
   for (v,i) in dict_active_cols_upper; list_active_cols_upper[i] = v end
   for (v,i) in dict_active_cols_lower; list_active_cols_lower[i] = v end

   discovered_rows       = zeros(list_active_rows,Bool)
   discovered_cols_upper = zeros(list_active_cols_upper,Bool)
   discovered_cols_lower = zeros(list_active_cols_lower,Bool)
   discovered_sets        = zeros(Bool, length(active_sets))

   set_active_rows       = Set{Int}()
   set_active_cols_upper = Set{Int}()
   set_active_cols_lower = Set{Int}()

   observed_active_constraints = Set{Int}()

   # Initialization: Add constraints from the first M active sets

   println("=== Identifying known active sets and adding their constraints ===")
   known_active_sets = unique(scenario_active_set[1:Minitial])
   discovered_sets[known_active_sets] = true

   println("ACTIVE ROWS")
   known_active_rows       = unique([active_sets[known_active_sets][i][1] for i in known_active_sets])
   set_active_rows         = add_new_constraints(list_active_rows[known_active_rows], set_active_rows)
   discovered_rows[known_active_rows] = true

   println("ACTIVE COLS UPPER")
   known_active_cols_upper = unique([active_sets[known_active_sets][i][2] for i in known_active_sets])
   set_active_cols_upper   = add_new_constraints(list_active_cols_upper[known_active_cols_upper], set_active_cols_upper)
   discovered_cols_upper[known_active_cols_upper] = true

   println("ACTIVE COLS LOWER")
   known_active_cols_lower = unique([active_sets[known_active_sets][i][3] for i in known_active_sets])
   set_active_cols_lower   = add_new_constraints(list_active_cols_lower[known_active_cols_lower], set_active_cols_lower)
   discovered_cols_lower[known_active_cols_lower] = true

   # Check if other sets are covered by the discovered ones
   println("=== Update discovery status for all collections ===")

   println("ACTIVE ROWS")
   discovered_rows = check_discovery(list_active_rows,
                                     set_active_rows,
                                     discovered_rows)

   println("ACTIVE COLS UPPER")
   discovered_cols_upper = check_discovery(list_active_cols_upper,
                                           set_active_cols_upper,
                                           discovered_cols_upper)

   println("ACTIVE COLS LOWER")
   discovered_cols_lower = check_discovery(list_active_cols_lower,
                                            set_active_cols_lower,
                                            discovered_cols_lower)

   #println(length(active_sets))
   #println(active_sets)
   #println("status discovered active sets $discovered_sets")

   # Update list of discovered active active sets
   # MAKE THIS PRETTIER WHEN YOU GET INTERNET...
   for ind_act=1:length(active_sets)
       all_discovered = discovered_rows[active_sets[ind_act][1]] == true &&
                        discovered_cols_upper[active_sets[ind_act][2]] == true &&
                        discovered_cols_lower[active_sets[ind_act][3]] == true
       if discovered_sets[ind_act] == false && all_discovered == true
           discovered_sets[ind_act] = true
           println("All constraints in $ind_act are discovered")
           println("This set is $(active_sets[ind_act])")
       end
   end

   # for ind_act=1:length(active_sets)
   #     if discovered_sets[ind_act] == false
   #         if discovered_rows[active_sets[ind_act][1]] == true
   #             if discovered_cols_upper[active_sets[ind_act][2]] == true
   #                 if discovered_cols_lower[active_sets[ind_act][3]] == true
   #                     discovered_sets[ind_act] = true
   #                     println("All constraints in $ind_act are discovered")
   #                     println("This set is $(active_sets[ind_act])")
   #                 end
   #             end
   #         end
   #     end
   # end

   # Gather the indices of the undiscovered sets
   for k=1:length(discovered_sets)
       if discovered_sets[k] == true
           push!(observed_active_constraints, active_index[active_sets[k]])
       end
   end

   #println(observed_active_constraints)

   # computing rate of discovery
   numerator_constr = 0
   denominator_constr = sum(l for (i,l) in windowcount)
   @assert denominator_constr > 0
   for (as,v) in windowcount
       if !in(as, observed_active_constraints)
           numerator_constr += v
       end
   end

   RoD_Constraint = numerator_constr/denominator_constr
   println("rate of discovery of new constraints: $(RoD_Constraint)")
   push!(iteration_RoD_constraints, RoD_Constraint)

   k_m=length(observed_active_sets)
   println("Iteration: 0 M: $Minitial W: $W K_M: $k_m")


   # Streaming

   M = Minitial;
   ProgressMeter.@showprogress 1 for jiter = 1:maxsamples
       M += 1;
       Wold = W;
       W = WindowSize(delta, epsilon, gamma, M)

       k_m=length(observed_active_sets)

       # add new samples to W
       for i = 1:W-Wold+1 #add at least one sample:
                          #replace the which will be moved in to M,
                          #plus add new samples to cover an increase in the size of W
          w_sample = Array{Float64}

          for nthtry=1:100
              # 1. Generate sample
              w_sample = rand(w,1)
              # 2. Fix NL parameters
              for (k,j) in enumerate(nonzeroload)
                  setvalue(nl_refs["u"][j], w_sample[k])
              end
              # 3. Solve problem
              status = solve(jm)
              if status == :Optimal
                  if nthtry > 1
                      println("It took $nthtry tries to solve the problem")
                  end
                  break
              else
                  infeasible_samples = infeasible_samples + 1
              end
          end

           # 3. Get active set
           active_set = find_active_set(jm, const_refs, var_refs, tol)

           # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)

           # If a new active row or column is observed, add them to the set and create a new index
           if !in(active_set["active_rows"], active_rows)
               push!(active_rows, active_set["active_rows"])
               dict_active_rows[active_set["active_rows"]] = length(active_rows)
           end
           if !in(active_set["active_cols_upper"], active_cols_upper)
               push!(active_cols_upper, active_set["active_cols_upper"])
               dict_active_cols_upper[active_set["active_cols_upper"]] = length(active_cols_upper)
           end
           if !in(active_set["active_cols_lower"], active_cols_lower)
               push!(active_cols_lower, active_set["active_cols_lower"])
               dict_active_cols_lower[active_set["active_cols_lower"]] = length(active_cols_lower)
           end

           # Put together the marker for the active set discovered now
           new_active_set = (
               dict_active_rows[active_set["active_rows"]],
               dict_active_cols_upper[active_set["active_cols_upper"]],
               dict_active_cols_lower[active_set["active_cols_lower"]]
           )

           # If the current active set has not been discovered
           if !in(new_active_set, all_active_sets)
               push!(all_active_sets, new_active_set)  # add the new active set to the (unordered) SET of discovered active sets
               push!(active_sets, new_active_set)      # add the new active set to the (ordered) VECTORS of active sets
               active_index[new_active_set] = length(active_sets)  # make a new index for the new active set
           end

           # Add the active set to the queue (the list of all active sets)
           DataStructures.enqueue!(q, active_index[new_active_set])

           push!(scenario_active_set, active_index[new_active_set])
           push!(scenario_realization, w_sample)

           # Update the count of how frequently the active sets appear in the queue
           windowcount[active_index[new_active_set]] = get(windowcount, active_index[new_active_set], 0) + 1
       end

       # pushing one sample from W into M
       curr_active_set = dequeue!(q)                   # remove the sample
       push!(observed_active_sets, curr_active_set)    # put active set into observed active set
       windowcount[curr_active_set] -= 1               # update windowcount
       @assert windowcount[curr_active_set] >= 0

       k_m=length(observed_active_sets)
       println("Iteration: $jiter M: $M W: $W K_M: $k_m")

       # computing rate of discovery
       numerator = 0
       denominator = sum(l for (i,l) in windowcount)
       @assert denominator > 0
       for (as,v) in windowcount
           if !in(as, observed_active_sets)
               numerator += v
           end
       end
       RoD = (numerator / denominator)

       println("rate of discovery: $(RoD)")

       push!(iteration_RoD_sets, RoD)



       # Updating the datastructures (this is very unelegant)
       if length(active_rows) > length(list_active_rows)
           for ind_over = length(list_active_rows)+1:length(active_rows)
               for (v,i) in dict_active_rows
                   if i == ind_over
                       println("Pushing into collection element $i")
                       push!(list_active_rows, v)
                       println("New length of collection: $(length(list_active_rows))")
                   end
               end
           end
           #for (v,i) in dict_active_rows; list_active_rows[i] = v end
           # old_length = length(list_active_rows)
           # for (v,i) in dict_active_rows
           #     if i>old_length
           #         println("Old length of collection: $(length(list_active_cols_upper))")
           #         println("Pushing into collection element $i")
           #         #push!(list_active_rows, v)
           #         list_active_rows[i]=v
           #     end
           # end
           newentry=length(active_rows)-length(discovered_rows)
           for ind = 1:newentry
               push!(discovered_rows, false)
           end
           @assert(length(list_active_rows)==length(discovered_rows))
       end

       if length(active_cols_upper) > length(list_active_cols_upper)
           println("size of collection: $(length(list_active_cols_upper))")
           println("size of dict: $(length(dict_active_cols_upper))")
           println("size of discovered vector: $(length(discovered_cols_upper))")
           for ind_over = length(list_active_cols_upper)+1:length(active_cols_upper)
               for (v,i) in dict_active_cols_upper
                   if i == ind_over
                       println("Pushing into collection element $i")
                       push!(list_active_cols_upper, v)
                       println("New length of collection: $(length(list_active_cols_upper))")
                   end
               end
           end
           # old_length = length(list_active_cols_upper)
           # for (v,i) in dict_active_cols_upper
           #     println("for element $i, is i>length(list_active_cols_upper)? $(i>length(list_active_cols_upper)) ")
           #     if i>old_length
           #         println("Old length of collection: $(length(list_active_cols_upper))")
           #         println("Pushing into collection element $i")
           #         #push!(list_active_cols_upper, v)
           #         list_active_cols_upper[i]=v
           #         println("New length of collection: $(length(list_active_cols_upper))")
           #     end
           # end
           newentry=length(active_cols_upper)-length(discovered_cols_upper)
           for ind = 1:newentry
               push!(discovered_cols_upper, false)
           end
           println("size of collection: $(length(list_active_cols_upper))")
           println("size of dict: $(length(dict_active_cols_upper))")
           println("size of discovered vector: $(length(discovered_cols_upper))")
           @assert(length(list_active_cols_upper)==length(discovered_cols_upper))
       end

       if length(active_cols_lower) > length(list_active_cols_lower)
           for ind_over = length(list_active_cols_lower)+1:length(active_cols_lower)
               for (v,i) in dict_active_cols_lower
                   if i == ind_over
                       println("Pushing into collection element $i")
                       push!(list_active_cols_lower, v)
                       println("New length of collection: $(length(list_active_cols_lower))")
                   end
               end
           end
           #for (v,i) in dict_active_cols_lower; list_active_cols_lower[i] = v end
           # old_length = length(list_active_cols_lower)
           # for (v,i) in dict_active_cols_lower
           #     if i>old_length
           #         println("Old length of collection: $(length(list_active_cols_lower))")
           #         println("Pushing into collection element $i")
           #         #push!(list_active_cols_lower, v)
           #         list_active_cols_lower[i]=v
           #         println("New length of collection: $(length(list_active_cols_lower))")
           #     end
           # end
           newentry = length(active_cols_lower)-length(discovered_cols_lower)
           for ind = 1:newentry
               push!(discovered_cols_lower, false)
           end
           @assert(length(list_active_cols_lower)==length(discovered_cols_lower))
       end

       if length(active_sets) > length(discovered_sets)
           for ind = 1:length(active_sets)-length(discovered_sets)
               push!(discovered_sets, false)
           end
       end


       # If the current active set was not discovered, make an update
       #if !in(curr_active_set, observed_active_constraints)

           # println(active_sets[curr_active_set])
           # println(curr_active_set)
           # println(active_sets[curr_active_set][1])
           # println(active_sets[curr_active_set][2])
           # println(active_sets[curr_active_set][3])

           # println("=== Adding constraints rom current active set ===")
           #println("ACTIVE ROWS")
           set_active_rows = add_new_constraints(list_active_rows[active_sets[curr_active_set][1]], set_active_rows)
           discovered_rows[active_sets[curr_active_set][1]] = true

           #println("ACTIVE COLS UPPER")
           set_active_cols_upper   = add_new_constraints(list_active_cols_upper[active_sets[curr_active_set][2]], set_active_cols_upper)
           discovered_cols_upper[active_sets[curr_active_set][2]] = true

           #println("ACTIVE COLS LOWER")
           set_active_cols_lower   = add_new_constraints(list_active_cols_lower[active_sets[curr_active_set][3]], set_active_cols_lower)
           discovered_cols_lower[active_sets[curr_active_set][3]] = true

           # Check if other sets are covered by the discovered ones
           # println("=== Update discovery status for all collections ===")

           #println("ACTIVE ROWS")
           discovered_rows = check_discovery(list_active_rows,
                                             set_active_rows,
                                             discovered_rows)

           #println("ACTIVE COLS UPPER")
           discovered_cols_upper = check_discovery(list_active_cols_upper,
                                                    set_active_cols_upper,
                                                    discovered_cols_upper)

           #println("ACTIVE COLS LOWER")
           discovered_cols_lower = check_discovery(list_active_cols_lower,
                                                    set_active_cols_lower,
                                                    discovered_cols_lower)

       #
       # #end
       #

       # Update list of discovered active active sets
       # MAKE THIS PRETTIER WHEN YOU GET INTERNET...
       for ind_act=1:length(active_sets)
           all_discovered = discovered_rows[active_sets[ind_act][1]] == true &&
                            discovered_cols_upper[active_sets[ind_act][2]] == true &&
                            discovered_cols_lower[active_sets[ind_act][3]] == true
           if discovered_sets[ind_act] == false && all_discovered == true
               discovered_sets[ind_act] = true
               println("All constraints in $ind_act are discovered")
               println("This set is $(active_sets[ind_act])")
           end
       end

       #
       # for ind_act=1:length(active_sets)
       #     if discovered_sets[ind_act] == false
       #         if discovered_rows[active_sets[ind_act][1]] == true
       #             if discovered_cols_upper[active_sets[ind_act][2]] == true
       #                 if discovered_cols_lower[active_sets[ind_act][3]] == true
       #                     discovered_sets[ind_act] = true
       #                     println("All constraints in $ind_act are discovered")
       #                     println("This set is $(active_sets[ind_act])")
       #                 end
       #             end
       #         end
       #     end
       # end

       # Gather the indices of the undiscovered sets
       for k=1:length(discovered_sets)
           if discovered_sets[k] == true
               push!(observed_active_constraints, active_index[active_sets[k]])
           end
       end

       println(observed_active_constraints)

       # computing rate of discovery
       numerator_constr = 0
       denominator_constr = sum(l for (i,l) in windowcount)
       @assert denominator_constr > 0
       for (as,v) in windowcount
           if !in(as, observed_active_constraints)
               numerator_constr += v
           end
       end

       RoD_Constraint = numerator_constr/denominator_constr
       println("rate of discovery of new constraints: $(RoD_Constraint)")
       push!(iteration_RoD_constraints, RoD_Constraint)

       println()



       if RoD <= threshold

           K_M = length(observed_active_sets)

           # # Create inverse mapping from the index to the active_rows, etc
           # list_active_rows = Array{Vector{Int}}(length(active_rows))
           # list_active_cols_upper = Array{Vector{Int}}(length(active_cols_upper))
           # list_active_cols_lower = Array{Vector{Int}}(length(active_cols_lower))
           # for (v,i) in dict_active_rows; list_active_rows[i] = v end
           # for (v,i) in dict_active_cols_upper; list_active_cols_upper[i] = v end
           # for (v,i) in dict_active_cols_lower; list_active_cols_lower[i] = v end

           results = Dict{String,Any}(
               "alpha" => alpha,
               "delta" => delta,
               "epsilon" => epsilon,
               "gamma" => gamma,
               "Minitial" => Minitial,
               "filename" => filename,

               "active_rows" => list_active_rows,
               "active_cols_upper" => list_active_cols_upper,
               "active_cols_lower" => list_active_cols_lower,

               "active_sets" => active_sets,
               "scenario_active_set" => scenario_active_set,
               "scenario_realization" => scenario_realization,

               "observed_active_sets" => observed_active_sets,

               "infeasible_samples" => infeasible_samples,
               "M" => M,
               "W" => W,
               "RoD" => RoD,
               "K_M" => K_M,

               "sigma" => sigma,

               "discovered_rows" => discovered_rows,
               "discovered_sets" => discovered_sets,
               "discovered_cols_upper" => discovered_cols_upper,
               "discovered_cols_lower" => discovered_cols_lower,

               "set_active_rows" => set_active_rows,
               "set_active_cols_lower" => set_active_cols_lower,
               "set_active_cols_upper" => set_active_cols_upper,

               "RoD_Constraint" => RoD_Constraint,

               "iteration_RoD_sets" => iteration_RoD_sets,
               "iteration_RoD_constraints" => iteration_RoD_constraints
           )

           return M, W, RoD, K_M, results
       end

       if mod(jiter,100)==0

           println("Printing intermediate results at iteration $jiter")

           K_M = length(observed_active_sets)

           # # Create inverse mapping from the index to the active_rows, etc
           # list_active_rows = Array{Vector{Int}}(length(active_rows))
           # list_active_cols_upper = Array{Vector{Int}}(length(active_cols_upper))
           # list_active_cols_lower = Array{Vector{Int}}(length(active_cols_lower))
           # for (v,i) in dict_active_rows; list_active_rows[i] = v end
           # for (v,i) in dict_active_cols_upper; list_active_cols_upper[i] = v end
           # for (v,i) in dict_active_cols_lower; list_active_cols_lower[i] = v end

           results = Dict{String,Any}(
               "alpha" => alpha,
               "delta" => delta,
               "epsilon" => epsilon,
               "gamma" => gamma,
               "Minitial" => Minitial,
               "filename" => filename,

               "active_rows" => list_active_rows,
               "active_cols_upper" => list_active_cols_upper,
               "active_cols_lower" => list_active_cols_lower,

               "active_sets" => active_sets,
               "scenario_active_set" => scenario_active_set,
               "scenario_realization" => scenario_realization,

               "observed_active_sets" => observed_active_sets,

               "infeasible_samples" => infeasible_samples,
               "M" => M,
               "W" => W,
               "RoD" => RoD,
               "K_M" => K_M,

               "sigma" => sigma,

               "discovered_rows" => discovered_rows,
               "discovered_sets" => discovered_sets,
               "discovered_cols_upper" => discovered_cols_upper,
               "discovered_cols_lower" => discovered_cols_lower,

               "set_active_rows" => set_active_rows,
               "set_active_cols_lower" => set_active_cols_lower,
               "set_active_cols_upper" => set_active_cols_upper,

               "RoD_Constraint" => RoD_Constraint,

               "iteration_RoD_sets" => iteration_RoD_sets,
               "iteration_RoD_constraints" => iteration_RoD_constraints
           )

           JLD.save("intermediate_results/$(network_data["name"])_iteration$jiter.jld", "results", results, "M", M, "W", W, "RoD", RoD, "K_M", K_M)

       end

   end

   println("Not enough samples available, algorithm did not terminate!")
   K_M = length(observed_active_sets)
   # Create inverse mapping from the index to the active_rows, etc
   list_active_rows = Array{Vector{Int}}(length(active_rows))
   list_active_cols_upper = Array{Vector{Int}}(length(active_cols_upper))
   list_active_cols_lower = Array{Vector{Int}}(length(active_cols_lower))
   for (v,i) in dict_active_rows; list_active_rows[i] = v   end
   for (v,i) in dict_active_cols_upper; list_active_cols_upper[i] = v end
   for (v,i) in dict_active_cols_lower; list_active_cols_lower[i] = v end

   results = Dict{String,Any}(
       "alpha" => alpha,
       "delta" => delta,
       "epsilon" => epsilon,
       "gamma" => gamma,
       "Minitial" => Minitial,
       "filename" => filename,

       "active_rows" => list_active_rows,
       "active_cols_upper" => list_active_cols_upper,
       "active_cols_lower" => list_active_cols_lower,

       "active_sets" => active_sets,
       "scenario_active_set" => scenario_active_set,
       "scenario_realization" => scenario_realization,

       "observed_active_sets" => observed_active_sets,

       "infeasible_samples" => infeasible_samples,
       "M" => M,
       "W" => W,
       "RoD" => RoD,
       "K_M" => K_M,
       "sigma" => sigma,

       "discovered_rows" => discovered_rows,
       "discovered_sets" => discovered_sets,
       "discovered_cols_upper" => discovered_cols_upper,
       "discovered_cols_lower" => discovered_cols_lower,

       "set_active_rows" => set_active_rows,
       "set_active_cols_lower" => set_active_cols_lower,
       "set_active_cols_upper" => set_active_cols_upper,

       "RoD_Constraint" => RoD_Constraint,

       "iteration_RoD_sets" => iteration_RoD_sets,
       "iteration_RoD_constraints" => iteration_RoD_constraints
   )
   return M, W, RoD, K_M, results

end
# ======================================================================
