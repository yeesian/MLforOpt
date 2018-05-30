# ======================================================================

"""
Run learning algorithm as a streaming algorithm for a given system and pre-computed set of samples

Inputs:\\
alpha   - maximum unobserved mass \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1\\

Outputs:\\
"""

function RunStreamingAlgorithmAC(alpha, delta, epsilon, gamma, Minitial, filename, NLsolver; maxsamples = 1000, tol = 1e-5)

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
    sigma = 0.05
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
    println("rate of discovery: $(numerator / denominator)")

    # Streaming

    M = Minitial;
    ProgressMeter.@showprogress 1 for j = 1:maxsamples
        M += 1;
        Wold = W;
        W = WindowSize(delta, epsilon, gamma, M)

        k_m=length(observed_active_sets)
        println("Iteration: $j M: $M W: $W K_M: $k_m")

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
        println("rate of discovery: $(numerator / denominator)")
        println()

        if RoD <= threshold

            K_M = length(observed_active_sets)

            # Create inverse mapping from the index to the active_rows, etc
            list_active_rows = Array(Vector{Int}, length(active_rows))
            list_active_cols_upper = Array(Vector{Int}, length(active_cols_upper))
            list_active_cols_lower = Array(Vector{Int}, length(active_cols_lower))
            for (v,i) in dict_active_rows; list_active_rows[i] = v end
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

                "infeasible_samples" => infeasible_samples,
                "M" => M,
                "W" => W,
                "RoD" => RoD,
                "K_M" => K_M
            )

            return M, W, RoD, K_M, results
        end
    end

    println("Not enough samples available, algorithm did not terminate!")
    K_M = length(observed_active_sets)
    # Create inverse mapping from the index to the active_rows, etc
    list_active_rows = Array(Vector{Int}, length(active_rows))
    list_active_cols_upper = Array(Vector{Int}, length(active_cols_upper))
    list_active_cols_lower = Array(Vector{Int}, length(active_cols_lower))
    for (v,i) in dict_active_rows; list_active_rows[i] = v end
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

        "infeasible_samples" => infeasible_samples,
        "M" => M,
        "W" => W,
        "RoD" => RoD,
        "K_M" => K_M
    )
    return M, W, RoD, K_M, results

end
# ======================================================================
