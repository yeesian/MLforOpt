


discovered_rows       = zeros(results["active_rows"],Bool)
discovered_cols_upper = zeros(results["active_cols_upper"],Bool)
discovered_cols_lower = zeros(results["active_cols_lower"],Bool)

set_active_rows       = Set{Int}()
set_active_cols_upper = Set{Int}()
set_active_cols_lower = Set{Int}()

discovered_sets = Set{Tuple{Int,Int,Int}}()

# Initialization: Add constraints from the first M active sets

println("=== Identifying known active sets and adding their constraints ===")
known_active_sets = unique(results["scenario_active_set"][1:M])

println("ACTIVE ROWS")
known_active_rows       = unique([results["active_sets"][known_active_sets][i][1] for i in known_active_sets])
set_active_rows         = add_new_constraints(results["active_rows"][known_active_rows], set_active_rows)
discovered_rows[known_active_rows] = true

println("ACTIVE COLS UPPER")
known_active_cols_upper = unique([results["active_sets"][known_active_sets][i][2] for i in known_active_sets])
set_active_cols_upper   = add_new_constraints(results["active_cols_upper"][known_active_cols_upper], set_active_cols_upper)
discovered_cols_upper[known_active_cols_upper] = true

println("ACTIVE COLS LOWER")
known_active_cols_lower = unique([results["active_sets"][known_active_sets][i][3] for i in known_active_sets])
set_active_cols_lower   = add_new_constraints(results["active_cols_lower"][known_active_cols_lower], set_active_cols_lower)
discovered_cols_lower[known_active_cols_lower] = true

# Check if other sets are covered by the discovered ones
println("=== Update discovery status for all collections ===")

println("ACTIVE ROWS")
set_active_rows, discovered_rows = check_discovery(results["active_rows"],
                                                   set_active_rows,
                                                   discovered_rows)

println("ACTIVE COLS UPPER")
set_active_cols_upper, discovered_cols_upper = check_discovery(results["active_cols_upper"],
                                                               set_active_cols_upper,
                                                               discovered_cols_upper)

println("ACTIVE COLS LOWER")
set_active_cols_lower, discovered_cols_lower = check_discovery(results["active_cols_lower"],
                                                                set_active_cols_lower,
                                                                discovered_cols_lower)




function check_discovery(collections_of_active_constraints, set_of_active_indices, discovery_status)
    # == collections_of_active_constraints:
    # the full set of all different collections of active constraints
    # == set_of_active_indices:
    # set of already discovered active constraints
    # == discovery_status:
    # sign of whether or not a particular set of active constraints has been discovered

    # Checking that the size of the collection has similar dimensions as the discovery status
    @assert length(collections_of_active_constraints)==length(discovery_status)

    # Updating discovery status
    for ind_vec = 1:length(collections_of_active_constraints)
        temporary_var = discovery_status[ind_vec]
        println("Collection number $ind_vec")
        println("- Discovery status: $temporary_var")
        if discovery_status[ind_vec] == false
            discovered = true
            for ind_con=1:length(collections_of_active_constraints[ind_vec])
                if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
                    println("- Discovered new constraint at element $ind_con")
                    discovered = false
                end
            end
            discovery_status[ind_vec]=discovered
            println("- Are all constraints discovered? $discovered")
            if discovered == true
                println("- New discovery status for $ind_vec: $discovered")
            end
        end
    end
    return set_of_active_indices, discovery_status
end

function add_new_constraints(collections_of_active_constraints, set_of_active_indices)
    # == collections_of_active_constraints:
    # collections of active constraints for which we want to add missing elements (if any)
    # == set_of_active_indices:
    # set of already discovered active constraints
    println("Adding missing constraints")
    for ind_vec = 1:length(collections_of_active_constraints)
        println("Collection number $ind_vec")
        for ind_con=1:length(collections_of_active_constraints[ind_vec])
            if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
                println("- Discovered new constraint at element $ind_con . Adding this constraint")
                push!(set_of_active_indices, collections_of_active_constraints[ind_vec][ind_con])
            end
        end
    end
    return set_of_active_indices
end



# for ind_vec = 1:length(results["active_cols_upper"])
#     println("active_cols_upper number $ind_vec")
#     discovered = true
#     for ind_con=1:length(results["active_cols_upper"][ind_vec])
#         if !in(results["active_cols_upper"][ind_vec][ind_con], set_active_cols_upper)
#             println("Discovered new constraint at element $ind_con")
#             discovered = false
#             push!(set_active_cols_upper, results["active_cols_upper"][ind_vec][ind_con])
#         end
#     end
#     if discovered == true
#         discovered_cols_upper[ind_vec]=true
#     end
#     println("What we have discovered so far: $discovered_cols_upper")
# end
#
#
# function check_and_add_new_constraints(collections_of_active_constraints, set_of_active_indices, discovery_status)
#     for ind_vec = 1:length(collections_of_active_constraints)
#         println("Collection number $ind_vec")
#         discovered = true
#         for ind_con=1:length(collections_of_active_constraints[ind_vec])
#             if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
#                 println("Discovered new constraint at element $ind_con")
#                 discovered = false
#                 push!(set_of_active_indices, collections_of_active_constraints[ind_vec][ind_con])
#             end
#         end
#         discovery_status[ind_vec]=discovered
#         println("What we have discovered so far: $discovered_cols_upper")
#     end
#     return set_of_active_indices, discovery_status
# end

# function check_discovery(collections_of_active_constraints, set_of_active_indices, discovery_status)
#     for ind_vec = 1:length(collections_of_active_constraints)
#         println("Collection number $ind_vec")
#         discovered = true
#         for ind_con=1:length(collections_of_active_constraints[ind_vec])
#             if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
#                 println("Discovered new constraint at element $ind_con")
#                 discovered = false
#             end
#         end
#         discovery_status[ind_vec]=discovered
#         println("What we have discovered so far: $discovered_cols_upper")
#     end
#     return set_of_active_indices, discovery_status
# end

# push!(set_active_cols_upper, results["active_cols_upper"][1][1])
#
# empty!(set_active_cols_upper)
#
# for i=1:length(results["active_cols_upper"][2])
#     push!(set_active_cols_upper, results["active_cols_upper"][2][i])
# end
