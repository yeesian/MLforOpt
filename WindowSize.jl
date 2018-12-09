# ======================================================================
"""
Function to calculate minimum window size

Inputs: \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1

Outputs:\\
W_min   - minimium window size\\
M_min   - number of samples for which we have the minimum window size

"""

function MinWindowSize(delta, epsilon, gamma)

    if gamma <= 1
        error("Gamma must be larger than 1!")
    end

    M_min = 1 + (gamma/(delta*(gamma-1)))^(1/(gamma-1));

    W_min = 2*gamma/epsilon^2*log(M_min);

    return W_min, M_min
end
# ======================================================================



# ======================================================================
"""
Function to calculate window size after M samples

Inputs: \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1 \\
M       - number of samples seen so far

Outputs:\\
W       - window size\\
"""

function WindowSize(delta, epsilon, gamma, M)

    if gamma <= 1
        error("Gamma must be larger than 1!")
    end

    W_min, M_min = MinWindowSize(delta, epsilon, gamma)

    W = ceil(2*gamma/epsilon^2*max(log(M_min), log(M)))

    W = Int(W)

    return W
end
# ======================================================================



# ======================================================================
"""
Calculate rate of discovery where the algorithm terminates

Inputs:\\
alpha - maximum probbaility of unobserved set\\
epsilon - difference between rate of discovery and probability of unobserved set

Outputs:\\
R_MW - rate of discovery where algorithm terminates
"""

function StoppingCriterion(alpha, epsilon)

    R_MW = alpha - epsilon

    return R_MW
end
# ======================================================================



# ======================================================================
"""
Determine the rate of discovery of new bases, given a set of already observed
bases and a set of new scenarios.whichbasis

Inputs:\\
uniqueM     unique bases already undiscovered\\
observedB   observed scenarios\\
testsize    length of the test window

Outputs:\\
RoD         Rate of Discovery of new bases in the test window

"""
function RateOfDiscovery(uniqueM, observedB, testsize)
        oldBasis = [sum([observedB[k,:]==uniqueM[j,:] for j in 1:size(uniqueM,1)]) for k in 1:size(observedB,1)]
        RoD = (testsize - sum(oldBasis))/testsize
    return RoD
end
# ======================================================================


# ======================================================================
"""
Find the number of feasible and infeasible scenarios
"""
function CheckFeasibility(scenarios, nsamples)

    n_feasible = size(scenarios.whichbasis)[1]
    n_infeasible = nsamples - size(scenarios.whichbasis)[1]

    return n_feasible, n_infeasible
end
# ======================================================================

"""
Out of sample test
"""
function OutOfSample(testsize, m, W, scenarios, uniqueM)
    # Check out-of-sample rate of discovery
    testsize = min(testsize, 50_000 - m - W)
    testsize = Int(ceil(testsize))
    observedTest = scenarios.whichbasis[end-testsize:end,:]
    R_OS = RateOfDiscovery(uniqueM, observedTest, testsize)

    return R_OS, testsize
end
# ======================================================================

"""
Generate independent uniform samples based on load
"""
function UniformSamples(stddev, data)
    wsamples = []
    for (i,bus) in ref[:bus]
        bus_std[i] = stddev*bus["Pd"]
        wsamples[i] = rand(Uniform(-3*bus_std[i], 3*bus_std[i]), nsamples)
    end

    return wsamples
end


# ======================================================================

"""
Add new constraints to the set of possible active constraints

Inputs: \\
collections_of_active_constraints   Collections of active constraints for which we want to add missing elements (if any)\\
set_of_active_indices               Set of already discovered active constraints\\

"""

function add_new_constraints(collections_of_active_constraints, set_of_active_indices)
    #println("Adding missing constraints")
    for ind_vec = 1:length(collections_of_active_constraints)
        #println("Collection number $ind_vec")
        for ind_con=1:length(collections_of_active_constraints[ind_vec])
            if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
                #println("- Discovered new constraint at element $ind_con . Adding this constraint")
                push!(set_of_active_indices, collections_of_active_constraints[ind_vec][ind_con])
            #else
            #    println("Element already discovered")
            end

        end
    end
    return set_of_active_indices
end

# ======================================================================

"""
Check discovery status for all of the collections of constraints

Inputs: \\
collections_of_active_constraints   All collections of active constraints\\
set_of_active_indices               Set of already discovered active constraints\\
discovery_status                    Whether or not a particular set of active constraints has been discovered\\

"""
function check_discovery(collections_of_active_constraints, set_of_active_indices, discovery_status)

    # Checking that the size of the collection has similar dimensions as the discovery status
    @assert length(collections_of_active_constraints)==length(discovery_status)

    # Updating discovery status
    for ind_vec = 1:length(collections_of_active_constraints)
        temporary_var = discovery_status[ind_vec]
        #println("Collection number $ind_vec")
        #println("- Discovery status: $temporary_var")
        if discovery_status[ind_vec] == false
            discovered = true
            for ind_con=1:length(collections_of_active_constraints[ind_vec])
                if !in(collections_of_active_constraints[ind_vec][ind_con], set_of_active_indices)
                    #println("- Discovered new constraint at element $ind_con")
                    discovered = false
                end
            end
            discovery_status[ind_vec]=discovered
            #println("- Are all constraints discovered? $discovered")
            # if discovered == true
            #     println("- New discovery status for collection number $ind_vec: $discovered")
            # end
        end
    end
    return discovery_status
end
