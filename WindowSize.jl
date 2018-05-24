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

function StoppingCriterion(alpha, delta)

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
    wsamples
    for (i,bus) in ref[:bus]
        bus_std[i] = stddev*bus["Pd"]
        wsamples[i] = rand(Uniform(-3*bus_std[i], 3*bus_std[i]), nsamples)
    end

return wsamples
