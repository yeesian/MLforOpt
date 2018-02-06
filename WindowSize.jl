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

    M_min = 1 + (gamma/(delta*(gamma-1)))^(1/(gamma-1));

    W = ceil(2*gamma/epsilon^2*max(log(M_min), log(M)));

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
Run learning algorithm for a given system and set of samples

Inputs:\\
alpha   - maximum unobserved mass \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1


Outputs:\\
"""

function RunLearningAlgorithm(alpha, delta, epsilon, gamma, scenarios, testsize)

    # Evaluate stopping criterion
    R_max = StoppingCriterion(alpha, delta)

    #Find maximum window length (to avoid going too far)
    W_max = WindowSize(delta, epsilon, gamma, length(scenarios.whichbasis[:,1]))

    offset = 0

    for m=1:length(scenarios.whichbasis[:,1])-W_max
        # i) Calculate window size for checking
        W = WindowSize(delta, epsilon, gamma, m)

        # ii) Find observed bases for s=1,...,m
        uniqueM = unique(scenarios.whichbasis[1+offset:m+offset,:],1)
        K_M = size(uniqueM,1)

        # iii) Find rate of discovery for s=m+1,...,m+W
        observedW = scenarios.whichbasis[m+1+offset:m+W+offset,:]
        R_MW = RateOfDiscovery(uniqueM, observedW, W)

        #oldBasis = [sum([observedW[k,:]==uniqueM[j,:] for j in 1:K_M]) for k in 1:size(observedW,1)]

        #R_MW = (W - sum(oldBasis))/W

        # iv) decide whether to terminate
        if R_MW <=R_max

            # Check out-of-sample rate of discovery
            testsize = min(testsize, 50_000 - m - W)
            testsize = Int(ceil(testsize))
            observedTest = scenarios.whichbasis[end-testsize:end,:]
            R_OS = RateOfDiscovery(uniqueM, observedTest, testsize)

            return m, W, R_MW, K_M, R_OS, testsize
        end

        if m==1000*ceil(m/1000)
            println(["M=", m])
        end

    end

    println()
    print("Not enough samples available, algorithm did not terminate!")
    m = 0
    W = 0
    R_MW = 0
    K_M = 0
    R_OS = 0
    testsize = 0
    return m, W, R_MW, K_M, R_OS, testsize
end
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
