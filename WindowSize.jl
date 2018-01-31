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

    M_min = 1 + (2*gamma/(delta*(gamma-1)))^(1/(gamma-1));

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

    M_min = 1 + (2*gamma/(delta*(gamma-1)))^(1/(gamma-1));

    W = 2*gamma/epsilon^2*max(log(M_min), log(M));

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

    R_WM = alpha - epsilon

    return R_WM
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

function RunLearningAlgorithm(alpha, delta, epsilon, gamma)

    # Evaluate stopping criterion
    R_max = StoppingCriterion(alpha, delta)

    for m=1:maxM
        # i) Calculate window size for checking
        W = WindowSize(delta, epsilon, gamma, M)

        # ii) Find observed bases for s=1,...,m

        # iii) Find rate of discovery for s=m+1,...,m+W

        # iv) decide whether to terminate
        if R_MW <=R_max
            #end loop
        end
    end

    return M, W, R_MW, K_M
end
# ======================================================================



    R_WM = alpha - epsilon

    return R_WM
end
# ======================================================================
