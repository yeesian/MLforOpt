# ======================================================================
"""
Run learning algorithm as a streaming algorithm for a given system and pre-computed set of samples

Inputs:\\
alpha   - maximum unobserved mass \\
delta   - confidence level \\
epsilon - difference between rate of discovery and probability of unobserved set\\
gamma   - constant > 1\\
scenarios - pre-computed set of scenarios\\
testsize  - number of scenarios used for out-of-sample test\\

Outputs:\\
"""

function RunStreamingAlgorithmAC(alpha, delta, epsilon, gamma, filename)

    # Evaluate stopping criterion
    R_max = StoppingCriterion(alpha, delta)

    #Find maximum window length (since the number of samples is limited)
    #W_max = WindowSize(delta, epsilon, gamma, length(scenarios.whichbasis[:,1]))

    # Initialize
    m = Mmin
    W = WindowSize(delta, epsilon, gamma, m)

    # Posting model with NL parameters for omega
    #UPDATE!!!!
    network_data = PowerModels.parse_file(filename)
    m = Model(solver = NLsolver)

    # Run OPF to get active sets
    for i = 1:m+W
        # 1. Generate sample
        # 2. Fix NL parameters
        # 3. Get active set
        active_set[i] = find_active_set(jm, const_refs, var_refs, tol)
        # 4. Create a simplified marker for the active set (similar to scenarios.whichbasis)
    end

    # WRITE SOME FUNCTION FOR THIS

    # find unique basis in M
    uniqueM = unique(scenarios.whichbasis[1:m,:],1)
    K_M = size(uniqueM,1)
    # determine which samples in W correspond to already observed basis
    basisW = scenarios.whichbasis[m+1:m+W,:]
    observed = [sum([basisW[k,:]==uniqueM[j,:] for j in 1:size(uniqueM,1)]) for k in 1:size(basisW,1)]
    # calculate rate of discovery
    R_MW = (W - sum(observed))/W

    for m=Mmin+1:length(scenarios.whichbasis[:,1])-W_max
        # Calculate window size
        Wold = W
        W = WindowSize(delta, epsilon, gamma, m)

        # Calculate number of new samples
        Nnew = W - Wold + 1

        for i = 1:Nnew
            # calculate new active sets
        end

        # Update the set of unique basis
        uniqueM = unique([uniqueM; scenarios.whichbasis[m,:]'],1)
        K_Mold = K_M
        K_M = size(uniqueM,1)

        # Update the basis which have been observed in window W
        if K_Mold == K_M # No new basis observed
            # Check if new samples have been observed
            basisNew = scenarios.whichbasis[Wold+1:Wold+Nnew,:]
            observedNew = [sum([basisNew[k,:]==uniqueM[j,:] for j in 1:size(uniqueM,1)]) for k in 1:size(basisNew,1)]
            # Construct new observed vector
            observed = [observed[2:end]; observedNew]
        else # New basis observed
            # Check if new samples have been observed
            basisNew = scenarios.whichbasis[Wold+1:Wold+Nnew,:]
            observedNew = [sum([basisNew[k,:]==uniqueM[j,:] for j in 1:size(uniqueM,1)]) for k in 1:size(basisNew,1)]
            # Update the observed vector
            for i=1:length(observed)
                if observed[i] == 0 # Only check if previously unobserved
                    observed[i] = sum([scenarios.whichbasis[m-1+i,:]==uniqueM[j,:] for j in 1:size(uniqueM,1)])
                end
            end
            # Construct new observed vector
            observed = [observed[2:end]; observedNew]
        end

        #observedW = scenarios.whichbasis[m+1+offset:m+W+offset,:]
        R_MW = (W - sum(observed))/W

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
