# ======================================================================
"""
Function to calculate window size
"""
# ======================================================================

function WindowSize(alpha, delta, epsilon, gamma)

    if gamma <= 1
        error("Gamma must be larger than 1!")
    end

    M_min = 1 + (2*gamma/(delta*(gamma-1)))^(1/(gamma-1));

    W_min = 2*gamma/epsilon^2*log(M_min);

    return W_min, M_min
end
