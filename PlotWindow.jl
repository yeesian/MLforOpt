# ==================================================
# Plotting window size for different constants
# ==================================================

# This file is generating plots to assess the affect of choice of constants
# on the window size for our algorithm

# === Inputs/outputs ===

# User defined inputs for generating the window size:
# alpha - fraction of undiscovered mass
# delty - confidence level

# Choices we need to make:
# epsilon - difference between empirically observed rate of discovery and the
# fraction of undicovered mass
# gamma - choice of constant

# Outputs:
# W_M - size of the test window after M observations
# M_min - minimum observation length for the test window

# === Equations ===

# Confidence:
# delta = 2*gamma/(gamma-1)*1/(M_min-1)^(gamma-1)

# Test window length:
# W = 2*gamma/epsilon^2*max(log(M_min), log(M))
# W_min = 2*gamma/epsilon^2*log(M_min)

# Stopping criterion:
# R_MW <= alpha - epsilon

# === To do ===

# 1. Define function to calculate W and M_min for given alpha, delta, epsilon and gamma

# 2. Plot W and M_min for different choices of epsilon and gamma

# =====================================================

using Plots
upscale = 1.3 #8x upscaling in resolution
default(size=(800*upscale,600*upscale)) #Plot canvas size

include("WindowSize.jl")


alpha = 0.1
delta = 0.01
epsilon = vec(collect(0.01:0.01:0.1))
gamma = vec(collect(1.5:0.5:5.0))
W_min = zeros(Float64,length(gamma),length(epsilon))
M_min = zeros(Float64,length(gamma),length(epsilon))
R_max = zeros(Float64,length(gamma),length(epsilon))

for j = 1:length(epsilon)
    for i = 1:length(gamma)
        # Test the window
        (W_min[i,j], M_min[i,j]) = WindowSize(alpha, delta, epsilon[j], gamma[i])
        R_max[i,j] = alpha - epsilon[j]
    end
end


# Plotting W_min
p1 = plot(
    gamma,
    [
        W_min[i,5]
        for i in 1:length(gamma)
    ],
    ylabel="Window Size",
    xlabel="gamma",
    label="epsilon = 0.05"
)

for k in 6:length(epsilon)-1
    plot!(p1,
        gamma,
        [
            W_min[i,k]
            for i in 1:length(gamma)
        ],
        label="epsilon = 0.0$k"
    )
end
plot!(p1,
    gamma,
    [
        W_min[i,length(epsilon)]
        for i in 1:length(gamma)
    ],
    label="epsilon = 0.1"
)

# Plotting M_min
p2 = plot(
    gamma,
    [
        M_min[i,5]
        for i in 1:length(gamma)
    ],
    ylabel="M bar",
    xlabel="gamma",
    label="disregard"
)

# Plotting log(M_bar)
p3 = plot(
    gamma,
    log.(M_min[:,5]),
    ylabel="log(M_bar)",
    xlabel="gamma",
    label="disregard"
)


# Plotting log(M)
x = vec(collect(1.0:1.0:500))
p4 = plot(
    x,
    log.(x),
    ylabel="log(M)",
    xlabel="M",
    label="M_bar"
)

#plot(p1,p2,p3,p4,layout=(4,1))
plot(p1,p2,p3,layout=(3,1))



M = vec(collect(1.0:1.0:100))
p = plot()
for k = 5#:length(epsilon)
    for j = 1:length(gamma)
        plot!(p,
            M,
            [
                2*gamma[j]/epsilon[k]^2*max(log(M_min[j,k]),log(M[i]))
                for i in 1:length(M)
            ],
            label = "gamma = $j"
            )
    end
end

plot(p)

#W = zeros(Float64,length(),1)
#for m = 1:length(M)
#    W = 2*gamma/epsilon^2*max(log(M_min),log(m))
#end
