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
gr()

upscale = 0.5#1.3 #8x upscaling in resolution
default(size=(800*2*upscale,600*upscale)) #Plot canvas size

include("WindowSize.jl")


alpha = 0.1
delta = 0.01
epsilon = vec(collect(0.03: 0.01 : 0.07))
linestyle = (:solid, :dash, :dot, :dashdot, :dashdotdot)
gamma = vec(collect(1.5 : 0.1 : 6))
M = vec(collect(1.0:1.0:800))
W_min = zeros(Float64,length(gamma),length(epsilon))
M_min = zeros(Float64,length(gamma),length(epsilon))
#R_max = zeros(Float64,length(gamma),length(epsilon))


# Calculating the initial window size for different epsilon, gamma
for j = 1:length(epsilon)
    for i = 1:length(gamma)
        # Test the window
        (W_min[i,j], M_min[i,j]) = MinWindowSize(delta, epsilon[j], gamma[i])
        #R_max[i,j] = alpha - epsilon[j]
    end
end

# Plotting the initial window size
start_eps = 1
p1 = plot(
    gamma,
    [
        W_min[i,start_eps]
        for i in 1:length(gamma)
    ],
    ylabel="Initial Window Size W(1)",
    xlabel="Hyperparameter \\gamma",
    legend = :topright,
    label="\\epsilon = $(epsilon[start_eps])",
    ylims = (0,30000)

)

for k in (start_eps+1):length(epsilon)
    plot!(p1,
        gamma,
        [
            W_min[i,k]
            for i in 1:length(gamma)
        ],
        label="\\epsilon = $(epsilon[k])",
        linestyle = linestyle[k],
        ylims = (0,35000)
#        xlims = (1.5,7)
    )
end

fn = "InitialWindowSize-delta0p01"
savefig(fn)





# Plotting M_min

# p2 = plot(
#     gamma,
#     [
#         M_min[i,5]
#         for i in 1:length(gamma)
#     ],
#     ylabel="M bar",
#     xlabel="gamma",
#     label="disregard"
# )
#
# # Plotting log(M_bar)
# p3 = plot(
#     gamma,
#     log.(M_min[:,5]),
#     ylabel="log(M_bar)",
#     xlabel="gamma",
#     label="disregard"
# )
#
#
# # Plotting log(M)
# x = vec(collect(1.0:1.0:500))
# p4 = plot(
#     x,
#     log.(x),
#     ylabel="log(M)",
#     xlabel="M",
#     label="M_bar"
# )
#
# plot(p1,p2,p3,p4,layout=(4,1))
# #plot(p1,p2,p3,layout=(3,1))
#
# fn = "AllPlots-delta0p01.png"
# savefig(fn)



gamma_new = vec(collect(1.5 : 0.5 : 3.5))
epsilon_new = 0.04
p = plot()

for k = 1:length(gamma_new)
    j = length(gamma_new)+1-k
    plot!(p,
        M,
        [
            WindowSize(delta, epsilon_new, gamma_new[j], M[i])
            for i in 1:length(M)
        ],
        label = "\\gamma = $(gamma_new[j])",
        ylabel="Window Size W(M)",
        xlabel="Iteration number M",
        ylims = (0,35000),
        #legend=:bottomright
        linestyle = linestyle[k],
        legend=:bottomright
        )
end

plot(p)

fn = "PerIteration-delta0p01-epsilon0p04.png"
savefig(fn)

#W = zeros(Float64,length(),1)
#for m = 1:length(M)
#    W = 2*gamma/epsilon^2*max(log(M_min),log(m))
#end

plot(p1,p,layout=(1,2))

fn = "Paper-plots-delta0p01-epsilon0p04.png"
savefig(fn)
