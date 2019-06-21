# Plotting DC results across systems

using Plots


k = 14 # number of systems considered

plotly()

nKM  = iteration_K_M[13]
count = collect(1:length(nKM))
plot(count,nKM)

for i=2:k
    nKM  = iteration_K_M[i]
    count = [1:1:length(nKM)]
    if length(nKM) > 1
        plot!(count,nKM)
    end
end

k=12
nKM  = iteration_K_M[k]
count = collect(1:length(nKM))
plot(count,nKM)

nRMW  = iteration_R_MW[k]
count = collect(1:length(nRMW))
plot!(count,nKM[end]*nRMW)

nRMW_diff = abs([1; nRMW[1:end-1]]-nRMW)
nonzeroind = find(nRMW_diff.>0.0001)
count_diff = collect(1:length(nRMW_diff))
bar(sort(nRMW_diff[nonzeroind], rev=true))
