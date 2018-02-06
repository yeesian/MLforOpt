
using JuMP, JLD, OPFRecourse, Clp, Distributions
#using OPFRecourse, JLD, Clp



# Running over OPF Benchmark Cases
for f in [
        "case3_lmbd",
        "case5_pjm",
        "case14_ieee",
        "case24_ieee_rts",
        #"case30_as",
        #"case30_fsr",
        "case30_ieee",
        "case39_epri",
        "case57_ieee",
        "case73_ieee_rts",
        # "case89_pegase", # (infeasible)
        "case118_ieee",
        "case162_ieee_dtc",
        "case200_pserc",
        "case240_pserc",
        "case300_ieee",
        # "case1354_pegase", # (infeasible)
         "case1888_rte",
         "case1951_rte",
         # "case2383wp_k" # (infeasible)
         #"case2736sp_k",
         "case2737sop_k", #seems fine
         #"case2746wop_k", # (infeasible)
         #"case2746wp_k", # (infeasible)
         "case2848_rte", #seems fine
         "case2868_rte", #seems fine
         "case2869_pegase", #seems fine (but seems to take longer)
         "case2869_pegase", #seems fine (but also longer)
         #"case3012wp_k", # (infeasible)
         #"case3120sp_k", # (infeasible)
         #"case3375wp_k", # (singluar)
         "case6468_rte", # seems fine (abou 1 min 30s for 10 samples)
         "case6495_rte", # seems fine (abou 1 min 37s for 10 samples)
         "case6515_rte", # seems fine (abou 1 min 40s for 10 samples)
         #"case9241_pegase", # infeasible
         "case13659_pegase"
    ]
    print("Working on $f: ")
    data_file = string("pglib-opf/pglib_opf_", f, ".m")
    #data_file = string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/pglib_opf_", f, ".m")
    for i=3 #in 1:5
        @time ref = OPFRecourse.NetworkReference(data_file, σscaling=0.01*i);

        #m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
        @time m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
        @time JuMP.solve(m.model)

        # Generate new uncertainties based on uniform distribution
        srand(1234)
        nsamples = 10
        ωsamples = zeros(ref.nbus, nsamples)
        nonzeroindices = (1:length(ref.stdω))[ref.stdω .> 1e-5]
        for b in nonzeroindices
            ωsamples[b,:] = rand(Uniform(-3*ref.stdω[b], 3*ref.stdω[b]), nsamples)
        end
        # solve OPF for uniform distribution realizations
        scenarios = OPFRecourse.OPFScenarios(ref, m, ωsamples);
        JLD.save("results_uniform/$(f)$(i).jld", "scenarios", scenarios)
        print("$i ")
    end
    println()
end
