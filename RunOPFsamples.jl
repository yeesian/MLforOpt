
using OPFRecourse, JLD, Clp



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
         "case2736sp_k",
         "case2737sop_k",
         "case2746wop_k",
         "case2746wp_k",
         "case2848_rte",
         "case2868_rte",
         "case2869_pegase",
         "case2869_pegase",
         "case3012wp_k",
         "case3120sp_k",
         "case3375wp_k",
         "case6468_rte",
         "case6495_rte",
         "case6515_rte",
         "case9241_pegase",
         "case13659_pegase"
    ]
    print("Working on $f: ")
    data_file = string("pglib-opf/pglib_opf_", f, ".m")
    #data_file = string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/pglib_opf_", f, ".m")
    for i in 1:5
        @time ref = OPFRecourse.NetworkReference(data_file, σscaling=0.01*i);

        m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
        #m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
        srand(1234)
        scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 10);
        JLD.save("results/$(f)$(i).jld", "scenarios", scenarios)
        print("$i ")
    end
    println()
end
