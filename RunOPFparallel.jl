
using OPFRecourse, JLD, Gurobi
#using OPFRecourse, JLD, Clp



# Running over OPF Benchmark Cases
f = ARGS[1]

print("Working on $f: ")
data_file = string("pglib-opf/pglib_opf_", f, ".m")
#data_file = string(Pkg.dir(),"/OPFRecourse/test/data/pglib-opf/pglib_opf_", f, ".m")
for i=3 #in 1:5
    @time ref = OPFRecourse.NetworkReference(data_file, σscaling=0.01*i);

    m = OPFRecourse.SingleScenarioOPF(ref, Gurobi.GurobiSolver());
    #m = OPFRecourse.SingleScenarioOPF(ref, Clp.ClpSolver());
    srand(1234)
    scenarios = OPFRecourse.OPFScenarios(ref, m, nsamples = 10);
    JLD.save("results/$(f)$(i).jld", "scenarios", scenarios)
    print("$i ")
end
println()
