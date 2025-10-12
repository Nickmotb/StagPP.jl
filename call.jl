using StagPP, BenchmarkTools

filename = "test_data/EW_rprof.dat"
# data, header = read_StagYY_timefile!(filename);
data, header = read_StagYY_timefile(filename)
Stag = aggregate_StagData("y20e20", "test_data/EW_parameters.dat", data)
# data = @btime readdlm($filename);