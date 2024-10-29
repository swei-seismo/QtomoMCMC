using JLD, HDF5
using Distributed 
using Glob, Plots, YAML

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

std_threshold = 5
transparency_threshold = 13

file_name = "inp.yml"

println("--------Loading Data-------")
@time par = load_par_from_yml(file_name)
@time (dataStruct, RayTraces) = load_data_Tonga(par)
make_dir(par, file_name)
println("--------Data Loaded-------")

if isdir(par["base_dir"] * "TongaAttenData") == false
    mkdir(par["base_dir"] * "TongaAttenData")
end

println("--------Plotting Convergence-------")
@time plot_convergence(par, 1)
@time plot_convergence(par, 100)
@time plot_convergence(par, 10000)