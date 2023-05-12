using JLD, HDF5
using Distributed 
using Glob
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
@everywhere include("TD_inversion_function.jl")

@time TD_parameters = define_TDstructrure()
@time (dataStruct, RayTraces) = load_data_Tonga(TD_parameters)
make_dir(TD_parameters)
println("--------Data Loaded-------")

@time models = pmap(x -> TD_inversion_function(TD_parameters, dataStruct, RayTraces, x), 1:TD_parameters["n_chains"])
# @time pmap(x -> TD_inversion_function(TD_parameters, dataStruct, x), 1:TD_parameters["n_chains"])
println("--------Inversion End-------")

# println("--------Plotting Models-------")
# @time plot_model_hist(models, dataStruct, TD_parameters, 30.0)
# @time plot_convergence(TD_parameters)

# @time save("model.jld","model",models)

# model_checkpoint_lists = glob("chain*jld")
# [rm(model_checkpoint_lists[i]) for i in 1:length(model_checkpoint_lists)]


