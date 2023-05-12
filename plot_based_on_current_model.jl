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

models = []
for chain in 1:TD_parameters["n_chains"]
    model_checkpoint_lists = glob("./models/chain" * string(chain) * "_*")
    if length(model_checkpoint_lists) == 0
        println("ERROR: Couldn't Find Models For Chain" * string(chain))
    else
        # load the newest model
        split1      = split.(model_checkpoint_lists,"%")
        split2      = [split1[i][1] for i in 1:length(split1)]
        split3      = split.(split2,"_")
        split4      = [split3[i][end] for i in 1:length(split3)]
        model_ind   =  findmax(parse.(Float64,split4))[2]
        load_model  = load(model_checkpoint_lists[model_ind])
        for irm in 1:length(model_checkpoint_lists)
            if irm == model_ind
                continue
            end
            rm(sort(model_checkpoint_lists)[irm])
        end
        push!(models,load_model["model_hist"])
        println("--------Chain loaded-------")
    end
end

println("--------Plotting Models-------")
#@time plot_model_hist(models, dataStruct, RayTraces, TD_parameters, TD_parameters["cmax"])
# @time plot_model_hist_weighted(models, dataStruct, RayTraces, TD_parameters, TD_parameters["cmax"])
#@time plot_convergence(TD_parameters)
@time plot_llh(models)



# @time save("model.jld","model",models)

# model_checkpoint_lists = glob("chain*jld")
# [rm(model_checkpoint_lists[i]) for i in 1:length(model_checkpoint_lists)]
