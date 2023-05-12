using JLD, HDF5
using Distributed 
using Glob

@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")

@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")

@time TD_parameters = define_TDstructrure()
@time dataStruct = load_data_Tonga(TD_parameters)
make_dir(TD_parameters)
println("--------Data Loaded-------")

for i in 1:TD_parameters["n_chains"]
    model_checkpoint_lists = glob("./models/chain" * string(i) * "_*")
    split1      = split.(model_checkpoint_lists,"%")
    split2      = [split1[i][1] for i in 1:length(split1)]
    split3      = split.(split2,"_")
    split4      = [split3[i][end] for i in 1:length(split3)]
    model_ind   =  findmax(parse.(Float64,split4))[2]
    chain       = load(model_checkpoint_lists[model_ind])

    # chain = load("./models/chain"*string(i)*"_iter100000_100.0%.jld")
    println("--------Model Loaded-------")
    n = 5
    ind_lst = Float64.(round.(rand(n)*length(chain["model_hist"])))
    for iter in ind_lst
        println(i," ",iter)
        plot_voronoi(chain["model_hist"][Int64(iter)],dataStruct,TD_parameters,i,iter)
    end
end