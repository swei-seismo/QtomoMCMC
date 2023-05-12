using JLD, HDF5
using Distributed 
using Glob, Plots
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")

@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
@everywhere include("TD_inversion_function.jl")

lat0 = -23.1000
lon0 = 174.6000
beta = 0.463647609
std_threshold = 5

println("--------Loading Data-------")
@time TD_parameters = define_TDstructrure()
@time (dataStruct, RayTraces) = load_data_Tonga(TD_parameters)
make_dir(TD_parameters)
println("--------Data Loaded-------")

println("--------Loading Models-------")
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
        # for irm in 1:length(model_checkpoint_lists)
        #     if irm == model_ind
        #         continue
        #     end
        #     rm(sort(model_checkpoint_lists)[irm])
        # end
        push_models = load_model["model_hist"]
        push!(models,push_models)
        println("--------Chain loaded-------")
    end
end
println("--------Models Loaded-------")

println("--------Interpolating Models-------")
open("Interpolated_Tonga_Atten_Model.txt", "w") do io
    write(io, "Longitude\t Latitude\t 1000/Q\t Uncertainty\t 1000/Q with Mask\n")
end
for z0 in dataStruct["zVec"]
    m_iz0 = []
    m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))

    for i = 1:length(models)

        for j = 1:length(models[i])

            m  = [ v_nearest(xs, ys, z0,
                models[i][j].xCell, models[i][j].yCell, models[i][j].zCell, models[i][j].zeta)
                for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
            append!(m_iz0,[m])

        end
    end

    model_mean   = mean(m_iz0)
    poststd      = std(m_iz0)
    mask         = ones(size(poststd))
    for i = 1:length(poststd)
        if poststd[i] > std_threshold
            mask[i] = NaN
        end
    end
    mask_model = mask .* model_mean
    open("Interpolated_Tonga_Atten_Model.txt", "a") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                write(io, string(lon) *"\t "* string(lat) *"\t " * string(z0) *"\t "
                    * string(model_mean[i,j]) *"\t "* string(poststd[i,j]) *"\t " 
                    * string(mask_model[i,j]) *"\n")
            end
        end
    end
end
println("--------Models Interpolated-------")


# println("--------Plotting Models-------")
# @time plot_model_hist(models, dataStruct, RayTraces, TD_parameters, TD_parameters["cmax"])
# @time plot_model_hist_weighted(models, dataStruct, RayTraces, TD_parameters, TD_parameters["cmax"])
# @time plot_convergence(TD_parameters)
# @time plot_llh(models)

# @time GMT_plot_model_hist(models, dataStruct, RayTraces, TD_parameters)

# println("--------Models Plotted-------")


# function GMT_plot_model_hist(
#     model_hist, 
#     dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
#     RayTraces::Dict{String,Array{Float64,2}},
#     TD_parameters::Dict{String,Any}, 
#     )
#     cmax = TD_parameters["cmax"]

# end