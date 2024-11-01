using JLD, HDF5
using Distributed, Statistics, StatsBase
using Glob, Plots, YAML
using DelimitedFiles

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

Total_mask = 999
Alpha_mask = 5

lat0 = -23.1000
lon0 = 174.6000
beta = 0.463647609
file_name = "inp.yml"

println("--------Loading Data-------")
@time par = load_par_from_yml(file_name)
@time (dataStruct, RayTraces) = load_data_Tonga(par)
make_dir(par, file_name)
println("--------Data Loaded-------")

if isdir(par["base_dir"] * "TongaAttenData") == false
    mkdir(par["base_dir"] * "TongaAttenData")
end

println("--------Loading Models-------")
model_hist = []
for chain in 1:par["n_chains"]
    model_checkpoint_lists = glob(par["base_dir"] * "models/chain" * string(chain) * "_*")
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
        push_models = load_model["model_hist"]
        push!(model_hist,push_models)
        println("--------Chain loaded-------")
    end
end
println("--------Models Loaded-------")

println("--------Interpolating Models-------")

if par["average_style"] == 1
    println("Using Equal Weights")
elseif par["average_style"] == 2
    println("Using Bayesian Model Averaging")
elseif par["average_style"] == 3
    println("Using Inverse Variance Weights")
elseif par["average_style"] == 4
    println("Using Central limit theorem likelihood")
end

# The whole model
open(par["base_dir"] * "TongaAttenData/Interpolated_Tonga_Atten_Model.txt", "w") do io
    write(io, "Longitude\t Latitude\t 1000/Q\t Uncertainty\t 1000/Q with Mask\n")
end
m_chain = []
w_chain = []
m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])), length(vec(dataStruct["zVec"])))
m_all = []
w_all = []
for i = 1:length(model_hist)

    m_i = []
    if par["average_style"] == 2
        w_i = []
    end
    if par["average_style"] == 4
        w_i = []
    end
    for j = 1:length(model_hist[i])

        m  = [ interp_nearest(par,xs, ys, zs, model_hist[i][j])
            for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]), zs in vec(dataStruct["zVec"]) ]
        append!(m_i,[m])
        append!(m_all,[m])
        if par["average_style"] == 2
            append!(w_i, exp(-model_hist[i][j].phi/length(dataStruct["tS"])))
        end
        if par["average_style"] == 4
            append!(w_i, model_hist[i][j].likelihood)
        end
    end
    
    if par["average_style"] == 1
        w_ichain = 1
    elseif par["average_style"] == 2
        w_ichain = mean(w_i)
    elseif par["average_style"] == 3
        w_ichain = 1 / mean(var(m_i))
    elseif par["average_style"] == 4
        w_ichain = mean(w_i)
    end
    append!(m_chain, [mean(m_i)])
    append!(w_chain, w_ichain)
    append!(w_all, ones(length(model_hist[i]))*w_ichain)
end
w_all = Float64.(w_all)
w_chain_norm = w_chain/sum(w_chain)
println("Weights for each chain: ", w_chain_norm)

model_mean = zeros(size(m_all[1]))
model_poststd = zeros(size(m_all[1]))
# model_skewness = zeros(size(m_all[1]))
# model_kurtosis = zeros(size(m_all[1]))
for i = 1:length(model_mean)
    model_mean[i] =  StatsBase.mean([m[i] for m in m_all],Weights(w_all))
    model_poststd[i] =  StatsBase.std([m[i] for m in m_all],Weights(w_all))
    # model_skewness[i] = StatsBase.skewness([m[i] for m in m_all],Weights(w_all))
    # model_kurtosis[i] =  StatsBase.kurtosis([m[i] for m in m_all],Weights(w_all))
end

model_mask         = ones(size(model_poststd))
for i = 1:length(model_poststd)
    if model_poststd[i] > Total_mask 
        model_mask[i] = NaN
    end
end
masked_model = model_mask .* model_mean

open(par["base_dir"] * "TongaAttenData/Interpolated_Tonga_Atten_Model.txt", "a") do io
    for k in 1:length(dataStruct["zVec"])
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if par["coordinates"] == 1
                    write(io, string(dataStruct["xVec"][i]) *"\t "* string(dataStruct["yVec"][j]) *
                        "\t " * string(dataStruct["zVec"][k]) *"\t " * string(model_mean[i,j,k]) *"\t "* 
                        string(model_poststd[i,j,k]) *"\t " * string(masked_model[i,j,k]) *"\n")
                else
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"\t "* string(lat) *"\t " * string(dataStruct["zVec"][k]) *"\t "
                    * string(model_mean[i,j,k]) *"\t "* string(model_poststd[i,j,k]) *"\t " 
                    * string(masked_model[i,j,k]) *"\n")
                end
            end
        end
    end
end



# Map View
for l0 in par["z0"]
    m_chain = []
    m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
    m_all = []
    
    for i = 1:length(model_hist)

        m_i = []
        for j = 1:length(model_hist[i])

            m  = [ interp_nearest(par,xs, ys, l0, model_hist[i][j])
                for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
            append!(m_i,[m])
            append!(m_all,[m])
        end
        append!(m_chain, [mean(m_i)])
    end

    model_mean = zeros(size(m_all[1]))
    model_poststd = zeros(size(m_all[1]))
    # model_skewness = zeros(size(m_all[1]))
    # model_kurtosis = zeros(size(m_all[1]))
    for i = 1:length(model_mean)
        model_mean[i] =  StatsBase.mean([m[i] for m in m_all],Weights(w_all))
        model_poststd[i] =  StatsBase.std([m[i] for m in m_all],Weights(w_all))
        # model_skewness[i] = StatsBase.skewness([m[i] for m in m_all],Weights(w_all))
        # model_kurtosis[i] =  StatsBase.kurtosis([m[i] for m in m_all],Weights(w_all))
    end

    model_mask         = ones(size(model_poststd))
    for i = 1:length(model_poststd)
        if model_poststd[i] > Total_mask 
            model_mask[i] = NaN
        end
    end
    masked_model = model_mask .* model_mean

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Mask_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if par["coordinates"] == 1
                    write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(masked_model[i,j]) *"\n")
                else
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *
                        "  " * string(masked_model[i,j]) *"\n")
                end
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Model_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if par["coordinates"] == 1
                    write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(model_mean[i,j]) *"\n")
                else
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *
                        "  " * string(model_mean[i,j]) *"\n")
                end
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Uncertainty_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if par["coordinates"] == 1
                    write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(model_poststd[i,j]) *"\n")
                else
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *
                        "  " * string(model_poststd[i,j]) *"\n")
                end
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Transparency_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if model_poststd[i,j]<Alpha_mask
                    transparency = 0
                elseif model_poststd[i,j]<Total_mask
                    # transparency = (model_poststd[i,j]-Total_mask)/(Alpha_mask-Total_mask)
                    transparency = 0.9
                else
                    transparency = -1
                end
                if par["coordinates"] == 1
                    write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                            "  " * string(transparency) *"\n")
                else
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *"  " * string(transparency) *"\n")
                end
            end
        end
    end
end

# Interpolate models for points locations in inp/loc.dat for posterior distribution
open(par["base_dir"] * "TongaAttenData/weights.txt", "w") do io
    for i in w_chain_norm
        write(io, string(i) *"\n")
    end
end
# read the locations
loc = readdlm(par["base_dir"] * "loc.dat")
for i in 1:size(loc,1)
    m = []
    iloc = loc[i,:]
    lon, lat, dep, n = iloc
    for i in 1:length(model_hist)
        for j in 1:length(model_hist[i])
            if par["coordinates"] == 1
                m_i = interp_nearest(par,lon, lat, dep, model_hist[i][j])
            else
                (x,y) = lonlat2xy(lon0,lat0,beta,lon,lat)
               m_i = interp_nearest(par,x, y, dep, model_hist[i][j])
            end
            append!(m, m_i)
        end
    end
    open(par["base_dir"] * "TongaAttenData/Posterior_distri_point"*string(i)*".txt", "w") do io
        for i in 1:length(m)
            write(io, string(m[i]) *"\n")
        end
    end
end