using JLD, HDF5
using Distributed, Statistics, StatsBase
using Glob, Plots, YAML
using DelimitedFiles
using Interpolations

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

############## INPUTS ##############
Total_mask = 999
Alpha_mask = 5
Add_ray_mask = true
raycount_threshold = 5
file_name = "inp.yml"
###################################




println("--------Loading Data-------")
@time par = load_par_from_yml(file_name)
@time (dataStruct, RayTraces) = load_data_Tonga(par)
lon0, lat0, beta = par["lon0"], par["lat0"], par["beta"]
make_dir(par, file_name)
if isdir(par["base_dir"] * "TongaAttenData") == false
    mkdir(par["base_dir"] * "TongaAttenData")
end
zVec = collect(dataStruct["zVec"])
insert!(zVec, searchsortedfirst(zVec, 50.0), 50.0)
dataStruct["zVec"] = zVec
#####
println("--------Data Loaded-------")


if isfile(joinpath(par["base_dir"], "TongaAttenData/Interpolated_Atten_Model_updated.jld"))
    println("--------Loading Interpolated Model-------")
    model_info = load(joinpath(par["base_dir"], "TongaAttenData/Interpolated_Atten_Model_updated.jld"))
    model_mean = model_info["model_mean"]
    model_poststd = model_info["model_poststd"]
    model_skeness = model_info["model_skeness"]
    model_kurtosis = model_info["model_kurtosis"]
    model_hist = model_info["model_hist"]
    w_chain_norm = model_info["w_chain_norm"]
    w_all = model_info["w_all"]
else
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

    model_mean = sum(m_chain .* w_chain_norm)
    model_poststd = zeros(size(m_all[1]))
    model_skeness = zeros(size(m_all[1]))
    model_kurtosis = zeros(size(m_all[1]))
    @time Threads.@threads for i = 1:length(model_mean)
        # model_mean[i] =  StatsBase.mean([m[i] for m in m_all],Weights(w_all))
        model_poststd[i] =  StatsBase.std([m[i] for m in m_all],Weights(w_all))
        # model_skeness[i] = StatsBase.skewness([m[i] for m in m_all],Weights(w_all))
        # model_kurtosis[i] =  StatsBase.kurtosis([m[i] for m in m_all],Weights(w_all))
    end

    println("--------Saving Interpolated Model-------")
    save(joinpath(par["base_dir"], "TongaAttenData/Interpolated_Atten_Model_updated.jld"), "model_mean", model_mean, "model_poststd", model_poststd, "model_skeness", model_skeness, "model_kurtosis", model_kurtosis, "model_hist", model_hist, "w_chain_norm", w_chain_norm,"w_all", w_all)
    println("--------Interpolated Model Saved-------")
end

println("--------Interpolating Attenuation Model-------")
itp_model = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(dataStruct["zVec"], digits=3)), model_mean, Gridded(Linear()))
itp_model_poststd = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(dataStruct["zVec"], digits=3)), model_poststd, Gridded(Linear()))
itp_model_skeness = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(dataStruct["zVec"], digits=3)), model_skeness, Gridded(Linear()))
itp_model_kurtosis = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(dataStruct["zVec"], digits=3)), model_kurtosis, Gridded(Linear()))
itp_model = extrapolate(itp_model, Flat())
itp_model_poststd = extrapolate(itp_model_poststd, Flat())
itp_model_skeness = extrapolate(itp_model_skeness, Flat())
itp_model_kurtosis = extrapolate(itp_model_kurtosis, Flat())

if Add_ray_mask == true
    if isfile(joinpath(par["base_dir"], "data/RayCount_update.jld"))
        println("--------Loading Interpolated Ray Count-------")
        raycount_info = load(joinpath(par["base_dir"], "data/RayCount_update.jld"))
        raycount = raycount_info["raycount"]
        node = raycount_info["node"]
    else
       println("ERROR: Couldn't Find RayCount_update.jld")
       println("Please run calc_raycount.jl first")
       exit()
    end
    zVec = 20.0:20.0:660.0
    itp_raycount = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(zVec, digits=3)), raycount[:,:,:], Gridded(Linear()))
    itp_raycount = extrapolate(itp_raycount, Flat())
end

# map view
for l0 in par["z0"]
    XVec = dataStruct["xVec"]
    YVec = dataStruct["yVec"]

    map_x = repeat(XVec', length(YVec), 1)      # Repeat xVec across rows
    map_y = repeat(YVec, 1, length(XVec))       # Repeat yVec across columns
    map_model = itp_model.(map_x, map_y, l0)
    map_model_poststd = itp_model_poststd.(map_x, map_y, l0)
    map_model_skeness = itp_model_skeness.(map_x, map_y, l0)
    map_model_kurtosis = itp_model_kurtosis.(map_x, map_y, l0)
    map_model_mask = ifelse.(map_model_poststd .> Total_mask, NaN, map_model)
    if Add_ray_mask == true
        map_model_alpha = zeros(size(map_model))
        for i in 1:length(map_model_alpha)
            if map_model_poststd[i] > Alpha_mask || itp_raycount(map_x[i], map_y[i], l0) < raycount_threshold
                map_model_alpha[i] = 0.9 
            end
        end
    else
        map_model_alpha = ifelse.(map_model_poststd .> Alpha_mask, 0.9, 0.0)
    end
    open(par["base_dir"] * "TongaAttenData/Atten_map_$(l0).txt", "w") do io
        for i in 1:size(map_model, 1), j in 1:size(map_model, 2)
            if par["coordinates"] == 1
                lon, lat = map_x[i,j], map_y[i,j]
            else
                (lon, lat) = xy2lonlat(lon0, lat0, beta, map_x[i,j], map_y[i,j])
            end
            write(io, "$(lon)\t$(lat)\t$(map_model[i,j])\t$(map_model_poststd[i,j])\t$(map_model_mask[i,j])\t$(map_model_alpha[i,j])\t$(map_model_skeness[i,j])\t$(map_model_kurtosis[i,j])\n")
        end
    end
    println("--------Map $(l0) Saved-------")
    println("lon:\n", dataStruct["xVec"])
    println("lat:\n", dataStruct["yVec"])
    println("depth:\n", dataStruct["zVec"])
end

# cross section
for ixsec in 1:4
    line_fl = "/mnt/home/yurong/data/Tonga/MCMC/gmt/line$(ixsec).dat"
    line = readdlm(line_fl)
    if par["coordinates"] == 1
        xsec_x = line[:,1]
        xsec_x = ifelse.(xsec_x .< 0, xsec_x .+ 360, xsec_x)
        xsec_y = line[:,2]
    else
        xsec_x = map(first, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
        xsec_y = map(last, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
    end
    
    dist = line[:,3]
    open(par["base_dir"] * "TongaAttenData/Atten_xsec_$(ixsec).txt", "w") do io
        for xsec_z in 20.0:20.0:660.0
            xsec_model = itp_model.(xsec_x, xsec_y, xsec_z)
            xsec_model_poststd = itp_model_poststd.(xsec_x, xsec_y, xsec_z)
            xsec_model_mask = ifelse.(xsec_model_poststd .> Total_mask, NaN, xsec_model)
            if Add_ray_mask == true
                xsec_model_alpha = zeros(size(xsec_model))
                for i in 1:length(xsec_model_alpha)
                   if xsec_model_poststd[i] > Alpha_mask || itp_raycount(xsec_x[i], xsec_y[i], xsec_z) < raycount_threshold
                        xsec_model_alpha[i] = 0.9 
                   end
                end
            else
                xsec_model_alpha = ifelse.(xsec_model_poststd .> Alpha_mask, 0.9, 0.0)
            end
            xsec_model_skeness = itp_model_skeness.(xsec_x, xsec_y, xsec_z)
            xsec_model_kurtosis = itp_model_kurtosis.(xsec_x, xsec_y, xsec_z)
            for i in 1:length(xsec_model)
                write(io, "$(dist[i])\t$(xsec_z)\t$(xsec_model[i])\t$(xsec_model_poststd[i])\t$(xsec_model_mask[i])\t$(xsec_model_alpha[i])\t$(xsec_model_skeness[i])\t$(xsec_model_kurtosis[i])\n")
            end
        end
    end
end

println("--------Calculating Posterior distribution of selected locations-------")
# Interpolate models for points locations in inp/loc.dat for posterior distribution
open(par["base_dir"] * "TongaAttenData/weights.txt", "w") do io
    for i in w_chain_norm
        write(io, string(i) *"\n")
    end
end
# read the locations
loc = readdlm("/mnt/home/yurong/data/Tonga/MCMC/gmt/loc.dat")
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

println("--------Done-------")