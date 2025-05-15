using JLD, HDF5
using Distributed
using Glob, Plots, YAML
using Interpolations
using DelimitedFiles

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

##########    INPUTS    ##########
file_name = "TongaLau_Spherical.yml"
raycount_threshold = 5
##################################

par = load_par_from_yml(file_name)
lon0, lat0, beta = par["lon0"], par["lat0"], par["beta"]

if isfile(joinpath(par["base_dir"], "data/RayCount_update.jld"))
    println("--------Loading Ray Count-------")
    raycount_info = load(joinpath(par["base_dir"], "data/RayCount_update.jld"))
    node = raycount_info["node"]
    par = raycount_info["par"]
    dataStruct = raycount_info["dataStruct"]
    raycount = raycount_info["raycount"]
else
    println("--------Loading Data-------")
    (dataStruct, RayTraces) = load_data_Tonga(par)
    make_dir(par, file_name)
    println("--------Data Loaded-------")

    if isdir(par["base_dir"] * "TongaAttenData") == false
        mkdir(par["base_dir"] * "TongaAttenData")
    end

    if par["coordinates"] == 1
        println("--------Converting Coordinates-------")
        result = lonlat2xy(lon0, lat0, beta, RayTraces["rayX"], RayTraces["rayY"])
        RayTraces["rayX"], RayTraces["rayY"] = result[1], result[2]
    end


    raycount = zeros(length(dataStruct["xVec"]), length(dataStruct["yVec"]), length(dataStruct["zVec"]))
    node = fill((0.0, 0.0, 0.0), length(dataStruct["xVec"]), length(dataStruct["yVec"]), length(dataStruct["zVec"]))

    for i in 1:length(dataStruct["xVec"]), j in 1:length(dataStruct["yVec"]), k in 1:length(dataStruct["zVec"])
        if par["coordinates"] == 1
            (x,y)= lonlat2xy(lon0, lat0, beta, dataStruct["xVec"][i], dataStruct["yVec"][j])
        else
            (x,y) = (dataStruct["xVec"][i], dataStruct["yVec"][j])
        end
        node[i,j,k] = (x, y, dataStruct["zVec"][k])
    end


    println("--------Calculating Ray Count-------")
    @time Threads.@threads for ix in 1:length(raycount)
        raycount[ix] = calc_raycount(ix, node, par, dataStruct, RayTraces)
    end

    println("--------Ray Count Calculated-------")


    println("--------Saving Ray Count-------")
    save(joinpath(par["base_dir"], "data/RayCount_update.jld"), "raycount", raycount, "node", node, "par", par, "dataStruct", dataStruct)
    println("--------Ray Count Saved-------")
end

println("--------Interpolating Ray Count-------")
itp_raycount = interpolate((round.(dataStruct["xVec"], digits=3), round.(dataStruct["yVec"], digits=3), round.(dataStruct["zVec"], digits=3)), raycount[:,:,:], Gridded(Linear()))
itp_raycount = extrapolate(itp_raycount, Flat())

# map view
for l0 in par["z0"]
    XVec = dataStruct["xVec"]
    YVec = dataStruct["yVec"]

    map_x = repeat(XVec', length(YVec), 1)      # Repeat xVec across rows
    map_y = repeat(YVec, 1, length(XVec))       # Repeat yVec across columns
    map_raycount = itp_raycount.(map_x, map_y, l0)
    map_raymask = ifelse.(map_raycount .> raycount_threshold, 0.0, 0.9)
    open(par["base_dir"] * "TongaAttenData/RayCount_map_$(l0).txt", "w") do io
        for i in 1:size(map_raycount, 1), j in 1:size(map_raycount, 2)
            if par["coordinates"] == 1
                (lon, lat) = (map_x[i,j], map_y[i,j])
            else
                lon,lat = xy2lonlat(lon0, lat0, beta, map_x[i,j], map_y[i,j])
            end
            write(io, "$(lon)\t$(lat)\t$(map_raycount[i,j])\t$(map_raymask[i,j])\n")
        end
    end
end

# cross section
for ixsec in 1:4
    line_fl = "/mnt/home/yurong/data/Tonga/MCMC/gmt/line$(ixsec).dat"
    line = readdlm(line_fl)
    if par["coordinates"] == 1
        xsec_x = line[:,1]
        xsec_x = ifelse.(xsec_x .< 0, xsec_x .+ 360.0, xsec_x)
        xsec_y = line[:,2]
    else
        xsec_x = map(first, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
        xsec_y = map(last, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
    end
    dist = line[:,3]
    open(par["base_dir"] * "TongaAttenData/RayCount_xsec_$(ixsec).txt", "w") do io
        for xsec_z in 0.0:20.0:660.0
            xsec_raycount = itp_raycount.(xsec_x, xsec_y, xsec_z)
            xsec_raymask = ifelse.(xsec_raycount .> raycount_threshold, 0.0, 0.9)
            for i in 1:length(xsec_raycount)
                write(io, "$(dist[i])\t$(xsec_z)\t$(xsec_raycount[i])\t$(xsec_raymask[i])\n")
            end
        end
    end
end