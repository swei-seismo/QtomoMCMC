using JLD, HDF5
using Distributed, Statistics, StatsBase
using Distributions
using Glob, Plots, YAML
using DelimitedFiles
using ScatteredInterpolation

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

############## input  ################
file_name = "inp.yml"
add_noise = true
interp4gmt = true
######################################

println("--------Loading Data-------")
@time par = load_par_from_yml(file_name)
lon0, lat0, beta = par["lon0"], par["lat0"], par["beta"]

@time (dataStruct, RayTraces) = load_data_Tonga(par)
make_dir(par, file_name)
println("--------Data Loaded-------")

if isdir(par["base_dir"] * "TongaAttenData") == false
    mkdir(par["base_dir"] * "TongaAttenData")
end

println("--------Loading Synthetic Q-------")
attin = readdlm("checker.in")  
attin = vec(attin) .* 1000 # 1000/Q

lines = readlines(open("Q.grid","r"))
Qgrid = [parse.(Float64, split(line)) for line in lines]
imax,jmax,kmax = Int.(Qgrid[1])
latitudes = [coord[1] for coord in Qgrid[3:end-1]]
longitudes = [coord[2] for coord in Qgrid[3:end-1]]
Qgrid_x, Qgrid_y = lonlat2xy(lon0, lat0, beta, longitudes, latitudes)
depths = parse.(Float64, split(lines[end]))
println("--------Synthetic Q Loaded-------")

##Longitude	 Latitude Depth 1000/Q
println("--------Interpolating Synthetic Q-------")
open("Interpolated_Syn_Atten_Model.txt", "w") do io
    write(io, "Longitude\t Latitude\t Depth\t 1000/Q\n")
end

allX, allY, allDep, allQinv = Float64[], Float64[], Float64[], Float64[]
for k in 1:kmax
    Qinv_lst = []
    for i in 1:imax
        for j in 1:jmax
            node=(i-1)*jmax*kmax+(j-1)*kmax+k
            append!(Qinv_lst,attin[node])
        end
    end
    depth = fill(depths[k], imax*jmax)  
    append!(allX, Qgrid_x)
    append!(allY, Qgrid_y)
    append!(allDep, depth)
    append!(allQinv, Qinv_lst)
end

# the input data points are 2D X-Y plane with increasing depth in the third dimension
# so we used scatter interpolation to interpolate the Q values
points = hcat(allX, allY, allDep)'
itp_Q = ScatteredInterpolation.interpolate(Multiquadratic(),points,allQinv)
#usage: ScatteredInterpolation.evaluate(itp_Q, samples) ([x;y;z]: [allX allY allDep]' or hcat(allX, allY, allDep)' )

if interp4gmt == true
    # map view Q
    for l0 in par["z0"]
        XVec = -93.14142445337606:20.0:986.8585755466239
        YVec = 421.9417419911798:20.0:981.9417419911798

        map_x = repeat(XVec', length(YVec), 1)      # Repeat xVec across rows
        map_y = repeat(YVec, 1, length(XVec))       # Repeat yVec across columns
        #flatten map_x and map_y
        map_x = map_x[:]
        map_y = map_y[:]
        samples = hcat(map_x, map_y, fill(l0, length(map_x)))'
        map_synattin = ScatteredInterpolation.evaluate(itp_Q, samples)
        
        open(par["base_dir"] * "TongaAttenData/Synattin_map_$(l0).txt", "w") do io
            for i in 1:length(map_x)
                lon,lat = xy2lonlat(lon0, lat0, beta, map_x[i], map_y[i])
                write(io, "$(lon)\t$(lat)\t$(map_synattin[i])\n")
            end
        end
    end

    # cross section
    for ixsec in 1:4
        line_fl = "/mnt/home/yurong/Data/MCMC/Tonga/Data/gmt/line$(ixsec).dat"
        line = readdlm(line_fl)
        xsec_x = map(first, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
        xsec_y = map(last, lonlat2xy.(lon0, lat0, beta, line[:,1], line[:,2]))
        dist = line[:,3]
        open(par["base_dir"] * "TongaAttenData/Synattin_xsec_$(ixsec).txt", "w") do io
            for xsec_z in 0.0:20.0:660.0
                samples = hcat(xsec_x, xsec_y, fill(xsec_z, length(xsec_x)))'
                xsec_synattin = ScatteredInterpolation.evaluate(itp_Q, samples)
                for i in 1:length(xsec_synattin)
                    write(io, "$(dist[i])\t$(xsec_z)\t$(xsec_synattin[i])\n")
                end
            end
        end
    end
end

println("--------Synthetic Q Interpolated-------")

println("--------Forward Modeling-------")
if par["coordinates"] == 1
    println("--------Converting Coordinates-------")
    result = lonlat2xy(lon0, lat0, beta, RayTraces["rayX"], RayTraces["rayY"])
    RayTraces["rayX"], RayTraces["rayY"] = result[1], result[2]
    dataStruct["xVec"] = -93.14142445337606:20.0:986.8585755466239
    dataStruct["yVec"] = 421.9417419911798:20.0:981.9417419911798
    println("--------Coordinates Converted-------")
end

(m, n) = size(RayTraces["rayX"])
ptS = zeros(n)
err = zeros(n)
observed_traveltime = zeros(n)
syn_aveatten = zeros(n)

# Based on the distriburion of the t* error, we determin to use a log normal distribution with log_mean -5.6146 and log_std 0.4023
log_mu, log_sigma = -5.6146, 0.4023
@time Threads.@threads for i in 1:n
    samples = hcat(RayTraces["rayX"][:,i], RayTraces["rayY"][:,i], RayTraces["rayZ"][:,i])'
    zeta0 = ScatteredInterpolation.evaluate(itp_Q, samples)
    rayzeta = 0.5 .* (zeta0[1:end-1] + zeta0[2:end])
    rayl = RayTraces["rayL"][:,i]
    index = findfirst(isnan, rayl)

    if index == nothing
        ptS[i] = sum(rayl .* RayTraces["rayU"][:,i] .* (rayzeta ./ 1000))
    else
        ptS[i] = sum(rayl[1:index-1] .* RayTraces["rayU"][1:index-1,i] .* (rayzeta[1:index-1] ./ 1000))
    end

    err[i] = round(rand(LogNormal(log_mu, log_sigma)),digits=6)

    if add_noise == true
        ptS[i] += rand(Normal(0, err[i]))
    end

    observed_traveltime[i] = 1000 * dataStruct["tS"][i] / dataStruct["allaveatten"][i]
    syn_aveatten[i] = 1000*ptS[i] / observed_traveltime[i]
end

tstar_dat = readdlm(par["base_dir"] * "data/p_tstar.dat")

open("synthetic_tstar.dat", "w") do io
    for i = 1:size(tstar_dat, 1)
        write(io, string(tstar_dat[i,1]) *"  "* 
                  string(tstar_dat[i,2]) *"  "* 
                  string(tstar_dat[i,3]) *"  "*
                  string(tstar_dat[i,4]) *"  "* 
                  string(tstar_dat[i,5]) *"  "* 
                  string(round(ptS[i], digits=6)) *"  "*
                  string(round(err[i], digits=6)) *"  "*
                  string(tstar_dat[i,8]) *"  "* 
                  string(round(syn_aveatten[i], digits=2))*"\n")
       
    end
end

println("--------Forward Modeling Done-------")