using Random, Distributions, Distributed, Base.Threads
using Plots, Glob, StatsBase, ColorSchemes
using HDF5, JLD


function make_dir(TD_parameters)
    FigDirLst = ["models", "figures", "./figures/nCells", "./figures/llh",
    "./figures/xzUncertainty", "./figures/xzMasked", "./figures/xzMean",
    "./figures/xyUncertainty", "./figures/xyMasked", "./figures/xyMean",
    "./figures/xyContour", "./figures/xzContour",
    "./figures/phi", "./figures/nCells"
        ]

    for i in TD_parameters["z0"]
        push!(FigDirLst,"./figures/xyVoronoi_"*string(i))
    end

    for i in TD_parameters["y0"]
        push!(FigDirLst,"./figures/xzVoronoi_"*string(i))
    end    

    for iDir in FigDirLst
        if isdir(iDir) == false
            mkdir(iDir)
        end
    end
end

function lonlat2xy( 
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    lon1, 
    lat1)
    # function lonlat2xy(lon0::Float64, lat0::Float64, beta::Float64, lon1::Array{Float64,2}, lat1::Array{Float64,2})
    # Convert lon, lat to x, y
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    # yurong 05/12/23 high latitude correction
    # xx = (lon1 .- lon0) .* re ./ r2d
    xx = cos.(lat1 ./ r2d) .* (lon1 .- lon0) .* re ./ r2d
    yy = (lat1 .- lat0) .* re ./ r2d
    x1 = (xx .- yy .* tan.(beta)) .* cos.(beta)
    y1 = x1 .* tan.(beta) .+ yy ./ cos.(beta)

    ###################CANNOT WORK!!!! HOW TO DEAL WITH IT?!!##############
    # if abs.(x1) < 0.01
    #     x1 = 0
    # end
    # if abs.(y1) < 0.01
    #     y1 = 0
    # end
    return x1, y1
end

function xy2lonlat(
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    x2, 
    y2
    )
    # Convert x, y to lon, lat 
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    yy = (y2 .- x2 .* tan.(beta)) .* cos.(beta)
    xx = yy .* tan.(beta) .+ x2 ./ cos.(beta)
    # yurong 05/12/23 high latitude correction
    # lon2 = xx .* r2d ./ re .+ lon0
    lat2 = yy .* r2d ./ re .+ lat0
    lon2 = xx .* r2d ./ re ./ cos.(lat2 ./ r2d) .+ lon0
    
    return lon2, lat2
end

# function interp1(
#     x::Array{Float64,1}, 
#     y::Array{Float64,1}, 
#     xx::Array{Float64,1}
#     )

#     xx_length = length(xx)
#     x_length = length(x)
#     # yy::Array{Float64,1} = zeros(xx_length)
#     yy = NaN .* xx
#     for i in 1:xx_length
#         for j in 1:x_length
#             if xx[i] >= x[j] && xx[i] < x[j + 1]
#                 # y1 = y[j] + (xx[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j])
#                 # println(xx[i], ' ', y1)
#                 yy[i] = y[j] + (xx[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j])
#             end
#         end
#     end
#     return yy
# end

function build_starting(
    TD_parameters::Dict{String,Any}, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}
    )

    # log uniform distribution, in Byrnes and Bezada, 2020, eq. 11

    xVec = dataStruct["xVec"]
    yVec = dataStruct["yVec"]
    zVec = dataStruct["zVec"]
    nCells = floor.(exp.(rand(1) * log(TD_parameters["max_cells"] / 
            TD_parameters["min_cells"]) .+ log(TD_parameters["min_cells"])))
    # nCells = floor.(rand(1)*(TD_parameters.max_cells- TD_parameters.min_cells)) .+ TD_parameters.min_cells

    xCell = min(xVec...) .+ (max(xVec...) - min(xVec...)) .* rand(Int(nCells[1]))
    yCell = min(yVec...) .+ (max(yVec...) - min(yVec...)) .* rand(Int(nCells[1]))
    zCell = min(zVec...) .+ (max(zVec...) - min(zVec...)) .* rand(Int(nCells[1]))


    if TD_parameters["prior"] == 1 
        # Uniform 
            # for absolute t*
        zeta = rand(Int(nCells[1])) .* TD_parameters["zeta_scale"]
            # for relative t*
        # zeta = rand(Int(nCells[1])) .* TD_parameters["zeta_scale"] * 2 .- TD_parameters["zeta_scale"]  
    elseif TD_parameters["prior"] == 2
        # Normal 
        zeta = rand(Normal(0, TD_parameters["zeta_scale"]), Int(nCells[1]))
    elseif TD_parameters["prior"] == 3
        # Exponential 
        zeta = -log.(rand(Int(nCells[1]))) .* TD_parameters["zeta_scale"]
    end

    model = Model(nCells[1],
        xCell, yCell, zCell,
        zeta,
        -1, -1, -1, -1, -1, -1
    )
    # setindex!(dataStruct, model["allSig"] .* ones(length(dataStruct["tS"])), "allSig")
    
    (model, dataStruct, valid) = evaluate(model, dataStruct, RayTraces, TD_parameters)
    valid = 1
    return model, dataStruct, valid
end

function evaluate(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    TD_parameters::Dict{String,Any}
    )
    valid       = 1
    likelihood  = 1
    phi         = 1
    model.phi = phi
    model.likelihood = likelihood

    if TD_parameters["debug_prior"] == 1
        return model, dataStruct, valid
    end

    (m, n) = size(RayTraces["rayX"])
    ptS = zeros(n)
    observed_traveltime = zeros(n,1)
    predicted_traveltime = zeros(n,1)
    Threads.@threads for i in 1:n
        zeta0 = Interpolation(TD_parameters, model, RayTraces["rayX"][:,i], RayTraces["rayY"][:,i], RayTraces["rayZ"][:,i]) 
        rayzeta = []
 
        endpoint = length(zeta0)
        rayzeta = 0.5 .* (zeta0[1:endpoint - 1] + zeta0[2:endpoint])
    
        rayl = RayTraces["rayL"][:,i]
        index = findall(x -> isnan(x), rayl)
        
        if isempty(index)
            ptS[i] = sum(rayl .* RayTraces["rayU"][:,i] .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* RayTraces["rayU"][:,i])
        else
            rayl = RayTraces["rayL"][1:index[1]-1,i]  
            rayu = RayTraces["rayU"][1:index[1]-1,i]
        # compare travel time
            ptS[i] = sum(rayl .* rayu .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* rayu)
        end
        observed_traveltime[i] = 1000*dataStruct["tS"][i]/dataStruct["allaveatten"][i]
    end

    (m, n) = size(RayTraces["rayZ"])

    tS = dataStruct["tS"]

    C = 0
    for k in 1:length(dataStruct["allSig"])
        C += (ptS - tS)[k].^2 .* 1.0 / dataStruct["allSig"][k][1]^2
    end
    model.phi = C
    # model.ptS = ptS
    # model.tS = tS
    # model.predicted_traveltime = predicted_traveltime
    # model.observed_traveltime = observed_traveltime

    # yurong 03/02/23 change llh fuction based on central limit theorem
    likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi))) - 
    sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)
    # likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi)) * length(tS)) - 
    # sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)


    model.likelihood = likelihood
    
    return model, dataStruct, valid
end

function IDW(
    TD_parameters::Dict{String,Any}, 
    xCell::Array{Float64,2}, 
    yCell::Array{Float64,2}, 
    zCell::Array{Float64,2}, 
    zetaCell::Array{Float64,2}, 
    X, #::Array{Float64,1},#StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, 
    Y, #::Array{Float64,1}, 
    Z, #::Array{Float64,1}, 
    nCells::Int64, 
    npoints::Int64)

    zeta = zeros(npoints, 1) 
    
    if nCells == 0
        return
    end

    X = X[1:npoints]

    xCell = xCell[1:nCells]
    yCell = yCell[1:nCells]
    zCell = zCell[1:nCells]
    zetaCell = zetaCell[1:nCells]

    if TD_parameters["add_yVec"] == 0
        for k in 1:npoints
            distance = (X[k] .- xCell).^2 + (Z[k] .- zCell).^2
            v_sum = sum(zetaCell ./ distance)
            inv_sum = sum(1 ./ distance)
            zeta[k] = v_sum / inv_sum
        end
    else
        for k in 1:npoints
            distance = sqrt.((X[k] .- xCell).^2 + (Y[k] .- yCell).^2 + (Z[k] .- zCell).^2)
            v_sum = sum(zetaCell ./ distance)
            inv_sum = sum(1 ./ distance)
            zeta[k] = v_sum / inv_sum
        end
    end
    return zeta
end

# function v_idw(x, y, z, mx, my, mz, mv)

#     v_sum   = zero(Float64)
#     inv_sum = zero(Float64)

#     @inbounds for i in 1:length(mx)

#         distance = (mx[i] - x)^2 + (my[i] - y)^2 + (mz[i] - z)^2
#         v_sum   += mv[i] / distance
#         inv_sum += 1 / distance

#     end

#     return v_sum / inv_sum

# end

function v_nearest(x, y, z, mx, my, mz, mv)

    v     = zero(Float64)
    mdist = 1e9

    @inbounds for i in 1:length(mx)

        distance = (mx[i] - x)^2 + (my[i] - y)^2 + (mz[i] - z)^2
        if distance < mdist
            mdist = distance
            v     = mv[i]
        end

    end
    return v

end

# function NearestInterpolation(
#     TD_parameters::parameters, 
#     xCell::Array{Float64,2}, 
#     yCell::Array{Float64,2}, 
#     zCell::Array{Float64,2}, 
#     zetaCell::Array{Float64,2}, 
#     X, #::Array{Float64,1}, 
#     Y, #::Array{Float64,1}, 
#     Z, #::Array{Float64,1}, 
#     nCells::Int64, 
#     npoints::Int64)
    
#     zeta = zeros(npoints, 1) 
    
#     if nCells == 0
#         return
#     end

#     X = X[1:npoints]

#     xCell = xCell[1:nCells]
#     yCell = yCell[1:nCells]
#     zCell = zCell[1:nCells]
#     zetaCell = zetaCell[1:nCells]


#     # find the nearest cell to the points in each ray
#     if TD_parameters.add_yVec == 0
#         for k in 1:npoints
#             tmp = (X[k] .- xCell).^2 + (Z[k] .- zCell).^2
#             zeta[k] = zetaCell[findmin(tmp)[2]]
#         end
#     else
#         for k in 1:npoints
#             tmp = (X[k] .- xCell).^2 + (Y[k] .- yCell).^2 + (Z[k] .- zCell).^2
#             zeta[k] = zetaCell[findmin(tmp)[2]]
#         end
#     end
#     return zeta
# end

function Interpolation(
    TD_parameters::Dict{String,Any},
    model::Model, 
    X, Y, Z
    )

    if .~isempty(findall(x -> isnan.(x), X))
        npoints = findall(x -> isnan.(x), X)[1] - 1
    else
        npoints = length(X)
    end
    if length(Y) == 1 # xzMap
        Y = Y .* ones(npoints)
    end
    if length(Z) == 1 # xyMap
        Z = Z .* ones(npoints)
    end
    # if length(X) == 1 # yzMap
    #     X = X .* ones(npoints)
    # end
    if TD_parameters["interp_style"] == 1  # Nearest Interpolation
        zeta = [ v_nearest(X[k], Y[k], Z[k], model.xCell, model.yCell, model.zCell, model.zeta) for k = 1:npoints ]
    elseif TD_parameters["interp_style"] == 2  # Inverse Distance Weighting(IDW) 
        # X[isnan.(X)] .= 99999
        # Y[isnan.(Y)] .= 99999
        # Z[isnan.(Z)] .= 99999
        zeta = IDW(TD_parameters, xCell, yCell, zCell, zeta, X, Y, Z, Int(model.nCells), npoints)
    end

    return zeta
end

function normalize(arr::Array, min_val::Number, max_val::Number)
    min_arr, max_arr = extrema(arr)
    return min_val .+ (max_val - min_val) .* ((max_arr .- arr) ./ (max_arr - min_arr))
end

# function CrossSectionInterpolation(
#     models::Array{Any,1}, # models in an individual chain
#     dataStruct::DataStruct, 
#     TD_parameters::parameters,
#     chain::Int64
#     )

#     makeplot = TD_parameters.plot_voronoi
#     ENV["GKSwstype"] = "nul"

#     for k in 1:length(models)
#         if TD_parameters.xzMap == 1
#             zeta_xz = []
#             for i in 1:length(dataStruct.zVec)
#                 tmp = Interpolation(TD_parameters, models[k], dataStruct.xVec, TD_parameters.y0, ones(length(dataStruct.xVec)) .* dataStruct.zVec[i])
#                 zeta_xz = vcat(zeta_xz, tmp)
#             end
#             models[k].zeta_xz = vec(zeta_xz)
#             if makeplot == 1
#                 p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(zeta_xz), c=:jet,linewidth=0.001, clims=(0,20), yflip=true)
#                 scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
#                 savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(k) * ("_xz_") * string(TD_parameters.y0) * (".png"))
#             end
#         end
#         if TD_parameters.xyMap == 1
#             zeta_xy = []
#             for i in 1:length(dataStruct.yVec)
#                 tmp = Interpolation(TD_parameters, models[k], dataStruct.xVec, ones(length(dataStruct.xVec)) .* dataStruct.yVec[i], TD_parameters.z0)
#                 zeta_xy = vcat(zeta_xy, tmp)
#             end
#             models[k].zeta_xy = vec(zeta_xy)
#             if makeplot == 1
#                 p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), vec(zeta_xy), linewidth=0.001, clims=(0,20), c=:jet)
#                 savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(k) * ("_xy_") * string(TD_parameters.z0) * (".png"))
#             end
#         end
#     end
#     return models
# end

function Plot_model(
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    model_mean::Array{Float64,2},  #used to be Array{Float64} since model_mean is a 1 dimensional vector, but it is now a 2D matrix
    TD_parameters::Dict{String,Any},
    cmax::Float64,
    l0::Int64,
    ModelType::String,
    CrossSection::String
    )

    ENV["GKSwstype"] = "nul"
    gr()
    closeenough = 2.0
    if ModelType == "Masked" 
        cmap = :jet
        cbtitle = "1000/Qp"
    elseif ModelType == "Uncertainty"
        cmap = :bone
        cbtitle = "σ"
    elseif ModelType == "Mean"
        cmap = :jet
        cbtitle = "1000/Qp"
    end

    dirname = "./figures/"*CrossSection*ModelType

    if TD_parameters["add_yVec"] == 0 # Tonga 2D
        println("****2D Plotting****")
        tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
        tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
        tmpy = Array{Float64,2}(tmpy)
        p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=cmap)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, dirname*"/model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, dirname*"/model_mean_2d_Tonga__0-300km_ray") 
    else # 3d
        if CrossSection == "xz"
            (m, n) = size(RayTraces["rayY"])
            nearrays = zeros(n)
            for i in 1:n
                yvec = RayTraces["rayY"][:,i] .- l0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
           
            if ModelType == "Uncertainty"
                title!("Model uncertainty on cross-section")
            else
                title!(ModelType*" model on cross-section")
            end
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap" * string(l0) * "km")

            title!(ModelType*" model on cross-section")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km" * string(TD_parameters["y0"][1]) * "km") 

            ylims!(p, (0, 660))
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_events" * string(l0) * "km") 

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            
            title!("Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")       
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")   
  
            # title!("model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km") 
        end

        if CrossSection == "xy"
            (m, n) = size(RayTraces["rayZ"])
            nearrays = zeros(n)
            for i in 1:n
                zvec = RayTraces["rayZ"][:,i] .- l0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
        
            p = contourf(dataStruct["xVec"], dataStruct["yVec"], vec(model_mean), xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
            title!("Model_"*ModelType*"_xyMap" * string(l0) * "km")  
            if ModelType == "Uncertainty"
                title!("Model uncertainty in map view (" * string(l0) * " km)")
            else
                title!(ModelType*" model in map view (" * string(l0) * " km)")
            end      
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap" * string(l0) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")
        
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")       
        end

    end


end

function Plot_model_with_uncertainty(
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    model_mean::Array{Float64,2},  #used to be Array{Float64} since model_mean is a 1 dimensional vector, but it is now a 2D matrix
    model_poststd::Array{Float64,2},
    TD_parameters::Dict{String,Any},
    cmax::Float64,
    l0::Int64,
    CrossSection::String
    )

    ENV["GKSwstype"] = "nul"
    gr()
    closeenough = 20.0
    ModelType = "Contour" 
    cmap = :jet
    cbtitle = "1000/Qp"
    threshold = 5.0
    dirname = "./figures/"*CrossSection*ModelType

    if TD_parameters["add_yVec"] == 0 # Tonga 2D
        println("****2D Plotting****")
        tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
        tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
        tmpy = Array{Float64,2}(tmpy)

        p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=cmap)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, dirname*"/model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, dirname*"/model_mean_2d_Tonga__0-300km_ray") 
    else # 3d
        if CrossSection == "xz"
            (m, n) = size(RayTraces["rayY"])
            nearrays = zeros(n)
            for i in 1:n
                yvec = RayTraces["rayY"][:,i] .- l0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
            
            alpha_mask = [if zi >= threshold zi else 1 end for zi in model_poststd]
            alpha_mask = normalize(alpha_mask, 0, 1)
            # alpha_mask = zeros(size(alpha_mask))

            p = contourf(
                vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), 
                linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", 
                yflip=true, clims=(0,cmax), c=cmap, alpha = alpha_mask, colorbar_title = cbtitle)
            
            contour!(
                p, vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_poststd), 
                levels = 20, linewidth=1, cbar = false, contour_labels = true)

            title!(ModelType*" model on cross-section")

            savefig(p, dirname*"/Model_"*ModelType*"_xzMap" * string(l0) * "km")

            title!(ModelType*" model on cross-section")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km" * string(TD_parameters["y0"][1]) * "km") 

            ylims!(p, (0, 660))
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_events" * string(l0) * "km") 

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            
            title!("Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")       
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")   
  
            # title!("model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km") 
        end

        if CrossSection == "xy"
            (m, n) = size(RayTraces["rayZ"])
            nearrays = zeros(n)
            for i in 1:n
                zvec = RayTraces["rayZ"][:,i] .- l0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
        
            alpha_mask = [if zi >= threshold zi else 1 end for zi in model_poststd]
            alpha_mask = normalize(alpha_mask, 0, 1)

            p = contourf(
                dataStruct["xVec"], dataStruct["yVec"], vec(model_mean), 
                xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, 
                clims=(0,cmax), c=cmap, alpha = alpha_mask, colorbar_title = cbtitle)
            contour!(
                p, vec(dataStruct["xVec"]), vec(dataStruct["yVec"]), vec(model_poststd), 
                levels = 10, linewidth=1, cbar = false, contour_labels = true)

            title!("Model_"*ModelType*"_xyMap" * string(l0) * "km")  
            if ModelType == "Uncertainty"
                title!("Model uncertainty in map view (" * string(l0) * " km)")
            else
                title!(ModelType*" model in map view (" * string(l0) * " km)")
            end      
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap" * string(l0) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")
        
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")       
        end

    end


end

# function PlotModelsOverIterations(
#     models, # models in an individual chain
#     dataStruct::DataStruct, 
#     TD_parameters::parameters,
#     chain::Int64,
#     cmax::Float64
#     )
#     ENV["GKSwstype"] = "nul"
#     gr()
#     max_std = 5 

    
#     n = length(models)

#     if TD_parameters.xzMap == true

#         voronoi_xz = vec(models[1].zeta_xz)
#         p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), voronoi_xz, linewidth=0.001, c=:jet, clims=(0,cmax), yflip=true)
#         ylims!(p, (0, 660))
#         scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
#         savefig(p, "./figures/xzVoronoi/chain#" * string(chain) * ("_1") * (".png"))

#         average_xz = vec(models[1].zeta_xz)
#         std_xz = zeros(length(average_xz))
#         mask = ones(length(std_xz))
#         for i = 1:length(std_xz)
#             if std_xz[i] > max_std
#                 mask[i] = NaN
#             end
#         end
#         mask_xz = mask .* average_xz
#         Plot_Contours(dataStruct, average_xz, TD_parameters, chain, 1, 1, 1, cmax)
#         Plot_Contours(dataStruct, std_xz, TD_parameters, chain, 1, 1, 2, cmax)
#         Plot_Contours(dataStruct, mask_xz, TD_parameters, chain, 1, 1, 3, cmax)
#     end

#     if TD_parameters.xyMap == true
#         voronoi_xy = vec(models[1].zeta_xy)
#         p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), voronoi_xy, clims=(0,cmax), c=:jet)
#         savefig(p, "./figures/xyVoronoi/chain#" * string(chain) * ("_1") * (".png"))

#         average_xy = vec(models[1].zeta_xy)
#         std_xy = zeros(length(average_xy))
#         mask = ones(length(std_xy))
#         for i = 1:length(std_xy)
#             if std_xy[i] > max_std
#                 mask[i] = NaN
#             end
#         end
#         # mask_xy = mask .* mean(average_xy)
#         mask_xy = mask .* average_xy
#         Plot_Contours(dataStruct, average_xy, TD_parameters, chain, 1, 2, 1, cmax)
#         Plot_Contours(dataStruct, std_xy, TD_parameters, chain, 1, 2, 2, cmax)
#         Plot_Contours(dataStruct, mask_xy, TD_parameters, chain, 1, 2, 3, cmax)
#     end

#     for i in 2:n
#         if TD_parameters.xzMap == true
#             voronoi_xz = vec(models[i].zeta_xz)
#             p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), voronoi_xz,linewidth=0.001, c=:jet,clims=(0,cmax),  yflip=true)
#             ylims!(p, (0, 660))
#             scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)    
#             savefig(p, "./figures/xzVoronoi/chain#" * string(chain) * ('_') * string(i) * (".png"))

#             all_xz = pmap(x -> models[x].zeta_xz, 1:i)
#             average_xz = mean(all_xz)
#             std_xz = std(all_xz)
#             mask = ones(length(std_xz))
#             for j = 1:length(std_xz)
#                 if std_xz[j] > max_std
#                     mask[j] = NaN
#                 end
#             end
#             mask_xz = mask .* average_xz
#             Plot_Contours(dataStruct, average_xz, TD_parameters, chain, i, 1, 1, cmax)
#             Plot_Contours(dataStruct, std_xz, TD_parameters, chain, i, 1, 2, cmax)
#             Plot_Contours(dataStruct, mask_xz, TD_parameters, chain, i, 1, 3, cmax)    
#         end
#         if TD_parameters.xyMap == true
#             voronoi_xy = vec(models[i].zeta_xy)
#             p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), voronoi_xy, clims=(0,cmax), c=:jet)
#             savefig(p, "./figures/xyVoronoi/chain#" * string(chain) * ('_') * string(i) * (".png"))

#             all_xy = pmap(x -> models[x].zeta_xy, 1:i)
#             average_xy = mean(all_xy)
#             std_xy = std(all_xy)
#             mask = ones(length(std_xy))
#             for j = 1:length(std_xy)
#                 if std_xy[j] > max_std
#                     mask[j] = NaN
#                 end
#             end
#             mask_xy = mask .* average_xy
#             Plot_Contours(dataStruct, average_xy, TD_parameters, chain, i, 2, 1, cmax)
#             Plot_Contours(dataStruct, std_xy, TD_parameters, chain, i, 2, 2, cmax)
#             Plot_Contours(dataStruct, mask_xy, TD_parameters, chain, i, 2, 3, cmax)    
#         end
#     end



# end

# function Plot_Contours(
#     dataStruct::DataStruct, 
#     plot_model, # ::Array{Float64,2}, 
#     TD_parameters::parameters,
#     chain::Int64,
#     k::Int64,
#     CrossSection::Int64,
#     type::Int64, # 1: model, 2: model uncertainty, 3: masked model
#     cmax::Float64
#     )

#     ENV["GKSwstype"] = "nul"
#     gr()

#     closeenough = 2.0

#     if TD_parameters.add_yVec == 0 # Tonga 2D
#         println("2dTonga")
#         tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
#         tmpy = vcat(zeros(size(dataStruct.dataX))', vec(dataStruct.edep)') 
#         tmpy = Array{Float64,2}(tmpy)
#         p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=:jet)

#         title!("model_mean_2d_Tonga")
#         savefig(p, "model_mean_2d_Tonga")

#         title!("model_mean_2d_Tonga_0-300km")
#         ylims!(p, (0, 300))
#         savefig(p, "model_mean_2d_Tonga_0-300km") 

#         ylims!(p, (0, 660))
#         plot!(p, tmpx, tmpy, color=:white, legend=false)

#         scatter!(p, dataStruct.dataX, zeros(size(dataStruct.dataX)), marker=:utri, color=:pink, label="station", markersize=6)
#         scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4)
#         title!("model_mean_2d_Tonga_ray")       
#         savefig(p, "model_mean_2d_Tonga_ray")   

#         title!("model_mean_2d_Tonga__0-300km_ray")
#         ylims!(p, (0, 300))
#         savefig(p, "model_mean_2d_Tonga__0-300km_ray") 

#     else # 3d
#         if CrossSection == 1 # xz

#             if type == 1
#                 appendname = "./figures/xzMap/model_mean_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_mean_chain#" * string(chain) * "_" * string(k)
#             elseif type == 2
#                 appendname = "./figures/xzUncertainty/model_uncertainty_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_uncertainty_chain#" * string(chain) * "_" * string(k)
#             elseif type == 3
#                 appendname = "./figures/xzMask/model_mask_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_mask_chain#" * string(chain) * "_" * string(k)
#             end

#             (m, n) = size(dataStruct.rayY)
#             y0 = TD_parameters.y0
#             nearrays = zeros(n)
#             for i in 1:n
#                 yvec = dataStruct.rayY[:,i] .- y0    
#                 if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
#                     nearrays[i] = 1
#                 end  
#                 if min((abs.(yvec) .- closeenough)...) < 1e-7
#                     nearrays[i] = 1
#                 end
#             end
#             if type == 2 # uncertainty map
#                 p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=:bone)
#             else # model, masked model 
#                 p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=:jet)
#             end 
#             ylims!(p, (0, 660))
#             scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
#             title!(titlename * "_xzMap")      
#             savefig(p, appendname * "_xzMap")   
#         end

#         if CrossSection == 2 # xy
#             if type == 1
#                 appendname = "./figures/xyMap/model_mean_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_mean_chain#" * string(chain) * "_" * string(k)
#             elseif type == 2
#                 appendname = "./figures/xyUncertainty/model_uncertainty_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_uncertainty_chain#" * string(chain) * "_" * string(k)
#             elseif type == 3
#                 appendname = "./figures/xyMask/model_mask_chain#" * string(chain) * "_" * string(k)
#                 titlename = "model_mask_chain#" * string(chain) * "_" * string(k)
#             end
#             (m, n) = size(dataStruct.rayZ)
#             z0 = TD_parameters.z0
#             nearrays = zeros(n)
#             for i in 1:n
#                 zvec = dataStruct.rayZ[:,i] .- z0    
#                 if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
#                     nearrays[i] = 1
#                 end  
#                 if min((abs.(zvec) .- closeenough)...) < 1e-7
#                     nearrays[i] = 1
#                 end
#             end
#             nearrays = Array{Bool,1}(nearrays)
#             tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
#             tmpy = vcat(vec(dataStruct.dataY)', vec(dataStruct.elatsY)') 
#             # tmpy = Array{Float64,2}(tmpy)
#             nearraysX = tmpx[:,nearrays]
#             nearraysY = tmpy[:,nearrays]
            
#             if type == 2 #uncertainty map
#                 p = contourf(dataStruct.xVec, dataStruct.yVec, vec(plot_model), xlabel="X(km)", ylabel="Y(km)", c=:bone)
#             else # model, masked model
#                 p = contourf(dataStruct.xVec, dataStruct.yVec, vec(plot_model), xlabel="X(km)", ylabel="Y(km)",clims=(0,cmax), c=:jet)
#             end
            
#             title!(titlename * "_xyMap" * string(TD_parameters.z0) * "km")        
#             savefig(p, appendname * "_xyMap" * string(TD_parameters.z0) * "km")

#             plot!(p, tmpx, tmpy, color=:white, legend=false)
#             plot!(p, nearraysX, nearraysY, color=:forestgreen)
#             scatter!(p, dataStruct.dataX, dataStruct.dataY, shape=:utri, color=:pink, label="station", markersize=6)
#             scatter!(p, dataStruct.elonsX, dataStruct.elatsY, shape=:o, color=:lightblue, label="events", markersize=4)
#             title!(titlename * "_xyMap_ray" * string(TD_parameters.z0) * "km")
        
#             savefig(p, appendname * "_xyMap_ray" * string(TD_parameters.z0) * "km")       
#         end

#     end
# end

function plot_model_hist(
    model_hist, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}},
    TD_parameters::Dict{String,Any}, 
    cmax::Float64
    )

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m_xz = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))

            mcount = 0

            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ v_nearest(xs, l0, zs,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]
                    append!(m_xz,[m])

                end
            end

            model_mean_xz   = mean(m_xz)
            poststd_xz      = std(m_xz)
            mask_xz         = ones(size(poststd_xz))
            for i = 1:length(poststd_xz)
                if poststd_xz[i] > 5
                    mask_xz[i] = NaN
                end
            end
            mask_model_xz = mask_xz .* model_mean_xz

            # @time Plot_model(dataStruct, RayTraces, model_mean_xz, TD_parameters, cmax, l0, "Mean", "xz")
            # @time Plot_model(dataStruct, RayTraces, poststd_xz, TD_parameters, maximum(poststd_xz), l0, "Uncertainty", "xz")
            # @time Plot_model(dataStruct, RayTraces, mask_model_xz, TD_parameters, cmax, l0, "Masked", "xz")
            @time Plot_model_with_uncertainty(
                dataStruct, RayTraces, model_mean_xz, poststd_xz, 
                TD_parameters, cmax, l0,  "xz" 
                )
        end
    end

    if TD_parameters["xyMap"] == true
        for l0 in TD_parameters["z0"]
            m_xy = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            mcount = 0

            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ v_nearest(xs, ys, l0,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                    append!(m_xy,[m])

                end
            end

            model_mean_xy   = mean(m_xy)
            poststd_xy      = std(m_xy)
            mask_xy         = ones(size(poststd_xy))
            for i = 1:length(poststd_xy)
                if poststd_xy[i] > 5
                    mask_xy[i] = NaN
                end
            end
            mask_model_xy = mask_xy .* model_mean_xy

            # @time Plot_model(dataStruct, RayTraces, model_mean_xy, TD_parameters, cmax, l0, "Mean", "xy")
            # @time Plot_model(dataStruct, RayTraces, poststd_xy, TD_parameters, maximum(poststd_xy), l0, "Uncertainty", "xy")
            # @time Plot_model(dataStruct, RayTraces, mask_model_xy, TD_parameters, cmax, l0, "Masked", "xy")
            @time Plot_model_with_uncertainty(
                dataStruct, RayTraces, model_mean_xy, poststd_xy, 
                TD_parameters, cmax, l0, "xy" 
                )
        end
    end

end

function plot_model_hist_weighted(
    model_hist, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}},
    TD_parameters::Dict{String,Any}, 
    cmax::Float64
    )

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m_xz = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))

            mcount = 0

            llh = []
            for i = 1:length(model_hist)
                
                for j = 1:length(model_hist[i])

                    m  = [ v_nearest(xs, l0, zs,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]
                    append!(m_xz,[m])

                    push!(llh, model_hist[i][j].likelihood)
                
                end

            end

            weight_llh      = llh / sum(llh)
            # yurong 03/13/23 could still try wsum in package StatsBase to calculate the weighted average
            # and try weighted std
            # weight_llh      = [ones(size(m_xz[i]))*weight_llh[i] for i in 1:length(weight_llh)]
            # model_mean_xz   = wsum(m_xz, weight_llh)
            model_mean_xz   = zeros(size(m_xz[1]))
            for i in 1:length(m_xz)
                model_mean_xz += m_xz[i]*weight_llh[i]
            end
            poststd_xz      = std(m_xz)
            mask_xz         = ones(size(poststd_xz))
            for i = 1:length(poststd_xz)
                if poststd_xz[i] > 5
                    mask_xz[i] = NaN
                end
            end
            mask_model_xz = mask_xz .* model_mean_xz

            @time Plot_model(dataStruct, RayTraces, model_mean_xz, TD_parameters, cmax, l0, "Mean", "xz")
            @time Plot_model(dataStruct, RayTraces, poststd_xz, TD_parameters, maximum(poststd_xz), l0, "Uncertainty", "xz")
            @time Plot_model(dataStruct, RayTraces, mask_model_xz, TD_parameters, cmax, l0, "Masked", "xz")
        end
    end

    if TD_parameters["xyMap"] == true
        
        for l0 in TD_parameters["z0"]
            m_xy = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            mcount = 0

            llh = []
            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ v_nearest(xs, ys, l0,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                    append!(m_xy,[m])

                    push!(llh, model_hist[i][j].likelihood)

                end

            end

            weight_llh      = llh / sum(llh)
            model_mean_xy   = zeros(size(m_xy[1]))
            for i in 1:length(m_xy)
                model_mean_xy += m_xy[i]*weight_llh[i]
            end
            poststd_xy      = std(m_xy)
            mask_xy         = ones(size(poststd_xy))
            for i = 1:length(poststd_xy)
                if poststd_xy[i] > 5
                    mask_xy[i] = NaN
                end
            end
            mask_model_xy = mask_xy .* model_mean_xy

            @time Plot_model(dataStruct, RayTraces, model_mean_xy, TD_parameters, cmax, l0, "Mean", "xy")
            @time Plot_model(dataStruct, RayTraces, poststd_xy, TD_parameters, maximum(poststd_xy), l0, "Uncertainty", "xy")
            @time Plot_model(dataStruct, RayTraces, mask_model_xy, TD_parameters, cmax, l0, "Masked", "xy")
        end
    end

end

function plot_voronoi(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    TD_parameters::Dict{String,Any}, 
    chain::Int64, 
    iter::Float64
    )
    ENV["GKSwstype"] = "nul"

    cmap    = :jet
    cbtitle = "1000/Qp"
    cmax    = TD_parameters["cmax"]

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))
            m  = [ v_nearest(xs, l0, zs,
                model.xCell, model.yCell, model.zCell, model.zeta)
                for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(m), 
                 linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, 
                 clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xzVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter)))

            # add rays, EQs, and stations
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            plot!(p, tmpx, tmpy, color=:white, legend=false)
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xzVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter))*"ray")

        end
    end

    if TD_parameters["xyMap"] == true
        for l0 in TD_parameters["z0"]
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            m  = [ v_nearest(xs, ys, l0,
                model.xCell, model.yCell, model.zCell, model.zeta)
                for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["yVec"]), vec(m),
                 xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, clims=(0,cmax), 
                 c=cmap, colorbar_title = cbtitle)
            title!("Map View Z="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xyVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter)))

            # add rays, EQs, and stations
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])')
            plot!(p, tmpx, tmpy, color=:white, legend=false)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xyVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter))*"ray")
        end
    end

end

function plot_convergence(TD_parameters::Dict{String,Any})
# plot the convergence of nCells and phi value over iterations for each chain
# saved in ./figures/nCells and ./figures/phi
    ENV["GKSwstype"] = "nul"

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
            # plot nCells and phi over iterations
            cellnumber_list = load_model["nCells"]
            phi_list        = load_model["phi"]

            p1 = plot(1:length(cellnumber_list),cellnumber_list)
            title!(p1,"nCells Convergence in Chain" *string(chain))
            xlabel!(p1,"Iterations")
            ylabel!(p1,"Number of Cells")
            p2 = plot(1:length(phi_list),phi_list)
            title!(p2,"Phi Convergence in Chain" *string(chain))
            xlabel!(p2,"Iterations")
            ylabel!(p2,"Phi")
            savefig(p1,"./figures/nCells/nCells_chain" * string(chain))
            savefig(p2,"./figures/phi/phi_chain" * string(chain))
 
        end
    end
end

function plot_llh(models)
    ENV["GKSwstype"] = "nul"
    
    llh_all = []
    for i in 1:length(models)
        llh = []
        for j in 1:length(models[i])
            push!(llh, models[i][j].likelihood)
        end
        push!(llh_all, llh)
        p = histogram(llh,bins=10,xlabel="likelihood", ylabel="normalized ratio",normalize=true)
        title!(p,"llh_chain"*string(i))
        savefig(p,"./figures/llh/chain"*string(i)*".png")
    end
    p = histogram(llh_all,xlabel="likelihood", ylabel="normalized ratio",normalize=true)

    title!(p,"llh_all")
    savefig(p,"./figures/llh/llh_all.png")

end
