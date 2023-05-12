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

if isdir("./TongaAttenData") == false
    mkdir("./TongaAttenData")
end

println("--------Loading Data-------")
@time TD_parameters = define_TDstructrure()
@time (dataStruct, RayTraces) = load_data_Tonga(TD_parameters)
make_dir(TD_parameters)
println("--------Data Loaded-------")

println("--------Loading Models-------")
model_hist = []
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
        push_models = load_model["model_hist"]
        push!(model_hist,push_models)
        println("--------Chain loaded-------")
    end
end
println("--------Models Loaded-------")

println("--------Interpolating Models-------")

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
            if poststd_xy[i] > std_threshold 
                mask_xy[i] = NaN
            end
        end
        mask_model_xy = mask_xy .* model_mean_xy
        open("./TongaAttenData/Tonga_Map_Mask_"*string(l0)*".txt", "w") do io
            for i in 1:length(dataStruct["xVec"])
                for j in 1:length(dataStruct["yVec"])
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *"  " * string(mask_model_xy[i,j]) *"\n")
                end
            end
        end

        open("./TongaAttenData/Tonga_Map_Model_"*string(l0)*".txt", "w") do io
            for i in 1:length(dataStruct["xVec"])
                for j in 1:length(dataStruct["yVec"])
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *"  " * string(model_mean_xy[i,j]) *"\n")
                end
            end
        end

        open("./TongaAttenData/Tonga_Map_Uncertainty_"*string(l0)*".txt", "w") do io
            for i in 1:length(dataStruct["xVec"])
                for j in 1:length(dataStruct["yVec"])
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    write(io, string(lon) *"  "* string(lat) *"  " * string(poststd_xy[i,j]) *"\n")
                end
            end
        end

        open("./TongaAttenData/Tonga_Map_Transparency_"*string(l0)*".txt", "w") do io
            for i in 1:length(dataStruct["xVec"])
                for j in 1:length(dataStruct["yVec"])
                    (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                    if poststd_xy[i,j]<std_threshold
                        transparency = 0
                    elseif poststd_xy[i,j]<13
                        transparency = (poststd_xy[i,j]-std_threshold)/(13-std_threshold)
                    else
                        transparency = 1
                    end
                    write(io, string(lon) *"  "* string(lat) *"  " * string(transparency) *"\n")
                end
            end
        end
    end
end

# ####false
# if TD_parameters["xzMap"] == true
#     for l0 in TD_parameters["y0"]
#         m_xz = []
#         m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))
#         mcount = 0

#         for i = 1:length(model_hist)

#             for j = 1:length(model_hist[i])

#                 m  = [ v_nearest(xs, l0, zs,
#                     model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
#                     for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]
#                 append!(m_xz,[m])

#             end
#         end

#         model_mean_xz   = mean(m_xz)
#         poststd_xz      = std(m_xz)
#         mask_xz         = ones(size(poststd_xz))
#         for i = 1:length(poststd_xz)
#             if poststd_xz[i] > 6
#                 mask_xz[i] = NaN
#             end
#         end
#         mask_model_xz = mask_xz .* model_mean_xz
#         open("Tongaxsec"*string(l0)*".txt", "w") do io
#             for i in 1:length(dataStruct["xVec"])
#                 for j in 1:length(dataStruct["zVec"])
#                     (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["zVec"][j])
#                     # for cross section BB'
#                     # newlon0 = -174.50
#                     # newlat0 = -21.50
#                     # (x,y) = lonlat2xy(newlon0,newlat0,beta,lon,lat)
#                     write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["zVec"][j]) *"  " * string(mask_model_xz[i,j]) *"\n")
#                 end
#             end
#         end
#     end
# end

