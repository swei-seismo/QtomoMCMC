using MAT,DelimitedFiles,Interpolations,Plots,JLD
include("MCsub.jl")

function load_3Dvel(TD_parameters::Dict{String,Any})
    vel = readlines("./Data/lau.vel")
    nnx, nny, nnz = parse(Int,split(vel[1])[1]), parse(Int,split(vel[1])[2]), parse(Int,split(vel[1])[3])
    lat0, lon0, beta = parse(Float64,split(vel[2])[1]), parse(Float64,split(vel[2])[2]), parse(Float64,split(vel[2])[3])
    lat, lon = Array{Float64,2}(undef,nnx,nny), Array{Float64,2}(undef,nnx,nny)
    for i = 1:nnx
        for j = 1:nny
            lat[i,j] = parse(Float64,split(vel[(i-1)*nny+j+2])[1])
            lon[i,j] = parse(Float64,split(vel[(i-1)*nny+j+2])[2])
        end
    end
    xx = [lon0 lon0]
    yy = [lat0 lat0]
    (dataX, dataY) = lonlat2xy(lon0, lat0, beta, lon, lat)
    (aa,bb) =  lonlat2xy(lon0, lat0, beta, xx, yy)

    z = [parse(Float64,i) for i in split(vel[nnx*nny+3])]
    sn, vps = Array{Float64,4}(undef,2,nnx,nny,nnz), Array{Float64,4}(undef,2,nnx,nny,nnz)
    for p = 1:2
        for i = 1:nnx
            for j = 1:nny
                for k = 1:nnz
                    vps[p,i,j,k] = parse(Float64,split(vel[(i-1+p*nnx)*nny+j+3])[k])
                    sn[p,i,j,k] = 1/vps[p,i,j,k]
                end
            end
        end
    end
    itp = interpolate((round.(dataX;digits=2)[:,1],round.(dataY;digits=2)[1,:],z),sn[1,:,:,:],Gridded(Linear()))
    return itp
end