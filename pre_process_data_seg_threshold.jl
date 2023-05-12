# load information from raypaths.p and p_tstar.dat, including stations.lst
# output raypaths.jld and traces.jld for MCMC input
using DelimitedFiles,Interpolations,JLD

include("DefStruct.jl")
include("define_TDstructure.jl")
include("MCsub.jl")

include("load_3Dvel.jl")
include("define_TDstructure.jl")

@time TD_parameters = define_TDstructrure()
@time itp = load_3Dvel(TD_parameters)

function load_raypath()
    loadraypath = readlines("./Data/raypaths.p")

    #load raypaths from raypaths
    #interpolate 3D slowness for each point in a ray
    X, Y, Z, U = [], [], [], []
    ix, iy, iz, iu = [], [], [], []
    n=0
    for i in 1:length(loadraypath)
        sign = split(loadraypath[i])[1]
        if sign == "1234567"
            n+=1
            if ix == []
                continue
            end
            iu = itp.(ix,iy,iz)
            push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)
            ix, iy, iz, iu = [], [], [], []
        else
            iix = parse(Float64,split(loadraypath[i])[1])
            iiy = parse(Float64,split(loadraypath[i])[2])
            iiz = parse(Float64,split(loadraypath[i])[3])
            append!(ix,iix); append!(iy,iiy); append!(iz,iiz)
        end
        
        if i == length(loadraypath)
            iu = itp.(ix,iy,iz)
            push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)
        end
    end

    #interpolate the ray if individual segment is shorter than seg_len [km]
    seg_len = 10.0
    for i in 1:length(U)
        for j in 2:length(U[i])
            ilength = sqrt((X[i][j]-X[i][j-1])^2+(Y[i][j]-Y[i][j-1])^2+(Z[i][j]-Z[i][j-1])^2)
            if ilength > seg_len
                #note: use round() instead of ceil() since there are many 10.5 km segments when 10 km is set to the threshold
                seg_num     = Int(ceil(ilength/seg_len))
                insert_x    = range(X[i][j-1],stop=X[i][j],length=seg_num+1)
                insert_y    = range(Y[i][j-1],stop=Y[i][j],length=seg_num+1)
                insert_z    = range(Z[i][j-1],stop=Z[i][j],length=seg_num+1)
                insert_u    = itp.(insert_x,insert_y,insert_z)
                for k in 2:seg_num
                    insert!(X[i],j+k-2,insert_x[k])
                    insert!(Y[i],j+k-2,insert_y[k])
                    insert!(Z[i],j+k-2,insert_z[k])
                    insert!(U[i],j+k-2,insert_u[k])
                end
            end
        end
    end

    #rearrange the raypaths and according slowness to be a matrix
    #fill the empty cell with NaN
    x = fill(NaN,(maximum(length.(U)),length(U)))
    y = fill(NaN,(maximum(length.(U)),length(U)))
    z = fill(NaN,(maximum(length.(U)),length(U)))
    u = fill(NaN,(maximum(length.(U)),length(U)))

    for i in 1:length(U)
        x[1:length(U[i]),i]=X[i]
        y[1:length(U[i]),i]=Y[i]
        z[1:length(U[i]),i]=Z[i]
        u[1:length(U[i]),i]=U[i]
    end

    #save as raypaths.jld
    save("./Data/raypaths.jld","x",x,"y",y,"z",z,"u",u)
    
    #release memory
    X, Y, Z, U, loadraypath = nothing, nothing, nothing, nothing, nothing
    x, y, z, u = nothing, nothing, nothing, nothing
end

function load_traceinfo()
    loadatten = readlines("./Data/p_tstar.dat")
    loadsta = readlines("./Data/stations.lst")

    #load station latitude and longitude
    stalat, stalon = Dict(), Dict()
    for line in loadsta
        stalat[split(line)[1]] = parse(Float64,split(line)[2])
        stalon[split(line)[1]] = parse(Float64,split(line)[3])
    end

    #load all the information from p_tstar.dat (from t* inversion or synthetic dataset)
    station, EventLatitude, EventLongitude, EventDepth = [], [], [], []
    latitude, longitude, tStar, error, aveatten = [], [], [], [], []
    for line in loadatten
        ista            = split(line)[1]
        iEventLatitude  = parse(Float64,split(line)[2])
        iEventLongitude = parse(Float64,split(line)[3])
        iEventDepth     = parse(Float64,split(line)[4])
        itStar          = parse(Float64,split(line)[5])
        ierror          = parse(Float64,split(line)[6]) + 0.03 # systematic error from t* inversion
        istd            = parse(Float64,split(line)[7])
        iaveatten       = parse(Float64,split(line)[8])
        
        station         = vcat(station, ista)
        EventLatitude   = vcat(EventLatitude, iEventLatitude)
        EventLongitude  = vcat(EventLongitude, iEventLongitude)
        EventDepth      = vcat(EventDepth, iEventDepth)
        latitude        = vcat(latitude, stalat[ista])
        longitude       = vcat(longitude, stalon[ista])
        tStar           = vcat(tStar, itStar)
        error           = vcat(error, ierror)
        aveatten        = vcat(aveatten, iaveatten)
    end

    #save as traces.jld
    save("./Data/traces.jld","station",reshape(station,(length(station),1)),
        "EventLatitude",reshape(EventLatitude,(length(EventLatitude),1)),
        "EventLongitude",reshape(EventLongitude,(length(EventLongitude),1)),
        "EventDepth",reshape(EventDepth,(length(EventDepth),1)),
        "latitude",reshape(latitude,(length(latitude),1)),
        "longitude",reshape(longitude,(length(longitude),1)),
        "tStar",reshape(tStar,(length(tStar),1)),
        "error",reshape(error,(length(error),1)),
        "aveatten",reshape(aveatten,(length(aveatten),1)))
end

load_raypath()
load_traceinfo()