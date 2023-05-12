using MAT,DelimitedFiles,Interpolations,Plots
include("MCsub.jl")

function load_data_Tonga(TD_parameters::Dict{String,Any})
    allLats, allLons    = [], [] #station latitudes and longitudes
    allTS, allSig, allaveatten  = [], [], []
    elats, elons, edep          = [], [], [] #event latitudes and longitudes
    
    traces    = load("./Data/traces.jld")
    elats       = Array{Float64}(traces["EventLatitude"])
    elons       = Array{Float64}(traces["EventLongitude"])
    edep        = Array{Float64}(traces["EventDepth"])
    allLats     = Array{Float64}(traces["latitude"])
    allLons     = Array{Float64}(traces["longitude"])
    allTS       = Array{Float64}(traces["tStar"])
    allaveatten = Array{Float64}(traces["aveatten"])
    allSig      = Array{Float64}(traces["error"])

    #drop to singleton dimensions
    allTS       = Vector{Float64}(vec(allTS))
    allLats     = Vector{Float64}(vec(allLats))
    allLons     = Vector{Float64}(vec(allLons))
    allaveatten = Vector{Float64}(vec(allaveatten))
    allSig      = Vector{Float64}(vec(allSig))

    lat0 = -23.1000
    lon0 = 174.6000
    # lat0 = -19.1400
    # lon0 = 175.7700
    beta = 0.463647609
    
    #stations and events coordicates in cartesian system
    (dataX, dataY)      = lonlat2xy(lon0, lat0, beta, allLons, allLats)
    (elonsX, elatsY)    = lonlat2xy(lon0, lat0, beta, elons, elats)

    #load coastline data
    #????not used yet
    # coastlines = load("./Data/coastlines.jld")
    # coastlon = coastlines["coastlon"]
    # coastlat = coastlines["coastlat"]
    # (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    #study area
    minX = minimum(dataX) - TD_parameters["buffer"]
    maxX = maximum(dataX) + TD_parameters["buffer"]
    minY = minimum(dataY) - TD_parameters["buffer"]
    maxY = maximum(dataY) + TD_parameters["buffer"]

    xVec = minX:TD_parameters["XYnodeSpacing"]:maxX
    yVec = minY:TD_parameters["XYnodeSpacing"]:maxY
    zVec = TD_parameters["min_depth"]:TD_parameters["ZnodeSpacing"]:TD_parameters["max_depth"]

    #load LAB discontinuity
    #??? not used yet
    LAB = load("./Data/LAB_discontinuity.jld")
    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

    #load raypaths
    raypath=load("./Data/raypaths.jld")
    x = raypath["x"]
    y = raypath["y"]
    z = raypath["z"]
    U = raypath["u"]

    #raylength and slowness for each segment
    rayl = sqrt.((x[1:end - 1,:] - x[2:end,:]).^2 +
    (y[1:end - 1,:] - y[2:end,:]).^2 +
    (z[1:end - 1,:] - z[2:end,:]).^2)
    rayu = 0.5 .* (U[1:end - 1,:] + U[2:end,:])

    dataStruct = Dict(
        "tS" => allTS   ::Array{Float64,1},
        "allaveatten" => allaveatten   ::Array{Float64,1},
        "allLats" => allLats        ::Array{Float64,1},    #Latitude for the stations in each event-station pair
        "allLons" => allLons        ::Array{Float64,1},    #Longitude for the stations in each event-station pair
        "allSig" => allSig      ::Array{Float64,1},      #ATTENTION!!! remain unknown!
        "dataX" => dataX        ::Array{Float64,1},        #Station position in each event-station pair
        "dataY" => dataY        ::Array{Float64,1},
        "xVec" => xVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "yVec" => yVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "zVec" => zVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "elonsX" => elonsX  ::Array{Float64,2},
        "elatsY" => elatsY  ::Array{Float64,2},
        "elons" => elons    ::Array{Float64,2},
        "elats" => elats    ::Array{Float64,2},
        "edep" => edep      ::Array{Float64,2},
        # "coastX" => coastX  ::Array{Float64,2},
        # "coastY" => coastY  ::Array{Float64,2},
    )

    RayTraces = Dict(
        "rayX" => x         ::Array{Float64,2},
        "rayY" => y         ::Array{Float64,2},
        "rayZ" => z         ::Array{Float64,2},
        "rayL" => rayl      ::Array{Float64,2},
        "rayU" => rayu      ::Array{Float64,2},
        "U" => U            ::Array{Float64,2}
    )

    return dataStruct, RayTraces
end

function load_synthetic_data_Tonga(TD_parameters::Dict{String,Any})
    allLats, allLons    = [], [] #station latitudes and longitudes
    allTS, allSig, allaveatten  = [], [], []
    elats, elons, edep          = [], [], [] #event latitudes and longitudes
    
    traces    = load("./Data/synthetic_traces.jld")
    elats       = Array{Float64}(traces["EventLatitude"])
    elons       = Array{Float64}(traces["EventLongitude"])
    edep        = Array{Float64}(traces["EventDepth"])
    allLats     = Array{Float64}(traces["latitude"])
    allLons     = Array{Float64}(traces["longitude"])
    allTS       = Array{Float64}(traces["tStar"])
    allaveatten = Array{Float64}(traces["aveatten"])
    allSig      = Array{Float64}(traces["error"])

    #drop to singleton dimensions
    allTS       = Vector{Float64}(vec(allTS))
    allLats     = Vector{Float64}(vec(allLats))
    allLons     = Vector{Float64}(vec(allLons))
    allaveatten = Vector{Float64}(vec(allaveatten))
    allSig      = Vector{Float64}(vec(allSig))

    lat0 = -23.1000
    lon0 = 174.6000
    beta = 0.463647609

    #stations and events coordicates in cartesian system
    (dataX, dataY)      = lonlat2xy(lon0, lat0, beta, allLons, allLats)
    (elonsX, elatsY)    = lonlat2xy(lon0, lat0, beta, elons, elats)

    #load coastline data
    #????not used yet
    # coastlines = load("./Data/coastlines.jld")
    # coastlon = coastlines["coastlon"]
    # coastlat = coastlines["coastlat"]
    # (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    #study area
    minX = minimum(dataX) - TD_parameters["buffer"]
    maxX = maximum(dataX) + TD_parameters["buffer"]
    minY = minimum(dataY) - TD_parameters["buffer"]
    maxY = maximum(dataY) + TD_parameters["buffer"]

    xVec = minX:TD_parameters["XYnodeSpacing"]:maxX
    yVec = minY:TD_parameters["XYnodeSpacing"]:maxY
    zVec = TD_parameters["min_depth"]:TD_parameters["ZnodeSpacing"]:TD_parameters["max_depth"]

    #load LAB discontinuity
    #??? not used yet
    LAB = load("./Data/LAB_discontinuity.jld")
    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

    #load raypaths
    raypath=load("./Data/synthetic_raypaths.jld")
    x = raypath["x"]
    y = raypath["y"]
    z = raypath["z"]
    U = raypath["u"]

    #raylength and slowness for each segment
    rayl = sqrt.((x[1:end - 1,:] - x[2:end,:]).^2 +
    (y[1:end - 1,:] - y[2:end,:]).^2 +
    (z[1:end - 1,:] - z[2:end,:]).^2)
    rayu = 0.5 .* (U[1:end - 1,:] + U[2:end,:])

    dataStruct = Dict(
        "tS" => allTS   ::Array{Float64,1},
        "allaveatten" => allaveatten   ::Array{Float64,1},
        "allLats" => allLats        ::Array{Float64,1},    #Latitude for the stations in each event-station pair
        "allLons" => allLons        ::Array{Float64,1},    #Longitude for the stations in each event-station pair
        "allSig" => allSig      ::Array{Float64,1},      #ATTENTION!!! remain unknown!
        "dataX" => dataX        ::Array{Float64,1},        #Station position in each event-station pair
        "dataY" => dataY        ::Array{Float64,1},
        "xVec" => xVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "yVec" => yVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "zVec" => zVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "elonsX" => elonsX  ::Array{Float64,2},
        "elatsY" => elatsY  ::Array{Float64,2},
        "elons" => elons    ::Array{Float64,2},
        "elats" => elats    ::Array{Float64,2},
        "edep" => edep      ::Array{Float64,2},
        # "coastX" => coastX  ::Array{Float64,2},
        # "coastY" => coastY  ::Array{Float64,2},
    )

    RayTraces = Dict(
        "rayX" => x         ::Array{Float64,2},
        "rayY" => y         ::Array{Float64,2},
        "rayZ" => z         ::Array{Float64,2},
        "rayL" => rayl      ::Array{Float64,2},
        "rayU" => rayu      ::Array{Float64,2},
        "U" => U            ::Array{Float64,2}
    )

    return dataStruct, RayTraces
end
