struct Ray
    x; y; z 
end
    
# DataStructs are in dictionary form
struct DataStruct
    tS          ::Array{Float64,1}
    allaveatten ::Array{Float64,1}
    allLats     ::Array{Float64,1}      #Latitude for the stations in each event-station pair
    allLons     ::Array{Float64,1}      #Longitude for the stations in each event-station pair
    allSig      ::Array{Float64,1}      #error terms from t* inversion
    dataX       ::Array{Float64,1}  #Station position in each event-station pair
    dataY       ::Array{Float64,1}
    xVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    yVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    zVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    elonsX      ::Array{Float64,2}
    elatsY      ::Array{Float64,2}
    elons       ::Array{Float64,2}
    elats       ::Array{Float64,2}
    edep        ::Array{Float64,2}
    # coastX      ::Array{Float64,2}
    # coastY      ::Array{Float64,2}
    # ray         ::Array{Any,1}

end

# RayTraces are in dictionary form 
struct RAY
    rayX        ::Array{Float64,2}
    rayY        ::Array{Float64,2}
    rayZ        ::Array{Float64,2}
    rayL        ::Array{Float64,2}
    rayU        ::Array{Float64,2}
    U           ::Array{Float64,2}
end

mutable struct Model
    nCells      ::Float64
    xCell       ::Vector{Float64}
    yCell       ::Vector{Float64}
    zCell       ::Vector{Float64}
    zeta        ::Vector{Float64}
    phi         ::Float64
    # ptS         ::Vector{Float64}
    # tS          ::Vector{Float64}
    # predicted_traveltime
    # observed_traveltime
    likelihood  ::Float64
    action      ::Int64
    accept      ::Int64
    zeta_xz     ::Float64
    zeta_xy     ::Float64
end