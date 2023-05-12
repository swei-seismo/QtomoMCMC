using JLD, HDF5
using Distributed 
using Glob
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")

@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")

function deg2rad(deg)
    return deg * π / 180
end

# Haversine formula to calculate great-circle distance between two points
function haversine_distance(lat1, lon1, lat2, lon2)
    R = 6371  # Earth's mean radius in km
    
    lat1_rad, lon1_rad = deg2rad(lat1), deg2rad(lon1)
    lat2_rad, lon2_rad = deg2rad(lat2), deg2rad(lon2)
    
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    a = sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
    c = 2 * atan(sqrt(a), sqrt(1-a))
    return R * c
end

# Calculate the distance between a point and a line segment
function point_line_distance(latp, lonp, lat1, lon1, lat2, lon2)
    d1 = haversine_distance(latp, lonp, lat1, lon1)
    d2 = haversine_distance(latp, lonp, lat2, lon2)
    d3 = haversine_distance(lat1, lon1, lat2, lon2)

    # Check if point is close to one of the end points
    if d1^2 >= d2^2 + d3^2 || d2^2 >= d1^2 + d3^2
        return min(d1, d2)
    end

    # Calculate the area of the triangle using Heron's formula
    s = (d1 + d2 + d3) / 2
    area = sqrt(s * (s - d1) * (s - d2) * (s - d3))

    # Calculate the distance between the point and the line segment
    return 2 * area / d3
end

lat0 = -23.1000
lon0 = 174.6000
beta = 0.463647609

println("--------Loading Data-------")
@time TD_parameters = define_TDstructrure()
@time (dataStruct, RayTraces) = load_data_Tonga(TD_parameters)
make_dir(TD_parameters)
println("--------Data Loaded-------")

# end2=( 178.10/-16.32 177.60/-17.32 177.10/-18.32 )
# end1=( 186.1/-20.31 185.6/-21.31 185.1/-22.31 )

line = [((178.10,-16.32),(186.1,-20.31)),
        ((177.60,-17.32),(185.6,-21.31)),
        ((177.10,-18.32),(185.1,-22.31))]

# Find points within 50 km of the line segment
threshold = 50  # Distance threshold in km

for i in 1:length(line)
    # Define the line segment with start and end points (lat1, lon1) and (lat2, lon2)
    lat1, lon1 = line[i][1][2], line[i][1][1]
    lat2, lon2 = line[i][2][2], line[i][2][1]

    nearby_points = []

    for j in 1:length(dataStruct["elats"])
        latp, lonp, depp = dataStruct["elats"][j], dataStruct["elons"][j], dataStruct["edep"][j]
        dist = point_line_distance(latp, lonp, lat1, lon1, lat2, lon2)
        if dist <= threshold
            push!(nearby_points, (lonp, latp, depp))
        end
    end
    # println(length(nearby_points),length(unique(nearby_points)))

    open("EQ_$(threshold)_line$(i).txt", "w") do io
        for i in unique(nearby_points)
            write(io,string(i[1])," ", string(i[2])," ", string(i[3]),"\n")
        end
    end

end