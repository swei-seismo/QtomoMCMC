using NearestNeighbors

function preprocess_interpolation(X, Y, Z, npoints::Int64)
    if length(Y) == 1 # xzMap
        Y = Y .* ones(npoints)
    end
    if length(Z) == 1 # xyMap
        Z = Z .* ones(npoints)
    end
    return X, Y, Z
end

# function cart_interpolation(
#     par::Dict{String,Any},
#     model::Model, 
#     X, Y, Z
#     )
#     npoints = calculate_npoints(X)
#     X, Y, Z = preprocess_interpolation(X, Y, Z, npoints)

#     zeta = par["interp_style"] == 1 ? 
#            [cart_v_nearest(X[k], Y[k], Z[k], model) for k = 1:npoints] :
#            IDW(par, model, X, Y, Z, npoints)

#     return zeta
# end

# function sph_interpolation(
#     par::Dict{String,Any},
#     model::Model, 
#     X, Y, Z
#     )
#     npoints = calculate_npoints(X)
#     X, Y, Z = preprocess_interpolation(X, Y, Z, npoints)

#     zeta = par["interp_style"] == 1 ? 
#            [sph_v_nearest(X[k], Y[k], Z[k], model) for k = 1:npoints] :
#            IDW(par, model, X, Y, Z, npoints)

#     return zeta
# end

function calculate_npoints(X)
    isempty(findall(isnan, X)) ? length(X) : findfirst(isnan, X) - 1
end

function interp_nearest(par,x, y, z, model::Model)
    if par["coordinates"] == 1
        return sph_v_nearest(x, y, z, model)
    else
        return cart_v_nearest(x, y, z, model)
    end
end

function cart_v_nearest(x, y, z, model::Model)
    v, mdist = zero(Float64), 1e9
    @inbounds for i in 1:length(model.xCell)
        distance = cart_dist(model.xCell[i], model.yCell[i], model.zCell[i], x, y, z)
        if distance < mdist
            mdist, v = distance, model.zeta[i]
        end
    end
    return v
end

function sph_v_nearest(x, y, z, model::Model)
    v, mdist = zero(Float64), 1e9
    @inbounds for i in 1:length(model.xCell)
        distance = sph_dist(model.xCell[i], model.yCell[i], model.zCell[i], x, y, z)
        if distance < mdist
            mdist, v = distance, model.zeta[i]
        end
    end
    return v
end

function IDW(
    par::Dict{String,Any}, 
    model::Model, 
    X, Y, Z, 
    npoints::Int64)
    zeta = zeros(npoints)
    if model.nCells == 0
        return zeta
    end

    for k in 1:npoints
        distances = par["add_yVec"] == 0 ? 
                    (X[k] .- model.xCell).^2 + (Z[k] .- model.zCell).^2 :
                    sqrt.((X[k] .- model.xCell).^2 + (Y[k] .- model.yCell).^2 + (Z[k] .- model.zCell).^2)

        v_sum = sum(model.zeta ./ distances)
        inv_sum = sum(1 ./ distances)
        zeta[k] = v_sum / inv_sum
    end
    return zeta
end


function geo_to_cart(lon, lat, dep)

    Re = 6371.0
    lon_rad = deg2rad.(lon)
    lat_rad = deg2rad.(lat)

    x = (Re .- dep) .* cos.(lat_rad) .* cos.(lon_rad)
    y = (Re .- dep) .* cos.(lat_rad) .* sin.(lon_rad)
    z = (Re .- dep) .* sin.(lat_rad)

    return x, y, z
end


function sph_interpolation(
        par::Dict{String,Any},
        model::Model, 
        x1, y1, z1)
    
    npoints = isempty(findall(isnan, x1)) ? length(x1) : findfirst(isnan, x1) - 1
    newx, newy, newz = geo_to_cart(model.xCell, model.yCell, model.zCell)
    data = [newx newy newz]
    newdata = transpose(data)
    tree = KDTree(newdata)
    zeta_lst = zeros(npoints)
    for k in 1:npoints
        newx1, newy1, newz1 = geo_to_cart(x1[k], y1[k], z1[k])
        query_point = [newx1, newy1, newz1]
        idx, dist = knn(tree, query_point, 1, true)
        zeta_lst[k] = model.zeta[idx[1]]
    end

    return zeta_lst
end

function cart_interpolation(
    par::Dict{String,Any},
    model::Model, 
    x1, y1, z1)
    npoints = isempty(findall(isnan, x1)) ? length(x1) : findfirst(isnan, x1) - 1
    newx, newy, newz = model.xCell, model.yCell, model.zCell
    data = [newx newy newz]
    newdata = transpose(data)
    tree = KDTree(newdata)
    zeta_lst = zeros(npoints)
    for k in 1:npoints
        newx1, newy1, newz1 = x1[k], y1[k], z1[k]
        query_point = [newx1, newy1, newz1]
        idx, dist = knn(tree, query_point, 1, true)
        zeta_lst[k] = model.zeta[idx[1]]
    end
    
    return zeta_lst
end