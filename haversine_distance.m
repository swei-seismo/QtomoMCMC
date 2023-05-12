function distance_km = haversine_distance(lat1, lon1, lat2, lon2)
    % Convert latitude and longitude from degrees to radians
    lat1_rad = deg2rad(lat1);
    lon1_rad = deg2rad(lon1);
    lat2_rad = deg2rad(lat2);
    lon2_rad = deg2rad(lon2);

    % Earth's mean radius in km
    R = 6371;

    % Haversine formula
    dlat = lat2_rad - lat1_rad;
    dlon = lon2_rad - lon1_rad;

    a = sin(dlat/2).^2 + cos(lat1_rad) .* cos(lat2_rad) .* sin(dlon/2).^2;
    c = 2 .* atan2(sqrt(a), sqrt(1-a));

    % Calculate distance in km
    distance_km = R * c;
end

