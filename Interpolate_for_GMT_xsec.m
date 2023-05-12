filename = 'Interpolated_Tonga_Atten_Model.txt'; 
data = readmatrix(filename, 'NumHeaderLines', 1);
longitude = data(:, 1);
latitude = data(:, 2);
depth = data(:, 3);
Model = data(:, 4); 
Poststd = data(:,5);
Mask = data(:,6);

Mask_threshold = 5;

F_Model = scatteredInterpolant(longitude, latitude, depth, Model, 'linear', 'none');
F_STD = scatteredInterpolant(longitude, latitude, depth, Poststd, 'linear', 'none');


for i = 1:3
    linename    = sprintf("line%s.dat",num2str(i));
    EQ_name     = sprintf("EQ_50_line%s.txt",num2str(i));
    linedata    = load(linename);
    EQdata      = load(EQ_name);
    new_lon     = linedata(:,1);
    new_lat     = linedata(:,2);
    EQ_lon      = EQdata(:,1);
    EQ_lat      = EQdata(:,2);
    EQ_dep      = EQdata(:,3);

    idx = find(new_lon<0);
    new_lon(idx) = new_lon(idx)+360;
    idx = find(EQ_lon<0);
    EQ_lon(idx) = EQ_lon(idx)+360;
    
    Model_name  = sprintf("./TongaAttenData/Model_xsec%s.dat",num2str(i));
    STD_name    = sprintf("./TongaAttenData/STD_xsec%s.dat",num2str(i));
    Mask_name   = sprintf("./TongaAttenData/Mask_xsec%s.dat",num2str(i));
    EQ_name     = sprintf("./TongaAttenData/EQ_xsec%s.dat",num2str(i));

    file_Model  = fopen(Model_name, 'w');
    file_STD    = fopen(STD_name, 'w');
    file_Mask   = fopen(Mask_name, 'w');
    file_EQ     = fopen(EQ_name, 'w');

    distance_km = haversine_distance(new_lat(1),new_lon(1),new_lat,new_lon);
    EQ_distance = haversine_distance(new_lat(1),new_lon(1),EQ_lat,EQ_lon);

    for j = 1:length(EQ_distance)
        fprintf(file_EQ, '%f %f\n', EQ_distance(j), EQ_dep(j));
    end

    for new_dep = 0:20:660
        new_Model = F_Model(new_lon,new_lat,ones(size(new_lon))*new_dep);
        new_STD = F_STD(new_lon,new_lat,ones(size(new_lon))*new_dep);
        new_Mask = size(new_Model);
        for iModel = 1:length(new_Model)
            if new_STD(iModel) > Mask_threshold
                new_Mask(iModel) = nan;
            else
                new_Mask(iModel) = new_Model(iModel);
            end
        end
        for j = 1:length(new_Model)
            fprintf(file_Model, '%f %f %f\n', distance_km(j), new_dep, new_Model(j));
            fprintf(file_STD, '%f %f %f\n', distance_km(j), new_dep, new_STD(j));
            fprintf(file_Mask, '%f %f %f\n', distance_km(j), new_dep, new_Mask(j));
        end

    end




    fclose(file_Model);
    fclose(file_STD);
    fclose(file_Mask);
    fclose(file_EQ);
        
end




