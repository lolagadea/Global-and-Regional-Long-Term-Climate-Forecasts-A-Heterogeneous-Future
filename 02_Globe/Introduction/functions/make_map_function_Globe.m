function make_map_function_Globe(data_coord)
% =========================================================================
% make_map_function_Globe.m
% -------------------------------------------------------------------------
% Purpose:
%   Plot selected CRUTEM stations on a world map.
%
% Input:
%   data_coord : matrix with columns
%       1. station ID
%       2. latitude
%       3. longitude
%       4. height
%
% Notes:
%   - This function preserves the visual style used in the paper:
%     blue sea, green land, and red station markers.
%   - It requires the MATLAB Mapping Toolbox.
% =========================================================================

n = size(data_coord,1);

figure('Color',[0.729411780834198 0.831372559070587 0.95686274766922]);

h = worldmap('World');
getm(h,'MapProjection');

geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15]);

[Stations(1:n).Geometry] = deal('Point');

for i = 1:n
    Stations(i).Lat = data_coord(i,2);
    Stations(i).Lon = data_coord(i,3);
    Stations(i).Name = num2str(data_coord(i,1));
end

geoshow(Stations, ...
    'Marker', '.', ...
    'Color', 'red', ...
    'MarkerFaceColor', 'red');

fig = gcf;
fig.InvertHardcopy = 'off';

end