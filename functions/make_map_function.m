function make_map_function(data_coord, map_region, title_text, colorCell)
% =========================================================================
% make_map_function.m
% -------------------------------------------------------------------------
% Purpose:
%   Plot selected CRUTEM stations on a map.
%
% Inputs:
%   data_coord : matrix with columns
%       1. station ID
%       2. latitude
%       3. longitude
%       4. height
%
%   map_region : map domain used by worldmap, e.g.
%       'World', 'Africa', 'Europe', 'North America', ...
%
%   title_text : optional figure title
%
%   colorCell : optional cell array with colors:
%       {land_color, sea_color, point_color}
%
% Notes:
%   - Requires the MATLAB Mapping Toolbox.
%   - Exports should be handled by the calling script, not inside this
%     function.
% =========================================================================

if nargin < 2 || isempty(map_region)
    map_region = 'World';
end

if nargin < 3 || isempty(title_text)
    title_text = '';
end

default_land  = [0.15 0.5 0.15];
default_sea   = [0.729411780834198 0.831372559070587 0.95686274766922];
default_point = [1 0 0];

if nargin == 4 && ~isempty(colorCell) && iscell(colorCell)
    land_color  = colorCell{1};
    sea_color   = colorCell{2};
    point_color = colorCell{3};
else
    land_color  = default_land;
    sea_color   = default_sea;
    point_color = default_point;
end

n = size(data_coord,1);

figure('Color',sea_color);

h = worldmap(map_region);
setm(h,'FFaceColor',sea_color);

land = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(land,'FaceColor',land_color);

[Stations(1:n).Geometry] = deal('Point');

for i = 1:n
    Stations(i).Lat  = data_coord(i,2);
    Stations(i).Lon  = data_coord(i,3);
    Stations(i).Name = num2str(data_coord(i,1));
end

geoshow(Stations, ...
    'Marker','.', ...
    'Color',point_color, ...
    'MarkerFaceColor',point_color);

if ~isempty(title_text)
    title(title_text,'FontSize',12);
end

fig = gcf;
fig.InvertHardcopy = 'off';

end