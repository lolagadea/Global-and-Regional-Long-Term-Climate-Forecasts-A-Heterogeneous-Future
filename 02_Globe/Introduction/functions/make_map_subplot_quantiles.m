function make_map_subplot_quantiles(coord_q05,coord_q50,coord_q95)

% =========================================================================
% make_map_subplot_quantiles
% -------------------------------------------------------------------------
% Plot stations contributing to q05, q50, and q95 of the annual
% cross-sectional temperature distribution.
% =========================================================================

sea_color = [0.7294 0.8314 0.9568];

% ========================================================================
% q05
% ========================================================================

subplot(3,1,1)

n = size(coord_q05,1);

h = worldmap('World');

setm(h,'FFaceColor',sea_color);

geoshow('landareas.shp', ...
    'FaceColor',[0.15 0.5 0.15]);

[Stations(1:n).Geometry] = deal('Point');

for i = 1:n

    Stations(i).Lat  = coord_q05(i,2);
    Stations(i).Lon  = coord_q05(i,3);
    Stations(i).Name = num2str(coord_q05(i,1));

end

geoshow(Stations, ...
    'Marker','o', ...
    'MarkerSize',3, ...
    'MarkerFaceColor','red', ...
    'MarkerEdgeColor','red');

title('q05','FontSize',14)

clear Stations

% ========================================================================
% q50
% ========================================================================

subplot(3,1,2)

n = size(coord_q50,1);

h = worldmap('World');

setm(h,'FFaceColor',sea_color);

geoshow('landareas.shp', ...
    'FaceColor',[0.15 0.5 0.15]);

[Stations(1:n).Geometry] = deal('Point');

for i = 1:n

    Stations(i).Lat  = coord_q50(i,2);
    Stations(i).Lon  = coord_q50(i,3);
    Stations(i).Name = num2str(coord_q50(i,1));

end

geoshow(Stations, ...
    'Marker','o', ...
    'MarkerSize',3, ...
    'MarkerFaceColor','red', ...
    'MarkerEdgeColor','red');

title('q50','FontSize',14)

clear Stations

% ========================================================================
% q95
% ========================================================================

subplot(3,1,3)

n = size(coord_q95,1);

h = worldmap('World');

setm(h,'FFaceColor',sea_color);

geoshow('landareas.shp', ...
    'FaceColor',[0.15 0.5 0.15]);

[Stations(1:n).Geometry] = deal('Point');

for i = 1:n

    Stations(i).Lat  = coord_q95(i,2);
    Stations(i).Lon  = coord_q95(i,3);
    Stations(i).Name = num2str(coord_q95(i,1));

end

geoshow(Stations, ...
    'Marker','o', ...
    'MarkerSize',3, ...
    'MarkerFaceColor','red', ...
    'MarkerEdgeColor','red');

title('q95','FontSize',14)

end