function lat_dec=convert_lat_gps_to_decimal(lat_gps)

hemis=lat_gps(7);

lat_dec=str2double(lat_gps(1:2))+str2double(lat_gps(3:4))/60+str2double(lat_gps(5:6))/3600;

if hemis=='N'
    lat_dec=lat_dec;
elseif hemis=='S'
    lat_dec=-lat_dec;
end