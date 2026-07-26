function long_dec=convert_long_gps_to_decimal(long_gps)
side=long_gps(8);
long_dec=str2double(long_gps(1:3))+str2double(long_gps(4:5))/60+str2double(long_gps(6:7))/3600;
if side=='E'
    long_dec=long_dec;
elseif side=='W'
    long_dec=-long_dec;
end