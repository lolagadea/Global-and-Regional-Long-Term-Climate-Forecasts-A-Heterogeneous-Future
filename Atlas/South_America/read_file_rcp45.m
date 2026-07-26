clear
ncdisp('t_CORDEX-SAM_rcp45_mon_200601-210012.nc')
gcm_model  = ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','gcm_model');
gcm_variant  = ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','gcm_variant');
t  = ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','t');
time  = ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','time');
lat= ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','lat');
lon= ncread('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','lon');
m = ncreadatt('t_CORDEX-SAM_rcp45_mon_200601-210012.nc','t','missing_value');
%Adjust the time
dt= datetime(time, 'ConvertFrom', 'datenum');
dt2=datevec(dt);
years=dt2(:,1)+1850; months=dt2(:,2);
years=[years(2:end);2100]; months=[months(2:end);12];

[nlon,nlat,nt,nmodels]=size(t);
N=nlon*nlat*nt*nmodels;

for i=1:nlon
    for j=1:nlat
        for k=1:nt
            for v=1:nmodels
    if t(i,j,k,v)==m
        t(i,j,k,v)=NaN;
    end
            end
        end
    end
end

w=cos(lat*pi/180);
w=w/sum(w);

tw=NaN(nlon,nt,nmodels);
for k=1:nt
    for v=1:nmodels
        for i=1:nlon
            a=squeeze(t(i,:,k,v));
            aw=mult_matrix_nan(a,w);
            tw(i,k,v)=aw;
        end
    end
end


for k=1:nt
    b=reshape(tw(:,k,:),nlon*nmodels,1);
    temp_month(k,1)=nanmean(b);
end

years2=unique(years);
nyears=length(years2);
temp_years=NaN(nyears,1);
for i=1:nyears
    i
    temp_years(i,:)=mean(temp_month((i-1)*12+1:i*12,:));
end

TEMP.data.montly=temp_month;
TEMP.data.yearly=temp_years;
TEMP.years.montly=years;
TEMP.years.years=years2;
TEMP.months=months;
TEMP.models=gcm_model;
TEMP.variant=gcm_variant;
TEMP.lat=lat;
TEMP.lon=lon;

save('TEMP_2006_2100_South_America_rcp45.mat','TEMP','-v7.3')



