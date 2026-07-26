
%June 2019

function [d,ar,sigma2,res,y_hat]=estima_arfima(y,ar_part)
p1 = genpath('C:\WORKA\Climate\PAPERS\New_papers\Forecast\Forecast_fall_2023\IJF_2025\CODES_Revision_2024\ARFIMA\arfima_est_v2');
path(path,p1)
%p2= genpath('c:\WORKA\Climate\PAPERS\New_papers\Forecast\Calculos\ARFIMA\MFE');
p2= genpath('C:\WORKA\Climate\PAPERS\New_papers\Forecast\Forecast_fall_2023\IJF_2025\CODES_Revision_2024\ARFIMA\MFEToolbox'); %new 2018 version
path(path,p2)

if ar_part==0
[whittle] = arfima_estimate(y,'FWHI',[0 0]);
ar=[];
elseif ar_part==1
[whittle] = arfima_estimate(y,'FWHI',[1 0]);
ar=whittle.AR(1);
end
d=whittle.d(1);
sigma2=whittle.sigma2;
res=whittle.errors;
y_hat=y-whittle.errors;




 