function [m,n]=compute_subplots(N)
if N==3
    m=3; n=1;
elseif round(sqrt(N))==sqrt(N)
    m=sqrt(N); n=sqrt(N);    
else
    m=round(sqrt(N))+1;
    n=m-1;
end