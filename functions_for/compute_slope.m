%this function computes the derivative of an polynomial model
%and build a t-test of the trend, 2012, revised January 2014

function[mslope]=compute_slope(beta,t)
trend=[1:1:t];
k=length(beta)-1;
if k==0
    mslope=0;
    test=NaN;
    ci_up=NaN;
    ci_low=NaN;
else
r=NaN(k,1);
r(1)=1;
for i=2:k
    r(i)=(i/t)*(sum(trend'.^(i-1)));
end
r=[0;r];
mslope=r'*beta;
end





