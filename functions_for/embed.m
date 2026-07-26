function X=embed(y,p)
%This function produces a matrix with yt and its p lags; similar to R
%function
t=length(y);
x=zeros(t-p,p);
%let's now fill the columns on x...
for j=1:p
    x(:,j)=y(p-j+1:t-j);
end
X=[y(p+1:t) x];