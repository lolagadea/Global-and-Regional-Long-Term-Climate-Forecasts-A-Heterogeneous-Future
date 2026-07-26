function b=array_NaN(a)
n=isnan(a);
if sum(n)>=1
    b=0;
else
    b=1;
end