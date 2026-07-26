function z=my_diff(y,d)
t=length(y);
z=NaN(t-d,1);
for i=1:t-d
    z(i)=y(d+i)-y(i);
end
