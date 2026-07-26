function A=array_to_triang_sup_matrix(a,l)
k=length(a);
i=1;
j=2;
m=1;
A=NaN(l,l);
while(m<=k)
    A(i,j)=a(m);  
    if j==l
        i=i+1;
        j=i;
    end
    j=j+1;
    m=m+1;
end