function A2=array_to_beat_matrix(a,l)
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
A2=A;
for i=1:l-1
    for j=i+1:l
        if A(i,j)==1
            A2(j,i)=0;
        elseif A(i,j)==0
            A2(j,i)=1;
        end
    end
end