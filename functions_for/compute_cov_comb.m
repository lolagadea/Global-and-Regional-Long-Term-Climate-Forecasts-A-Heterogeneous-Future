function var_comb=compute_cov_comb(L,Weig)
[t,m]=size(L);
pos=[];
%remove missing columns
for i=1:m
    if sum(isnan(L(:,i)))==t
        pos=[pos,i];
        m=m-1;
    end
    if m==0
        break
    end
end

    if m==0
        var_comb=NaN;
    elseif m==1
        var_comb=var(remove_NaN_col(L));
    else
L(:,pos)=[];
Weig(pos)=[];

if sum(sum(isnan(L)))==0
S=cov(L);
else
S=compute_cov_nan(L);
end
var_comb=Weig*diag(S);

if m>1
pair=combinator(m,2,'c');
n=size(pair,1);

for i=1:n
    var_comb=var_comb+2*Weig(pair(i,1))*Weig(pair(i,2))*S(pair(i,1),pair(i,2));
end
end
end

function B=compute_cov_nan(A)
[t,n]=size(A);
B=NaN(n,n);
for i=1:n
for j=1:n
if i==j
B(i,j)=nanvar(A(:,i));
else b=nancov(A(:,i),A(:,j));
B(i,j)=b(1,2);
end
end
end