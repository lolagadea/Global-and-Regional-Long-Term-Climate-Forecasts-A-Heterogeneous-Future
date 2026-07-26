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
