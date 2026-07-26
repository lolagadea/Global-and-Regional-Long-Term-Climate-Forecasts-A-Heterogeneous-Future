function z=mult_matrix_nan(a,w)
n=size(a,2);
if size(a,2)~=size(w,1)
    disp('Dimensions are no consistent')
    quit
end

if sum(isnan(a))==n
    z=NaN;
else
[t]=length(a);
A=[a',w];
A2=remove_NaN(A);
a2=A2(:,1);
w2=A2(:,2);
w2=w2/sum(w2); %rescale w
z=a2'*w2;
end