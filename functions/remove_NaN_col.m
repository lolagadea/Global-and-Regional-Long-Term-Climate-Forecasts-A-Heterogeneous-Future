function [Y,ind_c]=remove_NaN_col(X)
[index_r, index_c]=find(isnan(X));
index_c=unique(index_c);
n=size(X,2);
ind_c=setdiff(1:n,index_c);
Y=X(:,ind_c);