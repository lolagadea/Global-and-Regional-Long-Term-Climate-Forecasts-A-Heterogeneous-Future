function [Y,ind_r]=remove_NaN_row(X)
[index_r, index_c]=find(isnan(X));
index_r=unique(index_r);
n=size(X,1);
ind_r=setdiff(1:n,index_r);
Y=X(ind_r,:);