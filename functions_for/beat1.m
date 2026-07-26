function mod=beat1(X)
%by columns, 0 means is beaten, 1 beats
l=length(X);
mod=1;
for i=1:l
    if X(i)==0
        mod=NaN; break
    end
end

if sum(isnan(X))==l
    mod=NaN;
end