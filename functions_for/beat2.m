function mod=beat2(X)
%By rows, 0 means beats, 1 is beaten
l=length(X);
mod=0;
for i=1:l
    if X(i)==1
        mod=NaN; break
    end
end

if sum(isnan(X))==l
    mod=NaN;
end