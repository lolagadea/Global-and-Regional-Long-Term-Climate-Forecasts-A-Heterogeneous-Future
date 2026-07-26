function M=find_second_bestModel(BEAT)
a=max(BEAT(:,1));
pos=find(BEAT(:,1)==a);
for i=1:length(pos)
    if BEAT(pos(i),2)==0
        M(i)=pos(i);
    end
end

