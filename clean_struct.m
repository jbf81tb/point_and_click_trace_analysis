function clean_struct(sname)
tmp = load(sname);
tracest = tmp.tracest;
ismask = false;
if isfield(tmp,'mask')
    mask = tmp.mask;
    ismask = true;
end
for i = length(tracest):-1:1
    if isempty(tracest(i).frame)
        if ismask
            for j = i:length(tracest)-1
                mask(mask(:)==j) = 0;
                mask(mask(:)==j+1) = j;
            end
        end
    tracest(i) = [];
    continue;
    end
    if any(tracest(i).xpos==0)
        miss = find(tracest(i).xpos)==0;
        have = find(tracest(i).xpos);
        for j = miss
            [~,nearest] = min(abs(j-have));
            tracest(i).xpos(j) = tracest(i).xpos(have(nearest));
            tracest(i).ypos(j) = tracest(i).ypos(have(nearest));
        end
    end
end
save(sname,'tracest','mask')
end