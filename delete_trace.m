function delete_trace(ind, sname)
load(sname,'tracest','mask')
for i = ind:length(tracest)-1
    mask(mask(:)==i) = 0;
    mask(mask(:)==i+1) = i;
end
tracest(ind) = [];
save(sname,'tracest','mask','-v7.3')
end