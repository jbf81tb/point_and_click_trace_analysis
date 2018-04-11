function from488to560(fname)
    loadvar = load(fname,'tracest');
    tracest = loadvar.tracest;
%     mask = loadvar.mask;
    clear loadvar
    for i = 1:length(tracest)
            tracest(i).xpos = tracest(i).xpos-9;
            tracest(i).ypos = tracest(i).ypos-4;
%             for fr = tracest(i).frame
%                 [y,x] = find(mask(:,:,fr)==i);
%                 for j = 1:length(y)
%                     mask(y(j),x(j),fr) = 0;
%                     mask(y(j)-2,x(j)-8,fr) = i;
%                 end
%             end
    end
    ind = strfind(fname,'488');
    save([fname(1:ind-1) '560' fname(ind+3:end)],'tracest')
end