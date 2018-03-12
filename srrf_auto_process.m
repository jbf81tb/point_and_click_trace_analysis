function srrf_auto_process(reconmovnm,origmovnm,varargin)
% SRRF_AUTO_PROCESS: take results from cmeAnalysis and reprocess
%
% Josh Ferguson
% Kural Lab
% The Ohio State University
% ferguson.621@osu.edu
% https://github.com/jbf81tb/point_and_click_trace_analysis
switch nargin
    case 2
        tmpd = dir(reconmovnm);
        datafol = reconmovnm(1:end-(length(tmpd.name)+1));
        save_loc = [datafol filesep 'sm_tracest.mat'];
    case 3
        save_loc = varargin{1};
end
ml = length(imfinfo(reconmovnm));
ss = size(imread(reconmovnm));
rimg = zeros([ss ml],'uint16');
oimg = zeros([ss ml],'uint16');
for fr = 1:ml
    rimg(:,:,fr) = imread(reconmovnm,fr);
    oimg(:,:,fr) = imread(origmovnm,fr);
end
zrad = 5;
load_var = load(save_loc);
tracest = load_var.tracest;
disp(['Loaded file ' save_loc])
tic
ntrace = length(tracest);
for ind = 1:ntrace
    for ifr = 1:length(tracest(ind).frame)
        if tracest(ind).int(ifr)==0
            xpos = tracest(ind).xpos(ifr);
            ypos = tracest(ind).ypos(ifr);
            frame = tracest(ind).frame(ifr);
            timg = double(oimg(floor(ypos-zrad):floor(ypos+zrad),...
                floor(xpos-zrad):floor(xpos+zrad),frame));
            [tracest(ind).int(ifr),tracest(ind).SNR(ifr)] = twoDgaussianFitting_theta(timg);
            timg = double(rimg(floor(ypos-zrad):floor(ypos+zrad),...
                floor(xpos-zrad):floor(xpos+zrad),frame));
            timg = interpolate_image(timg);
            timg = sort(timg,'descend');
            tracest(ind).srrfint(ifr) = sum(double(timg(1:100)));
        end
    end
end
save(save_loc,'tracest')
toc
return

    function [integ, SNR] = twoDgaussianFitting_theta(img)
        % c(1) = background
        % c(2) = amplitude
        % c(3) = x center
        % c(4) = y center
        % c(5) = theta about z axis
        % c(6) = sd_x
        % c(7) = sd_y
        [ydata, xdata] = meshgrid(1:size(img,2), 1:size(img,1));
        F = @(back, amp, x0, y0, th, sx, sy, x, y)back+amp*exp(-( ...
            (cos(th)^2/(2*sx^2)+sin(th)^2/(2*sy^2))*(x-x0).^2 + ...
            (sin(2*th)/(4*sy^2)-sin(2*th)/(4*sx^2))*(x-x0).*(y-y0) + ...
            (sin(th)^2/(2*sx^2)+cos(th)^2/(2*sy^2))*(y-y0).^2));
        c0 = double([mean(min(img)) max(img(:))-mean(min(img))    ceil(size(img,2)/2) ceil(size(img,1)/2) pi   1             1]);
        low = double([min(img(:))   mean(img(:))-min(img(:))      c0(3)/2             c0(4)/2             0    0             0]);
        up = double([mean(img(:))   1.1*(max(img(:))-min(img(:))) 3*c0(3)/2           3*c0(4)/2           2*pi size(img,2)/4 size(img,1)/4]);
        xdata = double(xdata); ydata = double(ydata); img = double(img);
        gfit = fit([xdata(:), ydata(:)], img(:), F, 'StartPoint', c0, 'Lower', low, 'Upper', up);
        c = coeffvalues(gfit);
        SNR = c(2)/c(1);
        integ = quad2d(gfit,0.5,size(img,1)+.5,0.5,size(img,2)+.5)-c(1)*size(img,1)*size(img,2);
    end
    function int_img = interpolate_image(img)
        [xhave, yhave] = meshgrid(1:size(img,2),1:size(img,1));
        [xwant, ywant] = meshgrid(1:.1:size(img,2),1:.1:size(img,1));
        int_img = griddata(xhave(:),yhave(:),double(img(:)),xwant,ywant,'v4');
    end
end