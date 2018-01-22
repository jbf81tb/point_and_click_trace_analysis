function SRRF_bead_analysis(reconmovnm,origmovnm)
%COMBINED_ANALYSIS_CODE Point-and-click trace analysis
%
% possible updates: saving background, saving clip in structure, 
% saving both original and reconstructed clips as files

% Josh Ferguson
% Kural Lab
% The Ohio State University
% ferguson.621@osu.edu

global nsh gsh bsh afh
mod = [0,0,0];
ml = length(imfinfo(reconmovnm));
ff = 1; lf = ml;
ind = 0;
ss = size(imread(reconmovnm));
rxpos = NaN; rypos = NaN;
oxpos = NaN; oypos = NaN;
screen = get(groot,'ScreenSize');
rimg = zeros([ss ml],'uint16');
oimg = zeros([ss ml],'uint16');
for fr = 1:ml
    rimg(:,:,fr) = imread(reconmovnm,fr);
    oimg(:,:,fr) = imread(origmovnm,fr);
end
simg = sort(rimg(:),'ascend');
minrc = min(rimg(:)); maxrc = simg(ceil(0.9999*length(simg)));
minoc = min(oimg(:)); maxoc = max(oimg(:));
% imgy = .95;
% imgx = ss(2)/ss(1)*imgy*screen(4)/screen(3);
imgx = .5;
imgy = ss(1)/ss(2)*imgx*screen(3)/screen(4);
dimgy = .95;
spacer = (dimgy-imgy)/2;
close all
zimgy = 0.375;
zimgx = zimgy*screen(4)/screen(3);
graphy = dimgy-2*zimgy;
graphx = (1-imgx)/2;
cfr = 1;
cintrad = 7;
ozrad = 13;
zrad = ozrad;
tmpd = dir(reconmovnm);
datafol = reconmovnm(1:end-length(tmpd.name));
% if ~exist([datafol, 'movies'],'dir'), mkdir([datafol, 'movies']); end
if exist([datafol filesep 'tracest.mat'],'file')
    load([datafol filesep 'tracest.mat'])
    ntrace = length(tracest);
else
    tracest = struct('frame',[],'xpos',[],'ypos',[],'int',[],'snr',[],'area',[],'ishot',false,'ispair',false,'mask',[]);
    ntrace = 0;
end

fh_img = figure(...
    'units','normalized',...
    'WindowStyle','modal',...
    'Position',[0, (1-dimgy)+spacer, imgx, imgy],...
    'SelectionType','alt',...
    'Menubar','none',...
    'Toolbar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'KeyPressFcn',@key_fun);
ah_img = axes('units','normalized',...
    'position',[0 0 1 1]);
imagesc(rimg(:,:,cfr),[minrc maxrc])
scatter_points(cfr)
axis equal
axis off

figure(...
    'units','normalized',...
    'OuterPosition',[0 0 1 1-dimgy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Pointer','crosshair',...
    'Resize','off',...
    'Name','Scroll',...
    'NumberTitle','off',...
    'WindowButtonMotionFcn',@move_callback,...
    'KeyPressFcn',@key_fun);
ah_scroll = axes('units','normalized','Position',[0 0 1 1]);
imagesc(ah_scroll,1:ml)
colormap('cool')
axis off

figure(...
    'units','normalized',...
    'OuterPosition',[imgx 1-zimgy zimgx zimgy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Zoomed Recon Image',...
    'KeyPressFcn',@key_fun,...
    'WindowButtonDownFcn',@select_mask);
ah_zoom_img = axes('units','normalized',...
    'position',[0 0 1 1]);

figure(...
    'units','normalized',...
    'OuterPosition',[imgx 1-2*zimgy zimgx zimgy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Original Image',...
    'KeyPressFcn',@key_fun);
ah_orig_img = axes('units','normalized',...
    'position',[0 0 1 1]);

figure(...
    'units','normalized',...
    'OuterPosition',[imgx 1-dimgy graphx graphy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Area',...
    'KeyPressFcn',@key_fun);
ah_area_graph = axes('units','normalized');

figure(...
    'units','normalized',...
    'OuterPosition',[imgx+graphx 1-dimgy graphx graphy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Intensity',...
    'KeyPressFcn',@key_fun);
ah_int_graph = axes('units','normalized');

fh_int_fit_graph = figure(...
    'units','normalized',...
    'OuterPosition',[imgx+zimgx 1-2*zimgy 1-imgx-zimgx zimgy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Intensity',...
    'KeyPressFcn',@key_fun);
ah_int_fit_graph = axes('units','normalized');


fh_text = figure(...
    'Units','Normalized',...
    'OuterPosition',[imgx+zimgx 1-zimgy 1-imgx-zimgx zimgy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Shortcuts',...
    'KeyPressFcn',@key_fun);
uicontrol('Parent',fh_text,...
    'Style','Text',...
    'Units','Normalized',...
    'Position',[1/4 .7 1/2 .15],...
    'String','(q,w,e): first frame, last frame, save')
uicontrol('Parent',fh_text,...
    'Style','Text',...
    'Units','Normalized',...
    'Position',[1/4 .5 1/2 .15],...
    'String','(a,s): back and forward frames')
uicontrol('Parent',fh_text,...
    'Style','Text',...
    'Units','Normalized',...
    'Position',[1/4 .3 1/2 .15],...
    'String','(z,x,c): zoom in, zoom out, return to normal')
uicontrol('Parent',fh_text,...
    'Style','Text',...
    'Units','Normalized',...
    'Position',[1/4 .1 1/2 .15],...
    'String','(n,h,p): show index number, declare hotspot, declare pair')

upz = false;
disp_er = true;
while true
    try
        waitfor(fh_img,'SelectionType','normal');
        [area, int, snr, srrfint] = deal(zeros(1,ml));
        zp = fh_img.CurrentPoint;
        tx = zp(1)*ss(2)+.5;
        ty = ss(1)-zp(2)*ss(1)+.5;
        if tx>zrad && tx<=ss(2)-zrad && ty>zrad && ty<=ss(1)-zrad
            [rxpos, rypos] = cofint(rimg,tx,ty,cfr);
        end
        ind = already_found(rxpos,rypos,cfr);
        mask = false(2*zrad+1, 2*zrad+1, ml);
        if ind > 0
            mask(:,:,tracest(ind).frame) = tracest(ind).mask;
            tcfr = cfr;
            fff = tracest(ind).frame(1);
            flf = tracest(ind).frame(end);
            for fr = fff:flf-1
                update_int(fr,false)
                update_area(fr,false)
            end
            update_int(flf,true)
            update_area(flf,true)
            cfr = fff;
            sff;
            cfr = flf;
            slf;
            cfr = tcfr;
        else
            tracest(ntrace+1).ishot = false; %#ok<*AGROW>
            tracest(ntrace+1).ispair = false;
        end
        upz = true;
        zoom_in;
        frame_line(ah_scroll,cfr,[0 0 0]);
        fh_img.SelectionType = 'alt';
    catch ME
        if disp_er
            close all
            rethrow(ME)
        end
        try
            if ~isempty(tracest(ntrace+1).frame)
                save([datafol filesep 'tracest.mat'],'tracest')
            end
        catch
            save([datafol filesep 'tracest.mat'],'tracest')
        end
        close all
        return;
    end
end

    function move_callback(src,~)
        if strcmp(src.Name,'Scroll') && strcmp(src.SelectionType,'normal')
            cp = src.CurrentPoint;
            cfr = ceil(cp(1)*ml);
        end
        if cfr<1
            cfr = 1;
        elseif cfr>ml
            cfr = ml;
        end
        ind = already_found(rxpos,rypos,cfr);
        axes(ah_img)
        imagesc(rimg(:,:,cfr),[minrc maxrc]);
        scatter_points(cfr)
        axis off
        axis equal
        if upz
            if src~=fh_img, mod(3) = -inf; end
            zoom_in;
            mod(3) = 0;
        elseif ind~=0
            cf_ball
        end
        frame_line(ah_scroll,cfr,[.8 .8 .8])
    end
    function key_fun(~,event)
        %list of keys
        %escape, q,w,e, z,x,c, a,s, n, h, p
        if strcmp(event.Key,'escape')
            disp_er = false;
            close all
        end
        if strcmp(event.Key,'q')
            sff;
        elseif strcmp(event.Key,'w')
            slf;
        elseif strcmp(event.Key,'e')
            save_trace;
        elseif strcmp(event.Key,'r')
            save_all;
        end
        if strcmp(event.Key,'z') || strcmp(event.Key,'x') || strcmp(event.Key,'c')
            if strcmp(event.Key,'z'), mod(3) = -1; end
            if strcmp(event.Key,'x'), mod(3) =  1; end
            if strcmp(event.Key,'c'), mod(3) = ozrad-zrad; end
            zrad = zrad+mod(3);
            zoom_in;
            mod = [0,0,0];
        end
        if strcmp(event.Key,'a') || strcmp(event.Key,'s')
            if strcmp(event.Key,'a'), mod(3) = -1; end
            if strcmp(event.Key,'s'), mod(3) =  1; end
            cfr = cfr+mod(3);
            move_callback(fh_img);
            mod = [0,0,0];
        end
        if strcmp(event.Key,'n')
            if exist('gsh','var'), delete(gsh); end
            if exist('bsh','var'), delete(bsh); end
            x = []; y = []; num = {};
            for i = 1:length(tracest)
                frind = find(tracest(i).frame==cfr);
                if isempty(frind), continue; end
                x(end+1) = tracest(i).xpos(frind);
                y(end+1) = tracest(i).ypos(frind);
                num{end+1} = num2str(i);
            end
            axes(ah_img)
            hold on
            nsh = text(x,y,num,'Color',[.7 0 0],'HorizontalAlignment','center');
            hold off
        end
        if strcmp(event.Key,'h')
            if ind>0
                tracest(ind).ishot = true;
            else
                tracest(ntrace+1).ishot = true;
            end
            scatter_points(cfr);
        end
        if strcmp(event.Key,'p')
            if ind>0
                tracest(ind).ispair = true;
            else
                tracest(ntrace+1).ispair = true;
            end
            scatter_points(cfr);
        end
        if strcmp(event.Key,'delete')
            if ind>0
                afh = dialog('Units','Normalized',...
                    'Position',[.3 .3 .4 .4]);
                uicontrol('Parent',afh,...
                    'Units','Normalized',...
                    'Position',[.2 .6 .6 .15],...
                    'Style','Text',...
                    'FontSize',14,...
                    'String','Are you sure you want to delete this trace?')
                uicontrol('Parent',afh,...
                    'Units','Normalized',...
                    'Position',[.3 .35 .1 .1],...
                    'String','Yes',...
                    'Callback',@delete_trace);
                uicontrol('Parent',afh,...
                    'Units','Normalized',...
                    'Position',[.6 .35 .1 .1],...
                    'String','No',...
                    'Callback','close(gcf)');
                waitfor(afh)
            end
        end
        if strcmp(event.Key,'backspace')
            upz = false;
            rxpos = NaN; rypos = NaN;
            move_callback(fh_img)
        end
        if strcmp(event.Key,'i')
            [fname,fpath] = uigetfile('*.mat','Select the MAT file with x and y variables');
            try
                load([fpath filesep fname])
                for ix = 1:length(x)
                    rxpos = x(ix);
                    oxpos = x(ix);
                    rypos = y(ix);
                    oypos = y(ix);
                    save_all;
                end
            catch me
                close all
                rethrow(me)
            end
        end
                    
            
    end
    function delete_trace(~,~)
        tracest(ind) = [];
        ntrace = length(tracest);
        scatter_points(cfr);
        save([datafol filesep 'tracest.mat'],'tracest')
%         filenm = [datafol, 'movies', filesep, sprintf('%04u',ind),'.tif'];
%         if exist(filenm,'file'), delete(filenm); end
        ind = 0;
        close(afh)
    end

    function zoom_in
        persistent lh
        if all(mod==0)
            zp = fh_img.CurrentPoint;
            tmpx = zp(1)*ss(2)+.5;
            tmpy = ss(1)-zp(2)*ss(1)+.5;
            if tmpx>zrad && tmpx<=ss(2)-zrad && tmpy>zrad && tmpy<=ss(1)-zrad
                [rxpos, rypos] = cofint(rimg,tmpx,tmpy,cfr);
            end
        elseif mod(3)==0
            rxpos = rxpos + mod(1);
            rypos = rypos + mod(2);
        else
            if rxpos>zrad && rxpos<=ss(2)-zrad && rypos>zrad && rypos<=ss(1)-zrad
                [rxpos, rypos] = cofint(rimg,rxpos,rypos,cfr);
            end
        end
        if rxpos>zrad && rxpos<=ss(2)-zrad && rypos>zrad && rypos<=ss(1)-zrad
            axes(ah_zoom_img)
            zimg = rimg(floor(rypos-zrad):floor(rypos+zrad),...
                floor(rxpos-zrad):floor(rxpos+zrad),cfr);
            imagesc(zimg,[minrc maxrc])
            hold on
            line([zrad+1 zrad+1],[zrad-1 zrad+3],'linewidth',.5,'color','r')
            line([zrad-1 zrad+3],[zrad+1 zrad+1],'linewidth',.5,'color','r')
            line([zrad+.5-cintrad, zrad+1.5+cintrad, zrad+1.5+cintrad, zrad+0.5-cintrad, zrad+.5-cintrad],...
                [zrad+.5-cintrad, zrad+0.5-cintrad, zrad+1.5+cintrad, zrad+1.5+cintrad, zrad+.5-cintrad],...
                'color','k','linewidth',1);
            axis off
            hold off
            
            [oxpos, oypos] = cofint(oimg,rxpos,rypos,cfr);
            axes(ah_orig_img)
            oimgc = oimg(floor(oypos-zrad):floor(oypos+zrad),...
                floor(oxpos-zrad):floor(oxpos+zrad),cfr);
            imagesc(oimgc,[minoc maxoc]);
            hold on
            line([zrad+1 zrad+1],[zrad-1 zrad+3],'linewidth',.5,'color','r')
            line([zrad-1 zrad+3],[zrad+1 zrad+1],'linewidth',.5,'color','r')
            line([zrad+.5-cintrad, zrad+1.5+cintrad, zrad+1.5+cintrad, zrad+0.5-cintrad, zrad+.5-cintrad],...
                [zrad+.5-cintrad, zrad+0.5-cintrad, zrad+1.5+cintrad, zrad+1.5+cintrad, zrad+.5-cintrad],...
                'color','k','linewidth',1);
            axis off
            hold off
            
            axes(ah_img)
            hold on
            if exist('lh','var'), delete(lh); end
            lh = line([rxpos-zrad rxpos+zrad rxpos+zrad rxpos-zrad rxpos-zrad],...
                [rypos-zrad rypos-zrad rypos+zrad rypos+zrad rypos-zrad],...
                'color','r','linewidth',1);
            axis off
            hold off
            
            sizedif = size(mask,1)-size(zimg,1);
            if sizedif == 0
                if sum(sum(mask(:,:,cfr)))==0
                    initialize_area(cfr)
                end
            elseif sizedif > 0
                while sizedif > 0
                    mask = mask(2:end-1,2:end-1,:);
                    sizedif = sizedif-2;
                end
            elseif sizedif < 0
                while sizedif < 0
                    mask(end+1:end+2,:,:) = zeros(2,size(mask,2),size(mask,3));
                    mask(:,end+1:end+2,:) = zeros(size(mask,1),2,size(mask,3));
                    mask(2:end-1,2:end-1,:) = mask(1:end-2,1:end-2,:);
                    sizedif = sizedif+2;
                end
            end
            overlay_mask
            update_area(cfr,true)
            update_int(cfr,true)
            cf_ball
        end
    end
    function overlay_mask
        persistent msh
        axes(ah_zoom_img)
        delete(msh)
        hold on
        [iy,ix] = find(mask(:,:,cfr));
        msh = scatter(ix,iy,50,[1 0 1],'linewidth',1);
        axis off
        hold off
    end
    function [cx, cy] = cofint(img,tmpx,tmpy,frame)
        img = double(img(floor(tmpy-zrad):floor(tmpy+zrad),...
            floor(tmpx-zrad):floor(tmpx+zrad),frame));
        cx = ceil(size(img,2)/2);
        cy = ceil(size(img,1)/2);
        maxi = 0;
        for r = 0:ceil(sqrt(cx^2+cy^2))
            [x, y] = meshgrid(1:size(img,2),1:size(img,1));
            x = x-cx;
            y = y-cy;
            msk = round(sqrt(x.^2+y.^2))==r;
            if sum(msk(:))==0, break; end
            tmpimg = msk.*img;
            if max(tmpimg(:))>maxi
                maxi = max(tmpimg(:));
                [my,mx] = find(tmpimg==maxi,1);
            end
        end
        saved = false(size(img));
        cx = mx;
        cy = my;
        saved(cy,cx) = true;
        tmpimg = img;
        nx = cx; ny = cy;
        while true
            for i = -1:1
                for j = -1:1
                    if i==0 && j==0, continue; end
                    if cy+j>size(saved,1) || cy+j<=0, continue; end
                    if cx+i>size(saved,2) || cx+i<=0, continue; end
                    if saved(cy+j,cx+i), continue; end
                    try
                        dif = tmpimg(cy,cx)-tmpimg(cy+j,cx+i);
                        if dif>=0
                            nx(end+1) = cx+i;
                            ny(end+1) = cy+j;
                        end
                    catch
                    end
                end
            end
            saved(cy,cx) = true;
            rm_ind = nx==cx&ny==cy;
            nx(rm_ind) = [];
            ny(rm_ind) = [];
            if isempty(nx) && isempty(ny), break; end
            cx = nx(end);
            cy = ny(end);
        end
        tmp_mask = imfill(saved,'holes');
        
        tmpimg = img.*tmp_mask;
        if sum(tmpimg(:))==0
            cx = floor(tmpx);
            cy = floor(tmpy);
        else
            cx = floor(tmpx-zrad) + sum(tmpimg*(1:size(tmpimg,2))')/sum(tmpimg(:));
            cy = floor(tmpy-zrad) + sum((1:size(tmpimg,1))*tmpimg)/sum(tmpimg(:));
        end
    end
    function initialize_area(frame)
        tmp = double(rimg(floor(rypos-zrad):floor(rypos+zrad),...
            floor(rxpos-zrad):floor(rxpos+zrad),frame));
        mask(:,:,frame) = edge(tmp,'canny',[.5 .9])>0;
        mask(:,:,frame) = imfill(mask(:,:,frame),'holes')>0;
    end
    function update_area(frame,disp)
        persistent aph
        done = false;
        if ind > 0
            if any(tracest(ind).frame==frame)
                area(frame) = tracest(ind).area(tracest(ind).frame==frame);
                done = true;
            end
        end
        if ~done
            area(frame) = sum(sum(mask(:,:,frame)));
        end
        if disp
            axes(ah_area_graph)
            if exist('aph','var'), delete(aph); end
            aph = plot(area);
        end
    end
    function update_int(frame,disp)
        persistent iph
        done = false;
        if ind > 0
            if any(tracest(ind).frame==frame)
                int(frame) = tracest(ind).int(tracest(ind).frame==frame);
                snr(frame) = tracest(ind).snr(tracest(ind).frame==frame);
                srrfint(frame) = tracest(ind).srrfint(tracest(ind).frame==frame);
                done = true;
                ifgc = get(fh_int_fit_graph,'Children');
                for i = 1:length(ifgc)
                    ihgc = get(ifgc(i),'Children');
                    for j = 1:length(ihgc)
                        set(ihgc(j),'visible','off')
                    end
                    set(ifgc(i),'visible','off')
                end
            end
        end
        if ~done
            ifgc = get(fh_int_fit_graph,'Children');
            for i = 1:length(ifgc)
                ihgc = get(ifgc(i),'Children');
                for j = 1:length(ihgc)
                    set(ihgc(j),'visible','on')
                end
                set(ifgc(i),'visible','on')
            end
            tmp = double(oimg(floor(oypos-zrad):floor(oypos+zrad),...
                floor(oxpos-zrad):floor(oxpos+zrad),frame));
            [int(frame), snr(frame)] = twoDgaussianFitting_theta(tmp,false);
            tmp = double(rimg(floor(rypos-cintrad):floor(rypos+cintrad),...
                floor(rxpos-cintrad):floor(rxpos+cintrad),frame));
            tmp = interpolate_image(tmp);
            tmp = sort(tmp,'descend');
            srrfint(frame) = max(tmp(1:500));
        end
        if disp
            axes(ah_int_graph)
            if exist('iph','var'), delete(iph); end
            iph = plot(int);
        end
    end
    function int_img = interpolate_image(img)
        [xhave, yhave] = meshgrid(1:size(img,2),1:size(img,1));
        [xwant, ywant] = meshgrid(1:.1:size(img,2),1:.1:size(img,1));
        int_img = griddata(xhave(:),yhave(:),double(img(:)),xwant,ywant,'v4');
    end
    function [integ, SNR] = twoDgaussianFitting_theta(img, disp)
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
        up = double([mean(img(:))   1.1*(max(img(:))-min(img(:))) 3*c0(3)/2           3*c0(4)/2           2*pi size(img,2)/2 size(img,1)/2]);
        xdata = double(xdata); ydata = double(ydata); img = double(img);
        gfit = fit([xdata(:), ydata(:)], img(:), F, 'StartPoint', c0, 'Lower', low, 'Upper', up);
        if disp
            fh_gauss_fit = figure;
            plot(gfit, [xdata(:), ydata(:)], img(:)); 
            waitfor(fh_gauss_fit,'SelectionType','alt')
            close(fh_gauss_fit)
        end
        axes(ah_int_fit_graph)
        plot(gfit, [xdata(:), ydata(:)], img(:));
%         rotate3d on
%         title('Click another figure before pressing a key')
        c = coeffvalues(gfit);
        SNR = c(2)/c(1);
        integ = quad2d(gfit,0.5,size(img,1)+.5,0.5,size(img,2)+.5)-c(1)*size(img,1)*size(img,2);
    end
    function scatter_points(frame)
        if exist('nsh','var'), delete(nsh); end
        if exist('bsh','var'), delete(bsh); end
        if exist('gsh','var'), delete(gsh); end
        xb = []; xg = [];
        yb = []; yg = [];
        for i = 1:length(tracest)
            frind = find(frame==tracest(i).frame);
            if isempty(frind), continue; end
            if tracest(i).ishot || tracest(i).ispair
                xb(end+1) = tracest(i).xpos(frind);
                yb(end+1) = tracest(i).ypos(frind);
            else
                xg(end+1) = tracest(i).xpos(frind);
                yg(end+1) = tracest(i).ypos(frind);
            end
        end
        axes(ah_img)
        hold on
        bsh = scatter(xb,yb,100,'r','x');
        gsh = scatter(xg,yg,100,'g','o');
        hold off
    end
    function cf_ball
        persistent abh ibh
        axes(ah_area_graph)
        if exist('abh','var'), delete(abh); end
        hold on
        abh = plot(cfr,area(cfr),'ro');
        hold off
        axes(ah_int_graph)
        if exist('ibh','var'), delete(ibh); end
        hold on
        ibh = plot(cfr,int(cfr),'ro');
        hold off
    end
    function frame_line(src,loc,col)
        axes(src)
        tmpch = get(src,'Children');
        for i = 1:length(tmpch)
            if strcmp(tmpch(i).Type,'line')
                if all(tmpch(i).Color==col)
                    delete(tmpch(i))
                end
            end
        end
        hold on
        line([loc loc],[.5 1.5],'linewidth',1,'color',col);
        hold off
    end
    function select_mask(src,~)
        cp = src.CurrentPoint;
        cp(2) = 1-cp(2);
        cp = round(cp([2,1]).*(2*zrad+1)+.5);
        if strcmp(src.SelectionType,'normal')
            mask(cp(1),cp(2),cfr) = true;
        elseif strcmp(src.SelectionType,'alt')
            mask(cp(1),cp(2),cfr) = false;
        end
        mask(:,:,cfr) = imfill(mask(:,:,cfr),'holes')>0;
        update_area(cfr,true);
        overlay_mask;
        cf_ball;
    end
    function ind = already_found(xp,yp,frame)
        persistent uih_found
        ind = 0;
        delete(uih_found)
        for i = 1:length(tracest)
            frind = find(tracest(i).frame==frame);
            if isempty(frind), continue; end
            dist = sqrt((tracest(i).xpos(frind)-xp)^2+...
                (tracest(i).ypos(frind)-yp)^2);
            if dist<cintrad
                ind = i;
                uih_found = uicontrol('Parent',fh_text,...
                'Style','Text',...
                'FontSize',15,...
                'Units','Normalized',...
                'Position',[1/3 0 1/3 .1],...
                'String',sprintf('Index = %u',ind));
                return;
            end
        end
    end
    function [xf, yf] = get_positions(xp,yp,first,last)
        [xf,yf] = deal(zeros(1,ml));
        if ind>0
            for iframe = 1:length(tracest(ind).frame)
                xf(tracest(ind).frame(iframe)) = tracest(ind).xpos(iframe);
                yf(tracest(ind).frame(iframe)) = tracest(ind).ypos(iframe);
            end
        end
        i = cfr;
        if xf(i)==0 && yf(i)==0
            [xf(i),yf(i)] = cofint(rimg,xp,yp,i);
        end
        while i<last
            i = i+1;
            if xf(i)==0 && yf(i)==0
                [xf(i),yf(i)] = cofint(rimg,xf(i-1),yf(i-1),i);
            end
        end
        i = cfr;
        while i>first
            i = i-1;
            if xf(i)==0 && yf(i)==0
                [xf(i),yf(i)] = cofint(rimg,xf(i+1),yf(i+1),i);
            end
        end
    end
    function sff(~,~)
        ff = cfr;
        frame_line(ah_scroll,cfr,[0 .5 0])
    end
    function slf(~,~)
        lf = cfr;
        frame_line(ah_scroll,cfr,[.7 0 0])
    end
    function save_trace(~,~)
        [xot,yot] = get_positions(rxpos,rypos,ff,lf);
        if ind>0
            spt = ind;
        else
            spt = ntrace+1;
        end
            tracest(spt).frame = ff:lf;
            tracest(spt).xpos = xot(ff:lf);
            tracest(spt).ypos = yot(ff:lf);
            tracest(spt).int = int(ff:lf);
            tracest(spt).srrfint = srrfint(ff:lf);
            tracest(spt).area = area(ff:lf);
            tracest(spt).mask = mask(:,:,ff:lf);
            tracest(spt).snr = snr(ff:lf);
            uih_saved = uicontrol('Parent',fh_text,...
                'Style','Text',...
                'FontSize',15,...
                'Units','Normalized',...
                'Position',[1/3 0 1/3 .1],...
                'String','Saved!');
            save([datafol filesep 'tracest.mat'],'tracest')
            if ind==0
                ntrace = ntrace+1;
            end
%             ind = already_found(rxpos,rypos,cfr);
            ind = 0;
            scatter_points(cfr)
            upz = false;
            rxpos = NaN; rypos = NaN;
            pause(.5)
            delete(uih_saved)
    end
    function save_all(~,~)
        ff = 1;
        lf = ml;
        [xot,yot] = get_positions(rxpos,rypos,ff,lf);
        if ind>0
            spt = ind;
        else
            spt = ntrace+1;
        end
        for i = 1:ml
            initialize_area(i);
            update_area(i,false);
            update_int(i,false);
        end
        tracest(spt).frame = ff:lf;
        tracest(spt).xpos = xot(ff:lf);
        tracest(spt).ypos = yot(ff:lf);
        tracest(spt).int = int(ff:lf);
        tracest(spt).srrfint = srrfint(ff:lf);
        tracest(spt).area = area(ff:lf);
        tracest(spt).mask = mask(:,:,ff:lf);
        tracest(spt).snr = snr(ff:lf);
        tracest(spt).ishot = false;
        tracest(spt).ispair = false;
        uih_saved = uicontrol('Parent',fh_text,...
            'Style','Text',...
            'FontSize',15,...
            'Units','Normalized',...
            'Position',[1/3 0 1/3 .1],...
            'String','Saved!');
        save([datafol filesep 'tracest.mat'],'tracest')
        if ind==0
            ntrace = ntrace+1;
        end
%         ind = already_found(rxpos,rypos,cfr);
        ind = 0;
        scatter_points(cfr)
        upz = false;
        rxpos = NaN; rypos = NaN;
        pause(.5)
        delete(uih_saved)
    end
end
