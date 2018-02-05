function SRRF_analysis_code(reconmovnm,origmovnm,varargin)
%SRRF_ANALYSIS_CODE Point-and-click trace analysis
%
% possible updates: saving background, saving clip in structure,
% saving both original and reconstructed clips as files or in the structure

% Josh Ferguson
% Kural Lab
% The Ohio State University
% ferguson.621@osu.edu
switch nargin
    case 2
        tmpd = dir(reconmovnm);
        datafol = reconmovnm(1:end-length(tmpd.name));
        save_loc = [datafol filesep 'tracest.mat'];
    case 3
        save_loc = varargin{1};
end
global nsh gsh bsh afh
mod = [0,0,0];
ml = length(imfinfo(reconmovnm));
ff = 1; lf = ml;
ind = 0;
ss = size(imread(reconmovnm));
rxpos = NaN; rypos = NaN;
oxpos = NaN; oypos = NaN;
set(groot,'units','pixels');
screen = get(groot,'ScreenSize');
set(groot,'units','normalized');
rimg = zeros([ss ml],'uint16');
oimg = zeros([ss ml],'uint16');
for fr = 1:ml
    rimg(:,:,fr) = imread(reconmovnm,fr);
    oimg(:,:,fr) = imread(origmovnm,fr);
end
simg = sort(rimg(:),'ascend');
minrc = simg(ceil(0.0001*length(simg))); maxrc = simg(ceil(0.9999*length(simg)));
minoc = min(oimg(:)); maxoc = max(oimg(:));
imgy = .95;
imgx = ss(2)/ss(1)*imgy*screen(4)/screen(3);
close all
zimgy = 0.375;
zimgx = zimgy*screen(4)/screen(3);
graphy = imgy-2*zimgy;
graphx = (1-imgx)/2;
cfr = 1;
cintrad = 8;
ozrad = 12;
zrad = ozrad;
if exist(save_loc,'file')
    load_var = load(save_loc);
    tracest = load_var.tracest;
    disp(['Loaded file ' save_loc])
    ntrace = length(tracest);
else
    tracest = struct('frame',[],'xpos',[],'ypos',[],'int',[],'area',[],'ishot',false,'ispair',false,'mask',[]);
    ntrace = 0;
end

fh_img = figure(...
    'units','normalized',...
    'WindowStyle','modal',...
    'Position',[0, 1-imgy, imgx, imgy],...
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

fh_scroll = figure(...
    'units','normalized',...
    'OuterPosition',[0 0 1 1-imgy],...
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
    'OuterPosition',[imgx 1-imgy graphx graphy],...
    'MenuBar','none',...
    'ToolBar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Area',...
    'KeyPressFcn',@key_fun);
ah_area_graph = axes('units','normalized');

figure(...
    'units','normalized',...
    'OuterPosition',[imgx+graphx 1-imgy graphx graphy],...
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
        [area, int, SNR, srrfint] = deal(zeros(1,ml));
        zp = fh_img.CurrentPoint;
        tx = zp(1)*ss(2)+.5;
        ty = ss(1)-zp(2)*ss(1)+.5;
        if tx>zrad && tx<=ss(2)-zrad && ty>zrad && ty<=ss(1)-zrad
            [rxpos, rypos] = cofint(rimg,tx,ty,cfr);
        end
        ind = already_found(rxpos,rypos,cfr);
        if ind > 0
            mask = zeros([size(tracest(ind).mask(:,:,1)) ml]);
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
            mask = false(2*zrad+1, 2*zrad+1, ml);
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
                save(save_loc,'tracest')
            end
        catch
            save(save_loc,'tracest')
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
        axes(ah_img)
        imagesc(rimg(:,:,cfr),[minrc maxrc]);
        scatter_points(cfr)
        axis off
        axis equal
        if upz
            if src~=fh_img, mod(3) = -inf; end
            zoom_in;
            mod(3) = 0;
        end
        frame_line(ah_scroll,cfr,[.8 .8 .8])
    end
    function key_fun(~,event)
        %list of keys
        %escape, q,w,e, z,x,c, a,s, n, h, p, delete, backspace
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
        end
        if strcmp(event.Key,'l')
            key_test = 'l';
            while ~strcmp(key_test,'u')
                cp_tmp = get(groot,'PointerLocation');
                set(fh_scroll,'CurrentPoint',cp_tmp);
                move_callback(fh_scroll)
                pause(1/30)
                key_test = get(gcf,'CurrentCharacter');                
            end
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
            save(save_loc,'tracest');
            scatter_points(cfr);
        end
        if strcmp(event.Key,'p')
            if ind>0
                tracest(ind).ispair = true;
            else
                tracest(ntrace+1).ispair = true;
            end
            save(save_loc,'tracest');
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
            ind = 0;
            upz = false;
            rxpos = NaN; rypos = NaN;
            move_callback(fh_img)
        end
        if strcmp(event.Key,'g') || strcmp(event.Key,'b')
            if strcmp(event.Key,'g')
                goto_trace;
            elseif strcmp(event.Key,'b')
                goto_trace(ind+1);
            end
        end
        if strcmp(event.Key,'f')
            set(fh_img,'SelectionType','normal');
        end
    end
    function goto_trace(varargin)
        if nargin==0
            dh_goto = dialog(...
                'Units','Normalized',...
                'Position',[.4 .4 .2 .2],...
                'Name','Which trace?');
            uicontrol(...
                'Parent',dh_goto,...
                'Style','text',...
                'Units','Normalized',...
                'Position', [.3 .7 .4 .1],...
                'String','Which trace to go to?');
            qbox = uicontrol(...
                'Parent',dh_goto,...
                'Style','edit',...
                'Units','Normalized',...
                'Position',[.3 .3 .4 .2],...
                'Callback','set(gcf,''Windowstyle'',''normal'',''Visible'',''off'')',...
                'Selected','on');
            uicontrol(qbox)
            waitfor(dh_goto,'Visible','off')
            ind = str2double(qbox.String);
            close(dh_goto)
        else
            ind = varargin{1};
        end
        if ind>length(tracest)
            d = dialog(...
                'Units','Normalized',...
                'Position',[.4 .4 .2 .2],...
                'Name','Finished',...
                'KeyPressFcn',@key_fun);
            uicontrol('Parent',d,...
                'Units','Normalized',...
                'Style','text',...
                'Position',[.2 .5 .6 .1],...
                'String','Esc to close all, or click to continue.');
            uicontrol('Parent',d,...
                'Units','Normalized',...
                'Position',[.4 .2 .2 .1],...
                'String','Continue',...
                'Callback','delete(gcf)');
        else
            [area, int, SNR] = deal(zeros(1,ml));
            mask = false([size(tracest(ind).mask(:,:,1)), ml]);
            cfr = tracest(ind).frame(1);
            rxpos = tracest(ind).xpos(1);
            rypos = tracest(ind).ypos(1);
            mask(:,:,tracest(ind).frame) = tracest(ind).mask;
            tcfr = cfr;
            fff = tracest(ind).frame(1);
            flf = tracest(ind).frame(end);
            for ifr = fff:flf-1
                update_int(ifr,false)
                update_area(ifr,false)
            end
            update_int(flf,true)
            update_area(flf,true)
            cfr = fff;
            sff;
            cfr = flf;
            slf;
            cfr = tcfr;
            upz=true;
            mod(3) = -inf;
            move_callback(fh_img);
            mod = [0,0,0];
        end
    end
    function delete_trace(~,~)
        tracest(ind) = [];
        ntrace = length(tracest);
        scatter_points(cfr);
        save(save_loc,'tracest')
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
            %             line([zrad+1 zrad+1],[zrad-1 zrad+3],'linewidth',.5,'color','r')
            %             line([zrad-1 zrad+3],[zrad+1 zrad+1],'linewidth',.5,'color','r')
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
            %             line([zrad+1 zrad+1],[zrad-1 zrad+3],'linewidth',.5,'color','r')
            %             line([zrad-1 zrad+3],[zrad+1 zrad+1],'linewidth',.5,'color','r')
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
        tmpimg = double(img(floor(tmpy-zrad):floor(tmpy+zrad),...
            floor(tmpx-zrad):floor(tmpx+zrad),frame));
        tmp_mask = edge(tmpimg,'canny',[.5 .9])>0;
        tmp_mask = imfill(tmp_mask,'holes')>0;
        tmpimg = tmpimg.*tmp_mask;
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
        redo = false;
        if ind > 0
            if any(tracest(ind).frame==frame)
                area(frame) = tracest(ind).area(tracest(ind).frame==frame);
                done = true;
                if area(frame) == 0, redo = true; end
            end
        end
        if ~done || redo
            area(frame) = sum(sum(mask(:,:,frame)));
            if redo
                save(save_loc,'tracest')
            end
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
        redo = false;
        if ind > 0
            if any(tracest(ind).frame==frame)
                int(frame) = tracest(ind).int(tracest(ind).frame==frame);
                SNR(frame) = tracest(ind).SNR(tracest(ind).frame==frame);
                srrfint(frame) = tracest(ind).srrfint(tracest(ind).frame==frame);
                done = true;
                if int(frame)==0, redo = true; end
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
        if ~done || redo
            if redo, [oxpos, oypos] = cofint(oimg,rxpos,rypos,frame); end
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
            [int(frame),SNR(frame)] = twoDgaussianFitting_theta(tmp,false);
            tmp = double(rimg(floor(rypos-cintrad):floor(rypos+cintrad),...
                floor(rxpos-cintrad):floor(rxpos+cintrad),frame));
            tmp = interpolate_image(tmp);
            tmp = sort(tmp,'descend');
            srrfint(frame) = max(tmp(1:500));
            if redo
                tracest(ind).int(tracest(ind).frame==frame) = int(frame);
                tracest(ind).SNR(tracest(ind).frame==frame) = SNR(frame);
                tracest(ind).srrfint(tracest(ind).frame==frame) = srrfint(frame);
                save(save_loc,'tracest')
            end
        end
        if disp
            axes(ah_int_graph)
            if exist('iph','var'), delete(iph); end
            iph = plot(int);
        end
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
        up = double([mean(img(:))   1.1*(max(img(:))-min(img(:))) 3*c0(3)/2           3*c0(4)/2           2*pi size(img,2)/4 size(img,1)/4]);
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
        c = coeffvalues(gfit);
        SNR = c(2)/c(1);
        integ = quad2d(gfit,0.5,size(img,1)+.5,0.5,size(img,2)+.5)-c(1)*size(img,1)*size(img,2);
    end
    function int_img = interpolate_image(img)
        [xhave, yhave] = meshgrid(1:size(img,2),1:size(img,1));
        [xwant, ywant] = meshgrid(1:.1:size(img,2),1:.1:size(img,1));
        int_img = griddata(xhave(:),yhave(:),double(img(:)),xwant,ywant,'v4');
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
        if ind > 0
            if any(tracest(ind).frame==cfr)
                tracest(ind).area(tracest(ind).frame==cfr) = sum(sum(mask(:,:,cfr)));
                tracest(ind).mask(:,:,tracest(ind).frame==cfr) = mask(:,:,cfr);
            end
        end
        update_area(cfr,true);
        overlay_mask;
        cf_ball;
        save(save_loc,'tracest')
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
        tracest(spt).SNR = SNR(ff:lf);
        tracest(spt).area = area(ff:lf);
        tracest(spt).mask = mask(:,:,ff:lf);
        uih_saved = uicontrol('Parent',fh_text,...
            'Style','Text',...
            'FontSize',15,...
            'Units','Normalized',...
            'Position',[1/3 0 1/3 .1],...
            'String','Saved!');
        save(save_loc,'tracest')
        if ind==0
            ntrace = ntrace+1;
        end
        pause(.5)
        ind = already_found(rxpos,rypos,cfr);
        scatter_points(cfr)
        upz = false;
        rxpos = NaN; rypos = NaN;
        delete(uih_saved)
    end
end
