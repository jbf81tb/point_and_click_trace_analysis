function pcta(reconmovnm,origmovnm,varargin)
% PCTA: Point-and-click trace analysis
% pcta(reconstructed_movie, original_movie, save_file = './tracest.mat',...
%      reconstruction_type = 'sim', movie_type = 'tif'); 
%
% TODO: saving clips and linking clips to mask
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
        save_loc = [datafol filesep 'tracest.mat'];
        type = 'sim';
        ftype = 'tif';
    case 3
        save_loc = varargin{1};
        type = 'sim';
        ftype = 'tif';
    case 4
        save_loc = varargin{1};
        type = varargin{2};
        ftype = 'tif';
    case 5
        save_loc = varargin{1};
        type = varargin{2};
        ftype = varargin{3};
end
global nsh gsh bsh afh
ind = 0;
rxpos = NaN; rypos = NaN;
oxpos = NaN; oypos = NaN;
zoomQuad = true(2);
set(groot,'units','pixels');
screen = get(groot,'ScreenSize');
set(groot,'units','normalized');
if strcmp(ftype,'tif')
    ml = length(imfinfo(reconmovnm));
    ff = 1; lf = ml;
    ss = size(imread(reconmovnm));
    rimg = zeros([ss ml],'single');
    oimg = zeros([ss ml],'uint16');
    for fr = 1:ml
        rimg(:,:,fr) = imread(reconmovnm,fr);
        oimg(:,:,fr) = imread(origmovnm,fr);
    end
elseif strcmp(ftype,'mrc')
    [rimg,s] = ReadMRC(reconmovnm);
    ss = [s.ny, s.nx];
    ml = s.nz;
    oimg = zeros([ss ml],'uint16');
    placeholder = ReadMRC(origmovnm);
    ti = 1;
    while ti<=ml
        oimg(:,:,ti) = imresize(mean(placeholder(:,:,1:9),3),2,'bicubic');
        ti = ti+1;
        placeholder(:,:,1:9) = [];
    end
end

simg = sort(rimg(:),'ascend');
minrc = simg(ceil(0.0001*length(simg)));
medrc = simg(ceil(0.4*length(simg)));
maxrc = simg(ceil(0.9999*length(simg)));
% minoc = min(oimg(:)); maxoc = max(oimg(:));
imgy = .95;
imgx = ss(2)/ss(1)*imgy*screen(4)/screen(3);
close all
zimgy = 0.375;
zimgx = zimgy*screen(4)/screen(3);
graphy = imgy-2*zimgy;
graphx = (1-imgx)/2;
cfr = 1;
ozrad = 7;
zrad = ozrad;
if exist(save_loc,'file')
    load_var = load(save_loc);
    tracest = load_var.tracest;
    if isfield(load_var,'mask')
        mask = load_var.mask;
    else
        mask = zeros([ss ml],'uint16');
        save(save_loc,'tracest','mask','-v7.3');
    end
    disp(['Loaded file ' save_loc])
    ntrace = length(tracest);
else
    tracest = struct('frame',[],'xpos',[],'ypos',[],'int',[],'area',[],'ishot',false,'ispair',false);
    ntrace = 0;
    mask = zeros([ss ml],'uint16');
    save(save_loc,'tracest','mask','-v7.3');
end
mh = matfile(save_loc,'Writable',true);

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
image(ah_img,rimg(:,:,cfr),'CDataMapping','scaled')
set(ah_img,'CLim',[medrc maxrc],...
    'Xtick',[],'YTick',[]);
scatter_points(cfr)

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
set(ah_scroll,'XTick',[],'YTick',[])

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
    'Position',[0.05 .05 .9 .9],...
    'HorizontalAlignment','left',...
    'Fontsize',14,...
    'Fontweight','bold',...
    'String',{'(esc): Exit the program.',...
    '(q,w,e): Set first frame or last frame, then save the trace.',...
    '(a,s): Go back or forward by 1 frame.',...
    '(z,x): Zoom in or zoom out.',...
    '(h,p): Declare hotspot or declare pair.',...
    '(l,u): Lock or unlock frame scrolling.',...
    '(g,c): Goto a trace or goto the next trace.',...
    '(delete): Delete a trace.',...
    '(backspace): Disable the zoom window.',...
    '(4): Reduce the image to a single quadrant.'})
upz = false;
disp_er = true;
while true
    try
        waitfor(fh_img,'SelectionType','normal');
        [area, int, SNR, srrfint, xpos, ypos] = deal(zeros(1,ml));
        zp = fh_img.CurrentPoint;
        if all(zoomQuad(:))
            rxpos = ss(2)*zp(1)+.5;
            rypos = ss(1)*(1-zp(2))+.5;
        else
            rxpos = ss(2)/2*zp(1)+.5+ss(2)*sum(zoomQuad(:,2))/2;
            rypos = ss(1)/2*(1-zp(2))+.5+ss(1)*sum(zoomQuad(2,:))/2;
        end
        if rxpos>zrad && rxpos<=ss(2)-zrad && rypos>zrad && rypos<=ss(1)-zrad
            [rxpos, rypos] = cofint(rimg,rxpos,rypos,cfr);
            xpos(cfr) = rxpos;
            ypos(cfr) = rypos;
        end
        ind = already_found(rxpos,rypos,cfr);
        if ind > 0
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
            mask(mask==bitcmp(0,'uint16')) = 0;
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
        try tracest(ntrace+1) = []; catch, end
        mask(mask==bitcmp(0,'uint16')) = 0;
        save(save_loc,'tracest','mask')
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
        image(ah_img,rimg(:,:,cfr),'CDataMapping','scaled');
        scatter_points(cfr)
        set(ah_img,...
            'XLim',ceil([0, ss(2)]*(sum(zoomQuad(:))/6+1/3))+~any(zoomQuad(:,1))*floor(ss(2)/2)+.5,...
            'YLim',ceil([0, ss(1)]*(sum(zoomQuad(:))/6+1/3))+~any(zoomQuad(1,:))*floor(ss(1)/2)+.5,...
            'CLim',[medrc maxrc],...
            'Xtick',[],'YTick',[]);
        if upz
            zoom_in;
        end
        frame_line(ah_scroll,cfr,[.8 .8 .8])
    end
    function key_fun(~,event)
        %%%%% List of keys %%%%%
        %escape - save and escape from program
        %q,w,e - set frame window
        %l,u - Lock or Unlock mouse scrolling
        %z,x - Zoom level
        %a,s - individual frame movement
        %n - show trace Numbers
        %f - click in place
        %h,p - declare Hotspots and Pairs
        %delete - delete trace from structure
        %backspace - stop zooming and graphing
        %g,c - goto trace and next trace
        %k - save structure and mask ('Keep')
        %4 - select a quadrant
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
            cp_tmp = get(groot,'PointerLocation');
            set(groot,'PointerLocation',[cfr/ml cp_tmp(2)]);
            key_test = 'l';
            while ~strcmp(key_test,'u')
                try
                    cp_tmp = get(groot,'PointerLocation');
                    set(fh_scroll,'CurrentPoint',cp_tmp);
                    move_callback(fh_scroll)
                    pause(1/30)
                    key_test = lower(get(gcf,'CurrentCharacter'));
                catch
                    close all
                    return
                end
            end
        end
        if strcmp(event.Key,'z') || strcmp(event.Key,'x')
            if strcmp(event.Key,'z'), zrad = zrad-1; end
            if strcmp(event.Key,'x'), zrad = zrad+1; end
            zoom_in;
        end
        if strcmp(event.Key,'a') || strcmp(event.Key,'s')
            if strcmp(event.Key,'a'), cfr = cfr-1; end
            if strcmp(event.Key,'s'), cfr = cfr+1; end
            move_callback(fh_img);
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
            set(ah_img,'Nextplot','add')
            nsh = text(ah_img,x,y,num,...
                'Color',[.8 .8 .8],...
                'HorizontalAlignment','center',...
                'FontWeight','bold');
            set(ah_img,'Nextplot','replace')
        end
        if strcmp(event.Key,'f')
            upz = true;
            zoom_in
        end
        if strcmp(event.Key,'h')
            if ind>0
                if tracest(ind).ishot
                    tracest(ind).ishot = false;
                else
                    tracest(ind).ishot = true;
                end
                save(save_loc,'tracest','-append')
                upz = false;
                ind = already_found(rxpos,rypos,-1); %clear ind text
                move_callback(fh_img)
            end
        end
        if strcmp(event.Key,'p')
            if ind>0
                if tracest(ind).ispair
                    tracest(ind).ispair = false;
                else
                    tracest(ind).ispair = true;
                end
                save(save_loc,'tracest','-append')
                upz = false;
                ind = already_found(rxpos,rypos,-1); %clear ind text
                move_callback(fh_img)
            end
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
            ind = already_found(rxpos,rypos,-1); %clear ind text
            move_callback(fh_img)
        end
        if strcmp(event.Key,'g') || strcmp(event.Key,'c')
            if strcmp(event.Key,'g')
                goto_trace;
            elseif strcmp(event.Key,'c')
                upz = true;
                zoom_in
                goto_trace(ind+1);
            end
        end
        if strcmp(event.Key,'k')
            mask(mask==bitcmp(0,'uint16')) = 0;
            save(save_loc,'tracest','mask')
            ntrace = length(tracest);
        end
        if strcmp(event.Key,'4') || strcmp(event.Key,'numpad4')
            quad = 0;
            while ~any(quad==1:5)
                dh_goto = dialog(...
                    'Units','Normalized',...
                    'Position',[.4 .4 .2 .2],...
                    'Name','Which quad?');
                uicontrol(...
                    'Parent',dh_goto,...
                    'Style','text',...
                    'Units','Normalized',...
                    'Position', [.3 .6 .4 .3],...
                    'String',{'Which quadrant to zoom in on?','1|3','2|4','5 for full screen'});
                qbox = uicontrol(...
                    'Parent',dh_goto,...
                    'Style','edit',...
                    'Units','Normalized',...
                    'Position',[.3 .3 .4 .2],...
                    'Callback','set(gcf,''Windowstyle'',''normal'',''Visible'',''off'')',...
                    'Selected','on');
                uicontrol(qbox)
                waitfor(dh_goto,'Visible','off')
                quad = str2double(qbox.String);
                close(dh_goto)
            end
            if quad<5
                zoomQuad = false(2);
                zoomQuad(quad)=true;
            elseif quad==5
                zoomQuad = true(2);
            end
            move_callback(fh_img)
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
        if ind>length(tracest) || isempty(tracest(ind).frame)
            d = dialog(...
                'Units','Normalized',...
                'Position',[.4 .4 .2 .2],...
                'Name','Finished',...
                'KeyPressFcn',@key_fun);
            uicontrol('Parent',d,...
                'Units','Normalized',...
                'Style','text',...
                'Position',[.2 .5 .6 .1],...
                'String','Esc to save and close all, or click to continue.');
            uicontrol('Parent',d,...
                'Units','Normalized',...
                'Position',[.4 .2 .2 .1],...
                'String','Continue',...
                'Callback','delete(gcf)');
        else
            [area, int, SNR, xpos, ypos] = deal(zeros(1,ml));
            cfr = tracest(ind).frame(1);
            rxpos = tracest(ind).xpos(1);
            xpos(cfr) = rxpos;
            rypos = tracest(ind).ypos(1);
            ypos(cfr) = rypos;
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
            upz = true;
            move_callback(fh_img);
        end
    end
    function delete_trace(~,~)
        for i = ind:length(tracest)-1
            mask(mask(:)==i) = 0;
            mask(mask(:)==i+1) = i;
        end
        tracest(ind) = [];
        mask(mask==bitcmp(0,'uint16')) = 0;
        ntrace = length(tracest);
        scatter_points(cfr);
        save(save_loc,'tracest','mask')
        ind = already_found(rxpos,rypos,-1); %clear ind text
        close(afh)
    end

    function zoom_in
        persistent lh lzlh
        if rxpos>zrad && rxpos<=ss(2)-zrad && rypos>zrad && rypos<=ss(1)-zrad
            [rxpos, rypos] = cofint(rimg,rxpos,rypos,cfr);
            xpos(cfr) = rxpos;
            ypos(cfr) = rypos;
        end
        ind = already_found(rxpos,rypos,cfr);
        if rxpos>zrad && rxpos<=ss(2)-zrad && rypos>zrad && rypos<=ss(1)-zrad
            zimg = rimg(floor(rypos-zrad):floor(rypos+zrad),...
                floor(rxpos-zrad):floor(rxpos+zrad),cfr);
            image(ah_zoom_img,zimg,'CDataMapping','scaled');
            set(ah_zoom_img,'CLim',[minrc maxrc]);
            
%             [oxpos, oypos] = cofint(oimg,rxpos,rypos,cfr);
%             oimgc = oimg(floor(oypos-zrad):floor(oypos+zrad),...
%                 floor(oxpos-zrad):floor(oxpos+zrad),cfr);
%             image(ah_orig_img,oimgc,'CDataMapping','scaled');
%             set(ah_orig_img,'CLim',[minoc maxoc]);
            lzrad = 5*ozrad;
            if floor(rxpos-lzrad)<1
                lzrad = floor(rxpos-1);
            end
            if floor(rxpos+lzrad)>ss(2)
                lzrad = floor(ss(2)-rxpos);
            end
            if floor(rypos-lzrad)<1
                lzrad = floor(rypos-1);
            end
            if floor(rypos+lzrad)>ss(1)
                lzrad = floor(ss(1)-rypos);
            end
            [oxpos, oypos] = cofint(oimg,rxpos,rypos,cfr);
            oimgc = rimg(floor(rypos-lzrad):floor(rypos+lzrad),...
                floor(rxpos-lzrad):floor(rxpos+lzrad),cfr);
            image(ah_orig_img,oimgc,'CDataMapping','scaled');
            set(ah_orig_img,'CLim',[medrc maxrc]);
            
            set(ah_orig_img,'Nextplot','add')
            if exist('lzlh','var'), delete(lzlh); end
            lzlh = line(ah_orig_img,[lzrad+1-zrad-.5 lzrad+1+zrad+.5 lzrad+1+zrad+.5 lzrad+1-zrad-.5 lzrad+1-zrad-.5],...
                [lzrad+1-zrad-.5 lzrad+1-zrad-.5 lzrad+1+zrad+.5 lzrad+1+zrad+.5 lzrad+1-zrad-.5],...
                'color','r','linewidth',.5);
            set(ah_orig_img,'Nextplot','replace')
            
            set(ah_img,'Nextplot','add')
            if exist('lh','var'), delete(lh); end
            lh = line(ah_img,[rxpos-zrad-1 rxpos+zrad rxpos+zrad rxpos-zrad-1 rxpos-zrad-1],...
                [rypos-zrad-1 rypos-zrad-1 rypos+zrad rypos+zrad rypos-zrad-1],...
                'color','r','linewidth',1);
            set(ah_img,'Nextplot','replace')
            mask_clip = mask(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr);
            if ~any(reshape(mask_clip>0,[],1))
                tmp = edge(rimg(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr),'canny',[.5 .9])>0;
                tmp = imfill(tmp,'holes')>0;
                if ind>0
                    tmp = tmp*ind;
                else
                    tmp = tmp*bitcmp(0,'uint16');
                end
                mask(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr) = tmp;
            end
            overlay_mask
            update_area(cfr,true)
            update_int(cfr,true)
            cf_ball
        end
    end
    function overlay_mask
        persistent msh
        delete(msh)
        set(ah_zoom_img,'NextPlot','add')
        tmpmask = mask(floor(rypos-zrad):floor(rypos+zrad),...
            floor(rxpos-zrad):floor(rxpos+zrad),cfr);
        [iy,ix] = find(tmpmask);
        msh = scatter(ah_zoom_img,ix,iy,50,[1 0 1],'linewidth',1);
        set(ah_zoom_img,'NextPlot','replace','XTick',[],'YTick',[])
    end
    function [cx, cy] = cofint(img,tmpx,tmpy,frame)
        tmpimg = double(img(floor(tmpy-zrad):floor(tmpy+zrad),...
            floor(tmpx-zrad):floor(tmpx+zrad),frame));
        tmpmask = mask(floor(tmpy-zrad):floor(tmpy+zrad),...
            floor(tmpx-zrad):floor(tmpx+zrad),cfr);
        if any(reshape(tmpmask,[],1))
            tmp_mask = tmpmask;
        else
            tmp_mask = edge(tmpimg,'canny',[.5 .9])>0;
            tmp_mask = imfill(tmp_mask,'holes')>0;
        end
        tmpimg = tmpimg.*double(tmp_mask>0);
        if sum(tmpimg(:))==0
            cx = tmpx;
            cy = tmpy;
        else
            cx = floor(tmpx-zrad) + sum(tmpimg*(1:size(tmpimg,2))')/sum(tmpimg(:)) - 0.5;
            cy = floor(tmpy-zrad) + sum((1:size(tmpimg,1))*tmpimg)/sum(tmpimg(:)) - 0.5;
        end
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
            area(frame) = sum(sum(mask(floor(rypos-zrad):floor(rypos+zrad),...
                floor(rxpos-zrad):floor(rxpos+zrad),frame)>0));
            if redo
                tracest(ind).area(tracest(ind).frame==frame) = area(frame);
                save(save_loc,'tracest','-append')
            end
        end
        if disp
            if exist('aph','var'), delete(aph); end
            aph = plot(ah_area_graph,area);
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
                if strcmpi(type,'srrf')
                    srrfint(frame) = tracest(ind).srrfint(tracest(ind).frame==frame);
                end
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
            timg = double(oimg(floor(oypos-zrad):floor(oypos+zrad),...
                floor(oxpos-zrad):floor(oxpos+zrad),frame));
            [int(frame),SNR(frame)] = twoDgaussianFitting_theta(timg);
            if strcmpi(type,'srrf')
                timg = double(rimg(floor(rypos-ozrad):floor(rypos+ozrad),...
                    floor(rxpos-ozrad):floor(rxpos+ozrad),frame));
                timg = interpolate_image(timg);
                timg = sort(timg,'descend');
                srrfint(frame) = sum(double(timg(1:100)));
            end
            if redo
                tracest(ind).int(tracest(ind).frame==frame) = int(frame);
                tracest(ind).SNR(tracest(ind).frame==frame) = SNR(frame);
                if strcmpi(type,'srrf')
                    tracest(ind).srrfint(tracest(ind).frame==frame) = srrfint(frame);
                end
                save(save_loc,'tracest','-append')
            end
        end
        if disp
            if exist('iph','var'), delete(iph); end
            iph = plot(ah_int_graph,int);
        end
    end
    function [integ, SNR] = twoDgaussianFitting_theta(img)
        % c(1) = background
        % c(2) = amplitude
        % c(3) = x center
        % c(4) = y center
        % c(5) = theta about z axis
        % c(6) = sd_x
        % c(7) = sd_y
        integ = 0; SNR = 0;
        if zrad>0
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
            axes(ah_int_fit_graph)
            plot(gfit, [xdata(:), ydata(:)], img(:));
            c = coeffvalues(gfit);
            SNR = c(2)/c(1);
            integ = quad2d(gfit,0.5,size(img,2)+.5,0.5,size(img,1)+.5)-c(1)*size(img,1)*size(img,2);
        end
    end
    function int_img = interpolate_image(img)
        [xhave, yhave] = meshgrid(1:size(img,2),1:size(img,1));
        [xwant, ywant] = meshgrid(1:.1:size(img,2),1:.1:size(img,1));
        int_img = griddata(xhave(:),yhave(:),double(img(:)),xwant,ywant,'v4');
    end
    function scatter_points(frame)
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
        set(ah_img,'Nextplot','add')
        bsh = scatter(ah_img,xb-.5,yb-.5,100,'r','x');
        gsh = scatter(ah_img,xg-.5,yg-.5,150,'g','o');
        line(ah_img,[ss(2)/2 ss(2)/2],[0 ss(1)],'color',[.3 .3 .3],'linewidth',.5)
        line(ah_img,[0 ss(2)],[ss(1)/2 ss(1)/2],'color',[.3 .3 .3],'linewidth',.5)
        set(ah_img,'Nextplot','replace')
    end
    function cf_ball
        persistent abh ibh
        if exist('abh','var'), delete(abh); end
        set(ah_area_graph,'NextPlot','add')
        abh = plot(ah_area_graph,cfr,area(cfr),'ro');
        set(ah_area_graph,'NextPlot','replace')
        if exist('ibh','var'), delete(ibh); end
        set(ah_int_graph,'NextPlot','add')
        ibh = plot(ah_int_graph,cfr,int(cfr),'ro');
        set(ah_int_graph,'NextPlot','replace')
    end
    function frame_line(src,loc,col)
        tmpch = get(src,'Children');
        for i = 1:length(tmpch)
            if strcmp(tmpch(i).Type,'line')
                if all(tmpch(i).Color==col)
                    delete(tmpch(i))
                end
            end
        end
        set(src,'NextPlot','add')
        line(src,[loc loc],[.5 1.5],'linewidth',1,'color',col);
        set(src,'NextPlot','replace')
    end
    function select_mask(src,~)
        cp = src.CurrentPoint;
        cp(2) = 1-cp(2);
        cp = ceil(cp([2,1]).*(2*zrad+1));
        cp = cp-zrad-1;
        cp = cp + floor([rypos rxpos]);
        if strcmp(src.SelectionType,'normal')
            if ind>0
                mask(cp(1),cp(2),cfr) = ind;
            else
                mask(cp(1),cp(2),cfr) = bitcmp(0,'uint16');
            end
        elseif strcmp(src.SelectionType,'alt')
            mask(cp(1),cp(2),cfr) = 0;
        end
        mask(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr) =...
            imfill(mask(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr),'holes');
        if ind > 0
            if any(tracest(ind).frame==cfr)
                tracest(ind).area(tracest(ind).frame==cfr) = ...
                    sum(sum(mask(floor(rypos-zrad):floor(rypos+zrad),floor(rxpos-zrad):floor(rxpos+zrad),cfr)>0));
            end
        end
        zoom_in;
        if ind>0
            save(save_loc,'tracest','-append')
        end
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
            if dist<zrad
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
    function sff(~,~)
        ff = cfr;
        frame_line(ah_scroll,cfr,[0 .5 0])
    end
    function slf(~,~)
        lf = cfr;
        frame_line(ah_scroll,cfr,[.7 0 0])
    end
    function save_trace(~,~)
        persistent uih_saved
        if ind>0
            spt = ind;
        else
            spt = ntrace+1;
        end
        try
            tracest(spt).frame = ff:lf;
            tracest(spt).xpos = xpos(ff:lf);
            tracest(spt).ypos = ypos(ff:lf);
            tracest(spt).int = int(ff:lf);
            if strcmpi(type,'srrf')
                tracest(spt).srrfint = srrfint(ff:lf);
            end
            tracest(spt).SNR = SNR(ff:lf);
            tracest(spt).area = area(ff:lf);
            if ind==0
                tracest(spt).ishot = false;
                tracest(spt).ispair = false;
            end
            [i,j,k] = ind2sub(size(mask),find(mask==bitcmp(0,'uint16')));
            cond = i>=min(ypos(ff:lf))-zrad & i<=max(ypos(ff:lf))+zrad &...
                j>=min(xpos(ff:lf))-zrad & j<=max(xpos(ff:lf))+zrad &...
                k>=ff & k<=lf;
            for ijk = 1:length(cond)
                if cond(ijk)
                    mask(i(ijk),j(ijk),k(ijk)) = spt;
                end
            end
            uih_saved = uicontrol(...
                'Parent',fh_text,...
                'Style','Text',...
                'FontSize',15,...
                'Units','Normalized',...
                'Position',[1/3 0 1/3 .1],...
                'String',['Saved trace ' num2str(spt) '!']);
            mask(mask==bitcmp(0,'uint16')) = 0;
            mh.mask(floor(min(ypos(ff:lf))):ceil(max(ypos(ff:lf))),...
                floor(min(xpos(ff:lf))):ceil(max(xpos(ff:lf))),...
                ff:lf) = ...
                mask(floor(min(ypos(ff:lf))):ceil(max(ypos(ff:lf))),...
                floor(min(xpos(ff:lf))):ceil(max(xpos(ff:lf))),...
                ff:lf);
            save(save_loc,'tracest','-append')
            ntrace = length(tracest);
            pause(.5)
        catch
            uih_saved = uicontrol(...
                'Parent',fh_text,...
                'Style','Text',...
                'FontSize',16,...
                'FontWeight','bold',...
                'Units','Normalized',...
                'Position',[1/3 0 1/3 .1],...
                'String','Save failed!!!');
        end
        delete(uih_saved)
        upz = false;
        move_callback(fh_img)
        ind = 0;
    end
end