function [x,y,t1] = slide_threshold_spot_finder(img)
x = []; y = [];
if ischar(img)
    img = imread(img);
end
ss = size(img);
img = imgaussfilt(img,mean(ss)/300);
set(groot,'Units','Pixels');
screen = get(groot,'ScreenSize');
imgy = 1;
imgx = ss(2)/ss(1)*imgy*screen(4)/screen(3);
if imgx>0.9
    imgx = 0.9;
    imgy = ss(1)/ss(2)*imgx*screen(3)/screen(4);
end
minc = min(img(:)); maxc = max(img(:));
disp_er = true;
t1=0;


fh_img = figure(...
    'units','normalized',...
    'WindowStyle','modal',...
    'Position',[0, 1-imgy, imgx, imgy],...
    'Menubar','none',...
    'Toolbar','none',...
    'Resize','off',...
    'NumberTitle','off',...
    'KeyPressFcn',@key_fun);
ah_img = axes('units','normalized',...
    'position',[0 0 1 1]);
imagesc(img)
axis off

fh_thresh = figure(...
    'units','normalized',...
    'OuterPosition',[imgx 0 min(1-imgx,.1) 1],...
    'SelectionType','extend',...
    'MenuBar','none',...
    'ToolBar','none',...
    'Pointer','crosshair',...
    'Resize','off',...
    'Name','Scroll',...
    'NumberTitle','off',...
    'KeyPressFcn',@key_fun);
ah_thresh = axes('units','normalized','Position',[0 0 1 1]);
imagesc(ah_thresh,(minc:maxc)')
axis xy
axis off

while true
    try
        waitfor(fh_thresh,'SelectionType','normal')
        zp = fh_thresh.CurrentPoint;
        t1 = (maxc-minc)*zp(2)+.5+minc;
        draw_line(ah_thresh,t1-minc,[0 0 0]);
        B = bwboundaries((img>t1),'noholes');
        [x,y] = deal(zeros(1,length(B)));
        for iB = 1:length(B)
            x(iB) = mean(B{iB}(:,2));
            y(iB) = mean(B{iB}(:,1));
        end
        axes(ah_img)
        hold on
        if exist('sh','var'), delete(sh); end
        sh = scatter(x,y,100,'k');
        hold off
        set(fh_thresh,'SelectionType','extend')
    catch ME
        close all
        if disp_er
            rethrow(ME)
        end
        return;
    end
end

    function key_fun(~,event)
        if strcmp(event.Key,'escape')
            disp_er = false;
            close all
        end
    end
    function draw_line(src,loc,col)
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
        line([.5 1.5],[loc loc],'linewidth',1,'color',col);
        hold off
    end
end