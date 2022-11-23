
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot experimental setup, sample plate images, organism dominance      %%%  
%%% and intensity peak analysis of each plate (initial + cont experiments)%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Printing Figure 1 and Supplementary Figs 1,3\n');
tic
f = filesep;
%% file locations:
masks_main = ['..' f 'input_files' f 'coevolution' f 'masks_main.mat'];
masks_cont = ['..' f 'input_files' f 'coevolution' f 'masks_cont.mat'];
masks_noPhg = ['..' f 'input_files' f 'coevolution' f 'masks_SP.mat'];

sngPltMFile = ['..' f 'script_data' f 'sngPltM.mat'];
sngPltCFile = ['..' f 'script_data' f 'sngPltC.mat'];
sngPltSPFile = ['..' f 'script_data' f 'sngPltSP.mat'];

sngPltMeanMFile = ['..' f 'script_data' f 'sngPltMeanM.mat'];
sngPltMeanCFile = ['..' f 'script_data' f 'sngPltMeanC.mat'];
sngPltMeanSPFile = ['..' f 'script_data' f 'sngPltMeanSP.mat'];

movieParamsFile = ['..' f 'script_data' f 'movieParams.mat'];

ptsMainFile = ['..' f 'input_files' f 'coevolution' f 'intensityCoorMain.xlsx'];
ptsContFile = ['..' f 'input_files' f 'coevolution' f 'intensityCoorCont.xlsx'];
img_setup_file_main = ['..' f 'input_files' f 'exp_setup' f 'setup_Rev2-03.png'];
img_setup_file_cont = ['..' f 'input_files' f 'exp_setup' f 'setup_Rev2-04.png'];
%% file save locations:
figure_location =  'Figs' ;
calculations_location = 'output_files';
%% parameters
load(movieParamsFile,'N','A','ti','xyWin','tWin','plateSize','plateSizeSP');
% Chosen spots for intensity plots over time:
for pt = 4:-1:1
    ptsMain(pt).pts = zeros(3,2);
    tmpTbl = readtable(ptsMainFile,'sheet',pt);
    ptsMain(pt).pts = tmpTbl.Variables;
end
for pt = 4:-1:1
    ptsCont(pt).pts = zeros(3,2);
    tmpTbl = readtable(ptsContFile,'sheet',pt);
    ptsCont(pt).pts = tmpTbl.Variables;
end
%% Main make circle masks per mask 
maskM = load(masks_main,'mask');
% Main Average intensity on small windows: 
meanM = load(sngPltMeanMFile,'sngPltMean');
pltM = load(sngPltMFile,'sngPlt');
%% calculate median/mean peaks number
numPksMain = [];
pkCalcMain = [];
for pl = 4:-1:1
    % window the mask:
    maskWinXYcell = mat2cell(uint8(maskM.mask(pl).mask(1:plateSize,1:plateSize)),...
        xyWin*ones(plateSize/xyWin,1),xyWin*ones(plateSize/xyWin,1));
    maskWinXY = logical(cellfun(@(x) mean2(x),maskWinXYcell));
    peaksPlateMain = [meanM.sngPltMean(pl).pks.numPks];
    numPksMain(pl).np = peaksPlateMain(maskWinXY);
    peakProm = zeros(size(maskWinXY));
    % average peak prominence, -1 where mask is zero 
    for i = 1:size(maskWinXY,2)
        for j = 1:size(maskWinXY,1)
            if maskWinXY(i,j)
                peakProm(i,j) = mean(meanM.sngPltMean(pl).pks(1,sub2ind(size(maskWinXY),i,j)).prominence);
            else
                peakProm(i,j) = -1;
            end
        end
    end
    numPksMain(pl).prominence = peakProm(peakProm~=-1 &~isnan(peakProm));
end
pkCalcMain.allPks = [numPksMain.np];
pkCalcMain.PkMean = mean([numPksMain.np]);
pkCalcMain.PkStd = std([numPksMain.np]);
pkCalcMain.PkMedian = median([numPksMain.np]);
pkCalcMain.PkMeanProminence = vertcat(numPksMain.prominence);
pkCalcMain.PkMeanProminenceMean = mean(pkCalcMain.PkMeanProminence);
save([calculations_location f 'PeaksCalcMain'],'pkCalcMain');
%%  Dynamics of peaks and Max vs Final - Main
% Position figure:
y_cm = 18; 
space_x = 0.45;
margin_x = 0.5; 
margin_y = 1.2;
space_y = 0.45; 
letter_fs = 8;
fs_legend2D = 7;
setup_w = 7.5;
setup_h = 3.5;
frame_size = 2.5;
legend2D_x = frame_size*2/3;
legend2D_y = frame_size*2/3;
img_size = 2.25; 
x_cm = 18; 
bar_width = 0.3;

pl = 2;
margin_x_TL = 0.3;
margin_y_TL = 1;
space_x_imgs = 0.2; 
space_y_imgs = 0.2; 
intensify = 2;
reduce = 2;
timeClr = [252, 247, 135]/255;

peakMarSz = 2.5;

figure(1); clf;
set(gcf,'Name','Dynamics of peaks and Max vs Final - Main','Units',...
    'centimeters','Position',[1 1 x_cm y_cm])

% Load setup figure
img_setup = imread(img_setup_file_main);
ax_setup = axes('units','centimeters','position',[margin_x  margin_y+legend2D_y+3*(img_size+space_y)+space_y/2 setup_w setup_h]);
imshow(img_setup,'Parent',ax_setup)
text(ax_setup.XLim(1) - 30000/(ax_setup.XLim(2) - ax_setup.XLim(1)),0,...
    'a','Parent',ax_setup,'color','black','FontSize',letter_fs)

% Plot time lapse images of one plate (UR) for figure 1
timePointsMain = [8 36 63 105 137 165 207 240 298];
imgN(1) = size(pltM.sngPlt(1).imgs,3); 
tf(1) = imgN(1);
t_main = round(logspace(log10(ti(1)+A(1)),log10(tf(1)+A(1)),N(1))-A(1));
for t = 1:numel(timePointsMain)
    % plot image:
    fprintf('%d out of %d\n',t,numel(timePointsMain))
    if t<=3
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-1) ...
            margin_y+space_y+legend2D_y+(img_size+space_y_imgs)*2 img_size img_size]);
    elseif t<=6
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-3-1)...
            margin_y+space_y+legend2D_y+img_size+space_y_imgs img_size img_size]);
    else
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-6-1)...
            margin_y+space_y+legend2D_y img_size img_size]);
    end
    imshow((pltM.sngPlt(pl).imgs(1:reduce:end,1:reduce:end,t_main(timePointsMain(t)))-pltM.sngPlt(pl).bg(1:reduce:end,1:reduce:end))*intensify,'Parent',ax_TL);
    if t <= 3
        timeLabel = sprintf('%dhr', round(t_main(timePointsMain(t))*10/60));
    else
        timeLabel = sprintf('%0.1fd', t_main(timePointsMain(t))*10/60/24);
    end
    ysize = size(pltM.sngPlt(pl).imgs(1:reduce:end,1:reduce:end,t_main(timePointsMain(t))));
    text(7,ysize(1)-30,timeLabel,'Color',timeClr,'FontSize',letter_fs)
    if t==1
         text(ax_TL.XLim(1) - 30000/(ax_TL.XLim(2) - ax_TL.XLim(1)),0,'b',...
             'Parent',ax_TL,'color','black','FontSize',letter_fs)
    end
end

b_w = (space_x_imgs+img_size)*3;

% Plot overlay of maximal-final (red) and final (green) data:
for pl = 1:4
    maxImg = max(meanM.sngPltMean(pl).imgs,[],3);
    finalImg = meanM.sngPltMean(pl).imgs(:,:,end);
    overlay = cat(3,maxImg-finalImg,finalImg,maxImg-finalImg);
    overlay_adj = imadjust(overlay/255,[min(min(maxImg-finalImg)) min(min(finalImg)) min(min(maxImg-finalImg));...
        max(max(maxImg-finalImg)) max(max(finalImg)) max(max(maxImg-finalImg))]/255,[]);
    ax_maxFinal = axes('Units','centimeters','Position',[margin_x+b_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]); 
    imshow(overlay_adj)
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    if pl==1
        text(ax_maxFinal.XLim(1) - 1000/(ax_maxFinal.XLim(2) - ax_maxFinal.XLim(1)),0,...
            'c','Parent',ax_maxFinal,'color','black','FontSize',letter_fs)
    end
end
c_w = space_x+frame_size;

% Plot number of peak heatmap and intensity curves:
lnWdt = 1;
maxPks = 12;
minPks = 0;
ym = 210;
ptColor = [0.9 0 0.9; 
    0 0.9 0;
    0 0 0.9];
for pl = 1:4
    fprintf('plate %d\n',pl)
    ax_pks = axes('Units','centimeters','Position',[margin_x+b_w+c_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    imagesc(reshape([meanM.sngPltMean(pl).pks.numPks],[plateSize/xyWin,plateSize/xyWin]),'parent',ax_pks);
    pkPltSource = reshape([meanM.sngPltMean(pl).pks.numPks],[plateSize/xyWin,plateSize/xyWin]);
    save([calculations_location f 'Fig1d_' num2str(pl) '.mat'], 'pkPltSource')
    hold on
    ax_pks.Colormap = copper(maxPks);
    pbaspect([1,1,1])
    set(ax_pks,'Xtick',[],'ytick',[]);
    caxis(ax_pks,[0 maxPks]);
    if pl==1
        text(ax_pks.XLim(1) - 1000/(ax_pks.XLim(2) - ax_pks.XLim(1)),0,...
            'd','Parent',ax_pks,'color','black','FontSize',letter_fs)
    end
    
    ax_int = axes('Units','centimeters','Position',[margin_x+b_w+c_w+c_w+space_x*3,...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    
    for pt = 1:size(ptsMain(pl).pts,1)
        plot(ax_int,meanM.sngPltMean(pl).timeMean(:,sub2ind([plateSize/xyWin,plateSize/xyWin],...
            ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))),...
            'Color',ptColor(pt,:),'LineWidth',lnWdt,'LineStyle','-');
        intSource = meanM.sngPltMean(pl).timeMean(:,sub2ind([plateSize/xyWin,plateSize/xyWin],...
            ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1)))';
        hold on
        plot(ax_int,meanM.sngPltMean(pl).pks(sub2ind([plateSize/xyWin,plateSize/xyWin], ...
            ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))).locs,...
            meanM.sngPltMean(pl).pks(sub2ind([plateSize/xyWin,plateSize/xyWin], ...
            ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))).pks,'o',...
            'MarkerFaceColor',ptColor(pt,:),'MarkerEdgeColor',ptColor(pt,:),'MarkerSize',peakMarSz);
        pkSource = meanM.sngPltMean(pl).pks(sub2ind([plateSize/xyWin,plateSize/xyWin], ...
            ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))).locs;
        save([calculations_location f 'Fig1e_' num2str(pl) '_' num2str(pt) '.mat'], 'intSource','pkSource')
        timeLabel = 0:(6*24):size(meanM.sngPltMean(pl).timeMean,1);
        intLabel = 0:50:255;
        ylim([0 ym]);
        plot(ax_pks,ptsMain(pl).pts(pt,1),ptsMain(pl).pts(pt,2),'s','Color',ptColor(pt,:),'LineWidth',1.5);
    end
    ylabel('Intensity,\eta')
    set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',[]);
    if pl==4
        set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',...
            round(timeLabel*10/60/24,1));      
        xlabel('Time, [days]')
    end
    if pl==1
        text(ax_int.XLim(1) - 350000/(ax_int.XLim(2) - ax_int.XLim(1)),ax_int.YLim(2),'e',...
            'Parent',ax_int,'color','black','FontSize',letter_fs)
    end
end

% 2D legend for intensity
ax_2Dclrbr = axes('Units','centimeters','Position',...
    [margin_x+b_w+space_x*2, margin_y, legend2D_x,legend2D_y]);
r = repmat(0:255,256,1);
g = repmat(0:255,256,1)';
imagesc(flipud(cat(3,r,g,r)/255),'parent',ax_2Dclrbr);
set(gca,'XTick',50:100:250,'XTickLabel',50:100:250,'YTick',fliplr(255-(50:100:250)),...
    'YTickLabel',250:-100:50,'fontSize',letter_fs)
pbaspect([1 1 1])
ylabel('\eta_{final}')
xl2D = xlabel('\eta_{max} - \eta_{final}');
ax2D = gca;
ax2D.XRuler.TickLabelGapOffset = 0;
xl2D.Position(2) = xl2D.Position(2) - xl2D.Position(2) * 0.1;
text(100,50,'\bfBacteria\newlineDominate','HorizontalAlignment','center','FontSize',fs_legend2D,'color','k');
text(160,200,'\bfPhage\newlineDominate','HorizontalAlignment','center','FontSize',fs_legend2D,'color','k');

% Colorbar for peaks
ax_clrbr = axes('Units','centimeters','Position',[margin_x+b_w+c_w+space_x,...
    margin_y+3*space_y, frame_size,bar_width]);
imagesc(flipud((minPks:maxPks)),'Parent',ax_clrbr);
ax_clrbr.Colormap = copper(maxPks);
set(ax_clrbr,'ytick',[],'xtick',(minPks:3:maxPks)+1,'XTickLabel',minPks:3:maxPks,'YAxisLocation','left','YDir','reverse'  );
xlPeaks = xlabel('Number of Peaks');
xlPeaks.Position(2) = xlPeaks.Position(2) - xlPeaks.Position(2) * 0.1;

print([figure_location f 'Figure1'],'-dpng','-r300')
print([figure_location f 'Figure1'],'-depsc2')
%% CONT make circle masks per mask 
maskC = load(masks_cont,'mask');
%% CONT Average intensity on small windows: 
meanC = load(sngPltMeanCFile,'sngPltMean');
pltC = load(sngPltCFile,'sngPlt');
%% calculate median/mean peaks number
numPksCont = [];
pkCalcCont = [];
for pl = 4:-1:1
    maskWinXYcell = mat2cell(uint8(maskC.mask(pl).mask(1:plateSize,1:plateSize)),...
        xyWin*ones(plateSize/xyWin,1),xyWin*ones(plateSize/xyWin,1));
    maskWinXY = logical(cellfun(@(x) mean2(x),maskWinXYcell));
    peaksPlateCont = [meanC.sngPltMean(pl).pks.numPks];
    numPksCont(pl).np = peaksPlateCont(maskWinXY);
    peakProm = zeros(size(maskWinXY));
    % average peak prominence, -1 where mask is zero 
    for i = 1:size(maskWinXY,2)
        for j = 1:size(maskWinXY,1)
            if maskWinXY(i,j)
                peakProm(i,j) = mean(meanC.sngPltMean(pl).pks(1,sub2ind(size(maskWinXY),i,j)).prominence);
            else
                peakProm(i,j) = -1;
            end
        end
    end
    numPksCont(pl).prominence = peakProm(peakProm~=-1 &~isnan(peakProm));
end
pkCalcCont.allPks = [numPksCont.np];
pkCalcCont.PkMean = mean([numPksCont.np]);
pkCalcCont.PkStd = std([numPksCont.np]);
pkCalcCont.PkMedian = median([numPksCont.np]);
pkCalcCont.PkMeanProminence = vertcat(numPksCont.prominence);
pkCalcCont.PkMeanProminenceMean = mean(pkCalcCont.PkMeanProminence);
save([calculations_location f 'PeaksCalcCont'],'pkCalcCont');
%%  Dynamics of peaks and Max vs Final - Cont
figure(102); clf;
set(gcf,'Name','Dynamics of peaks and Max vs Final - Cont','Units','centimeters','Position',[1 1 x_cm y_cm])

% Load setup figure
img_setup = imread(img_setup_file_cont);
ax_setup = axes('units','centimeters','position',[margin_x  margin_y+legend2D_y+3*(img_size+space_y)+space_y/2 setup_w setup_h]);
imshow(img_setup,'Parent',ax_setup)
text(ax_setup.XLim(1) - 30000/(ax_setup.XLim(2) - ax_setup.XLim(1)),0,...
    'a','Parent',ax_setup,'color','black','FontSize',letter_fs)

% Plot time lapse images of one plate for figure 1
pl = 1;
intensify = 1.5;
imgN(2) = size(pltC.sngPlt(1).imgs,3); 
tf(2) = imgN(2);
t_cont = round(logspace(log10(ti(2)+A(2)),log10(tf(2)+A(2)),N(2))-A(2));    
timePointsCont =  [47 95 135 165 210 245 275 300 340]; 
for t = 1:numel(timePointsCont)
    % plot image:
    fprintf('%d out of %d\n',t,numel(timePointsCont))
    if t<=3
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-1) ...
            margin_y+space_y+legend2D_y+(img_size+space_y_imgs)*2 img_size img_size]);
    elseif t<=6
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-3-1)...
            margin_y+space_y+legend2D_y+img_size+space_y_imgs img_size img_size]);
    else
        ax_TL = axes('units','centimeters','position',[margin_x+(space_x_imgs+img_size)*(t-6-1)...
            margin_y+space_y+legend2D_y img_size img_size]);
    end
    imshow((pltC.sngPlt(pl).imgs(1:reduce:end,1:reduce:end,t_cont(timePointsCont(t)))-pltC.sngPlt(pl).bg(1:reduce:end,1:reduce:end))*intensify,'Parent',ax_TL);
    if t <= 2
        timeLabel = sprintf('%dhr', round(t_cont(timePointsCont(t))*10/60));
    else
        timeLabel = sprintf('%0.1fd', t_cont(timePointsCont(t))*10/60/24);
    end
    ysize = size(pltC.sngPlt(pl).imgs(1:reduce:end,1:reduce:end,t_cont(timePointsCont(t))));
    text(7,ysize(1)-30,timeLabel,'Color',timeClr,'FontSize',10)
    if t==1
        text(ax_TL.XLim(1) - 30000/(ax_TL.XLim(2) - ax_TL.XLim(1)),0,'b',...
             'Parent',ax_TL,'color','black','FontSize',letter_fs)
    end
end

% Plot overlay of maximal (red) and final (green) data:
for pl = 1:4
    maxImg = max(meanC.sngPltMean(pl).imgs,[],3);
    finalImg = meanC.sngPltMean(pl).imgs(:,:,end);
    overlay = cat(3,maxImg-finalImg,finalImg,maxImg-finalImg);
    overlay_adj = imadjust(overlay/255,[min(min(maxImg-finalImg)) min(min(finalImg)) min(min(maxImg-finalImg));...
        max(max(maxImg-finalImg)) max(max(finalImg)) max(max(maxImg-finalImg))]/255,[]);
    ax_maxFinal = axes('Units','centimeters','Position',[margin_x+b_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]); 
    imshow(overlay_adj)
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    if pl==1
         text(ax_maxFinal.XLim(1) - 1000/(ax_maxFinal.XLim(2) - ax_maxFinal.XLim(1)),0,...
            'c','Parent',ax_maxFinal,'color','black','FontSize',letter_fs)
    end
end

% Plot number of peaks in each time window:
lnWdt = 1;
maxPks = 10;
minPks = 0;
for pl = 1:4
    fprintf('plate %d\n',pl)
    ax_pks = axes('Units','centimeters','Position',[margin_x+b_w+c_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    imagesc(reshape([meanC.sngPltMean(pl).pks.numPks],[plateSize/xyWin,plateSize/xyWin]),'parent',ax_pks);
    hold on
    ax_pks.Colormap = copper(maxPks);
    pbaspect([1,1,1])
    set(ax_pks,'Xtick',[],'ytick',[]);
    caxis(ax_pks,[0 maxPks]);
    if pl==1
        text(ax_pks.XLim(1) - 1000/(ax_pks.XLim(2) - ax_pks.XLim(1)),0,...
            'd','Parent',ax_pks,'color','black','FontSize',letter_fs)
    end
    ax_int = axes('Units','centimeters','Position',[margin_x+b_w+c_w+c_w+space_x*3,...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    for pt = 1:size(ptsCont(pl).pts,1)
        plot(ax_int,meanC.sngPltMean(pl).timeMean(:,sub2ind([plateSize/xyWin,plateSize/xyWin], ptsCont(pl).pts(pt,2), ptsCont(pl).pts(pt,1))),...
            'Color',ptColor(pt,:),'LineWidth',lnWdt,'LineStyle','-');
        hold on
        plot(ax_int,meanC.sngPltMean(pl).pks(sub2ind([plateSize/xyWin,plateSize/xyWin], ptsCont(pl).pts(pt,2), ptsCont(pl).pts(pt,1))).locs,...
            meanC.sngPltMean(pl).pks(sub2ind([plateSize/xyWin,plateSize/xyWin], ptsCont(pl).pts(pt,2), ptsCont(pl).pts(pt,1))).pks,'o',...
            'MarkerFaceColor',ptColor(pt,:),'MarkerEdgeColor',ptColor(pt,:),'MarkerSize',peakMarSz);
        timeLabel = 0:(6*24):size(meanC.sngPltMean(pl).timeMean,1);
        intLabel = 0:50:255;
        ylim([0 210]);
        plot(ax_pks,ptsCont(pl).pts(pt,1),ptsCont(pl).pts(pt,2),'s','Color',ptColor(pt,:),'LineWidth',1.5);
    end
    ylabel('Intensity,\eta')
    set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',[]);
    if pl==4
        set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',round(timeLabel*10/60/24,1));
        xlabel('Time, [days]')
    end
    if pl==1
         text(ax_int.XLim(1) - 450000/(ax_int.XLim(2) - ax_int.XLim(1)),ax_int.YLim(2),'e',...
            'Parent',ax_int,'color','black','FontSize',letter_fs)
    end
end

% 2D legend for intensity
ax_2Dclrbr = axes('Units','centimeters','Position',...
    [margin_x+b_w+space_x*2, margin_y, legend2D_x,legend2D_y]);
r = repmat(0:255,256,1);
g = repmat(0:255,256,1)';
imagesc(flipud(cat(3,r,g,r)/255),'parent',ax_2Dclrbr);
set(gca,'XTick',50:100:250,'XTickLabel',50:100:250,'YTick',fliplr(255-(50:100:250)),...
    'YTickLabel',250:-100:50,'fontSize',10)
pbaspect([1 1 1])
ylabel('\eta_{final}')
xl2D = xlabel('\eta_{max} - \eta_{final}');
ax2D = gca;
ax2D.XRuler.TickLabelGapOffset = 0;
xl2D.Position(2) = xl2D.Position(2) - xl2D.Position(2) * 0.1;
text(100,50,'\bfBacteria\newlineDominate','HorizontalAlignment','center','FontSize',fs_legend2D,'color','k');
text(160,200,'\bfPhage\newlineDominate','HorizontalAlignment','center','FontSize',fs_legend2D,'color','k');
% Colorbar for peaks
ax_clrbr = axes('Units','centimeters','Position',[margin_x+b_w+c_w+space_x,...
    margin_y+3*space_y, frame_size,bar_width]);
imagesc(flipud((minPks:maxPks)),'Parent',ax_clrbr);
ax_clrbr.Colormap = copper(maxPks);
set(ax_clrbr,'ytick',[],'xtick',(minPks:3:maxPks)+1,'XTickLabel',minPks:3:maxPks,'YAxisLocation','left','YDir','reverse'  );
xlabel('Number of Peaks')


print([figure_location f 'SuppFigure1'],'-dpng','-r300')

%% Supplementary Figure 3 - No-Phage control and peak distribution
maskSP = load(masks_noPhg,'mask');
% Main Average intensity on small windows: 
load(sngPltMeanSPFile,'sngSPMean');
load(sngPltSPFile,'sngSP');
%% calculate median/mean peaks number
numPksSP = [];
pkCalcSP = [];
for pl = 4:-1:1
    maskWinXYcell = mat2cell(uint8(maskSP.mask(pl).mask(1:plateSizeSP,1:plateSizeSP)),...
        xyWin*ones(plateSizeSP/xyWin,1),xyWin*ones(plateSizeSP/xyWin,1));
    maskWinXY = logical(cellfun(@(x) mean2(x),maskWinXYcell));
    peaksPlateSP = [sngSPMean(pl).pks.numPks];
    numPksSP(pl).np = peaksPlateSP(maskWinXY);
    peakProm = zeros(size(maskWinXY));
    % average peak prominence, -1 where mask is zero 
    for i = 1:size(maskWinXY,2)
        for j = 1:size(maskWinXY,1)
            if maskWinXY(i,j)
                peakProm(i,j) = mean(sngSPMean(pl).pks(1,sub2ind(size(maskWinXY),i,j)).prominence);
            else
                peakProm(i,j) = -1;
            end
        end
    end
    numPksSP(pl).prominence = peakProm(peakProm~=-1 &~isnan(peakProm));
end
pkCalcSP.allPks = [numPksSP.np];
pkCalcSP.PkMean = mean([numPksSP.np]);
pkCalcSP.PkStd = std([numPksSP.np]);
pkCalcSP.PkMedian = median([numPksSP.np]);
pkCalcSP.PkMeanProminence = vertcat(numPksSP.prominence);
pkCalcSP.PkMeanProminenceMean = mean(pkCalcSP.PkMeanProminence);
save([calculations_location f 'PeaksCalcSP'],'pkCalcSP');

%% Plot Supplementary Figure 3
x_cm_sup3 = 24;
y_cm_sup3 = 16;
margin_x_sup3 = 1.2;
margin_y_sup3 = 1;
histSize = 6;
sp_sup3 = 1.8;
xmin = -1;
xmax = 10;
figure(103); clf;
set(gcf,'Name','Peak Distribution','Units','centimeters','Position',[1 1 x_cm_sup3 y_cm_sup3])

pl = 1;
lineClr = colormap('gray');
ax_int_SP = axes('units','centimeters','position',[margin_x_sup3  ...
    margin_y_sup3+(histSize+sp_sup3) histSize histSize]);
timeLabel = 0:(6*24):size(sngSPMean(pl).timeMean,1);
for i = 1:200:size(sngSPMean(pl).timeMean,2)
    plot(ax_int_SP,sngSPMean(pl).timeMean(:,i),'color',lineClr(ceil(i/200)*4,:))
    hold on
    plot(sngSPMean(pl).pks(i).locs,sngSPMean(pl).pks(i).pks,'o',...
            'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3);
end
xlim([0 size(sngSPMean(pl).timeMean,1)])
set(ax_int_SP,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',round(timeLabel*10/60/24,1));      
xlabel('Time, [days]')
ylabel('Intensity,\eta')
annotation('textbox',[margin_x_sup3/x_cm_sup3,(margin_y_sup3+(histSize+sp_sup3)+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','a','color','black','EdgeColor','none','FontSize',letter_fs+3)


ax_pkSP = axes('units','centimeters','position',[margin_x_sup3+(histSize+sp_sup3)  ...
    margin_y_sup3+(histSize+sp_sup3)  histSize histSize]);
histogram([numPksSP.np],'Parent',ax_pkSP)
xlim([xmin xmax])
xline(pkCalcSP.PkMedian,'color','r','linewidth',2)
xlabel('Detected peaks')
ylabel('Count');
title('Bacterial Migration')
annotation('textbox',[(margin_x_sup3+(histSize+sp_sup3))/x_cm_sup3,(margin_y_sup3+(histSize+sp_sup3)+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','b','color','black','EdgeColor','none','FontSize',letter_fs+3)

ax_promSP = axes('units','centimeters','position',[margin_x_sup3+2*(histSize+sp_sup3)  ...
    margin_y_sup3+(histSize+sp_sup3)  histSize histSize]);
histogram(pkCalcSP.PkMeanProminence,'Normalization','probability','Parent',ax_promSP)
xlim([0 50])
xlabel('Mean Peak Prominence')
ylabel('Density');
title('Bacterial Migration')
annotation('textbox',[(margin_x_sup3+2*(histSize+sp_sup3))/x_cm_sup3,(margin_y_sup3+(histSize+sp_sup3)+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','c','color','black','EdgeColor','none','FontSize',letter_fs+3)

ax_pkM = axes('units','centimeters','position',[margin_x_sup3  ...
    margin_y_sup3 histSize histSize]);
histogram([numPksMain.np],'Parent',ax_pkM)
xlim([xmin xmax])
xline(pkCalcMain.PkMedian,'color','k','linewidth',2)
xline(pkCalcSP.PkMedian,'color','r','linewidth',2)
xlabel('Detected peaks')
ylabel('Count');
title('Coevolution: Initial')
annotation('textbox',[(margin_x_sup3)/x_cm_sup3,(margin_y_sup3+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','d','color','black','EdgeColor','none','FontSize',letter_fs+3)

ax_pkC = axes('units','centimeters','position',[margin_x_sup3+(histSize+sp_sup3)  ...
    margin_y_sup3 histSize histSize]);
histogram([numPksCont.np],'Parent',ax_pkC)
xlim([xmin xmax])
xline(pkCalcCont.PkMedian,'color','k','linewidth',2)
xline(pkCalcSP.PkMedian,'color','r','linewidth',2)
xlabel('Detected peaks')
ylabel('Count');
title('Coevolution: Continual')
annotation('textbox',[(margin_x_sup3+(histSize+sp_sup3))/x_cm_sup3,(margin_y_sup3+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','e','color','black','EdgeColor','none','FontSize',letter_fs+3)

ax_prominence = axes('units','centimeters','position',[margin_x_sup3+2*(histSize+sp_sup3)  ...
    margin_y_sup3 histSize histSize]);
hold on
histogram(pkCalcMain.PkMeanProminence,'Normalization','pdf');
histogram(pkCalcCont.PkMeanProminence,'Normalization','pdf');
annotation('textbox',[(margin_x_sup3+2*(histSize+sp_sup3))/x_cm_sup3,(margin_y_sup3+histSize)/y_cm_sup3 0.07 0.07] ,...
    'string','f','color','black','EdgeColor','none','FontSize',letter_fs+3)
set(gca,'Box','on')
legend({'Initial Coevolution','Continual Coevolution'})
xlim([0 120])
xlabel('Mean Peak Prominence')
ylabel('Density');
title('Coevolution')
fprintf('mean prominence Initial:%0.1f\n mean prominence Continual:%0.1f\n mean prominence bacterial migration:%0.1f\n',...
    pkCalcMain.PkMeanProminenceMean,pkCalcCont.PkMeanProminenceMean,pkCalcSP.PkMeanProminenceMean);
print([figure_location f 'SuppFigure3'],'-dpng','-r300')

%% Full figure for no-phage control (without saving)

figure(1003); clf;
set(gcf,'Name','Dynamics of peaks and Max vs Final - No Phage','Units','centimeters','Position',[1 1 x_cm y_cm])
for pl = 1:4
    maxImg = max(sngSPMean(pl).imgs,[],3);
    finalImg = sngSPMean(pl).imgs(:,:,end);
    overlay = cat(3,maxImg-finalImg,finalImg,maxImg-finalImg);
    overlay_adj = imadjust(overlay/255,[min(min(maxImg-finalImg)) min(min(finalImg)) min(min(maxImg-finalImg));...
        max(max(maxImg-finalImg)) max(max(finalImg)) max(max(maxImg-finalImg))]/255,[]);
    ax_maxFinal = axes('Units','centimeters','Position',[margin_x+b_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]); 
    imshow(overlay_adj)
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    set(gca,'xtick',[],'ytick',[]);
    pbaspect([1,1,1])
    if pl==1
         text(ax_maxFinal.XLim(1) - 1000/(ax_maxFinal.XLim(2) - ax_maxFinal.XLim(1)),0,...
            'c','Parent',ax_maxFinal,'color','black','FontSize',letter_fs)
    end
end

% Plot number of peaks in each time window:
leftFig = space_x+frame_size;
lnWdt = 1;
maxPks = 12;
minPks = 0;
ym_noPhg = 255;
ptColor = [0.9 0 0.9; 
    0 0.9 0;
    0 0 0.9];
for pl = 1:4
    fprintf('plate %d\n',pl)
    ax_pks = axes('Units','centimeters','Position',[margin_x+b_w+c_w+space_x, ...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    imagesc(reshape([sngSPMean(pl).pks.numPks],[plateSizeSP/xyWin,plateSizeSP/xyWin]),'parent',ax_pks);
    hold on
    ax_pks.Colormap = copper(maxPks);
    pbaspect([1,1,1])
    set(ax_pks,'Xtick',[],'ytick',[]);
    caxis(ax_pks,[0 maxPks]);
    if pl==1
        text(ax_pks.XLim(1) - 1000/(ax_pks.XLim(2) - ax_pks.XLim(1)),0,...
            'd','Parent',ax_pks,'color','black','FontSize',letter_fs)
    end
    ax_int = axes('Units','centimeters','Position',[margin_x+b_w+c_w+c_w+space_x*3,...
        (space_y + frame_size)*(4-pl)+margin_y+legend2D_y+space_y, frame_size,frame_size]);
    for pt = 1:size(ptsMain(pl).pts,1)
        plot(ax_int,sngSPMean(pl).timeMean(:,sub2ind([plateSizeSP/xyWin,plateSizeSP/xyWin], ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))),...
            'Color',ptColor(pt,:),'LineWidth',lnWdt,'LineStyle','-');
        hold on
        plot(ax_int,sngSPMean(pl).pks(sub2ind([plateSizeSP/xyWin,plateSizeSP/xyWin], ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))).locs,...
            sngSPMean(pl).pks(sub2ind([plateSizeSP/xyWin,plateSizeSP/xyWin], ptsMain(pl).pts(pt,2), ptsMain(pl).pts(pt,1))).pks,'o',...
            'MarkerFaceColor',ptColor(pt,:),'MarkerEdgeColor',ptColor(pt,:),'MarkerSize',peakMarSz);
        timeLabel = 0:(6*24):size(sngSPMean(pl).timeMean,1);
        intLabel = 0:50:255;
        ylim([0 ym_noPhg]);
        plot(ax_pks,ptsMain(pl).pts(pt,1),ptsMain(pl).pts(pt,2),'s','Color',ptColor(pt,:),'LineWidth',1.5);
    end
    ylabel('Intensity,\eta')
    set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',[]);
    if pl==4
        set(ax_int,'ytick',intLabel,'YTickLabel',intLabel,'xtick',timeLabel,'XTickLabel',round(timeLabel*10/60/24,1));      
        xlabel('Time, [days]')
    end
    if pl==1
        text(ax_int.XLim(1) - 300000/(ax_int.XLim(2) - ax_int.XLim(1)),ax_int.YLim(2),'e',...
            'Parent',ax_int,'color','black','FontSize',letter_fs)
    end
end