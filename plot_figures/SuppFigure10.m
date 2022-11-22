%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print plaque assay results of opgG and   %%%
%%% mlaA engineered mutants                  %%%
%%% (mutations in mCherry strai              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Printing Supplmentary Figure 10 mlaA opgG dye swap\n');
f = filesep;
anaFolder_opgG = ['..' f 'script_data' f 'opgGmlaA' f 'opgG_analysis'];
imgsOpgGfn = 'imgsMagenta.mat';
imgFolder_opgG = ['..' f 'input_files' f 'opgGmlaA' f '20210129_opgG_fluor' f 'newImaging'];
figFolder = 'Figs';
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
%%
if ~exist(['..' f 'script_data' f 'opgGmlaA'],'dir')
    mkdir(['..' f 'script_data' f 'opgGmlaA'])
end
if ~exist(['..' f 'script_data' f 'opgGmlaA' f 'opgG_analysis'],'dir')
    mkdir(['..' f 'script_data' f 'opgGmlaA' f 'opgG_analysis'])
end
%% Read images:
if ~exist([anaFolder_opgG f imgsOpgGfn],'file')
    imgsy = dir([imgFolder_opgG f '*_YFP.CR2']);
    imgsm = dir([imgFolder_opgG f '*_RFP.CR2']);
    bgY = imread([imgFolder_opgG f 'BG_YFP.CR2']);
    bgY_sm = imgaussfilt(bgY,10);
    bgR = imread([imgFolder_opgG f 'BG_RFP.CR2']);
    bgR_sm = imgaussfilt(bgR,10);
    bg = cat(3,bgR_sm(:,:,1),bgY_sm(:,:,2),bgR_sm(:,:,1));

    for im = numel(imgsy):-1:1
        fprintf('Loading img %d\n',im)
        imgy = imread([imgFolder_opgG f imgsy(im).name]);
        imgm = imread([imgFolder_opgG f imgsm(im).name]);
        imgs(:,:,:,im) = flipud(cat(3,imgm(:,:,1),imgy(:,:,2),imgm(:,:,1))-bg);

        imgNms{im} = imgsy(im).name(1:end-8);
    end
    clear imgy imgm bgY bgY_sm bgR bgR_sm bg imgsy imgsm im
    save([anaFolder_opgG f imgsOpgGfn],'imgs','imgNms')
else
    load([anaFolder_opgG f imgsOpgGfn]);
end
%% Get drop centers from drop cropping:
if ~exist([anaFolder_opgG f 'drops_centers.mat'],'file')
    get_drops_centers([imgFolder_opgG f],[anaFolder_opgG f],1)
else 
    load([anaFolder_opgG f 'drops_centers.mat'])
end

% Chosen dilution from each phage:
if ~exist([anaFolder_opgG f 'dilutionVec.mat'],'file')
    dilutionVec = [5,2,1,1,3,0,1,2,2,2,1,3]; 
    save([anaFolder_opgG f 'dilutionVec.mat'],'dilutionVec')
else
    load([anaFolder_opgG f 'dilutionVec.mat']);
end

%%
allPhgNms = {'Phage1_UR_12-B'
'Phage2_Cont-21-A '
'Phage1_UR-S05-B'
'Phage2_Cont-22-A' 
'Phage1_Cont-16-A'
'LB'
'Phage1_cont-08-A'
'Phage1_DL-38-A'
'Phage1_Cont-07-A ' 
'Phage1_Cont-07-B'
'Phage2_DL-20-A'
'Phage_T7_WT_1'};
phgIsoTbl = readtable(PhgIsolatesTable);
%% Plot both colors in one figure, in centimeters: 
fig = figure(110); clf; 
set(gcf,'Units','centimeters', 'position',[1 2 ...
    8.8 13])
r = 180; % drop 'radius' cropped from plate
maxMinPxl = 100;
letter_fs = 14;
imgNmsVecWT = {'WmWy'};
imgNmsVecM = {'GmWy','JmWy','JGmWy'}; 
imgTitsWT = {'WT','WT'};
imgTitsM = {{'WT','opgG*'},{'WT','waaJ*'},{'WT','waaJ*opgG*'}}; 
margin_x_left = 80;
margin_x_right = 20;
plotPhg = [4,11,12];
[~,phgNums] = ismember(allPhgNms(plotPhg(1:2)),erase(phgIsoTbl.SequenceNms,'Sample_'));
phgTitles = [phgIsoTbl.phgNms(phgNums);{'WT T7'}];
[~,relImgsWT] = ismember(imgNmsVecWT,imgNms);
[~,relImgsM] = ismember(imgNmsVecM,imgNms);
w = 1.6;
h = 1.6; 
gap_x = 0.2; 
gap_y = 1.5;
margin_x = 1.5;
margin_y = 0.5;
gap_mlaAopgG = (h+gap_y)*numel(imgNmsVecM);
titClr = [0.9 0 0.9; 0 0.8 0.2];
letters = [-250, -60];
axesXY = [margin_x margin_y w h];
textX = 15;
textY = 325;
fntSz = 8;
numClr = 'k';

% Plot WTs
for phg = 1:numel(plotPhg)
    alldrops = imgs(drops_centers(1,2)-r:drops_centers(8,2)+r,drops_centers(1,1)-r:drops_centers(96,1)+r,:,relImgsWT);
    maxR = double(median(maxk(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
    maxY = double(median(maxk(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
    minR = double(median(mink(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
    minY = double(median(mink(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
    crpImg = imgs(drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)-r:...
        drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)+r,...
        drops_centers((plotPhg(phg)-1)*8+1,1)-r:drops_centers((plotPhg(phg)-1)*8+1,1)+r,:,relImgsWT);
    crpAdjImg = imadjust(crpImg,[minR,minY,minR; maxR,maxY,maxR]);

    ax = axes('Units','centimeters','position',axesXY+[0 (h+gap_y)*(phg-1) 0 0]);
    imagesc(crpAdjImg,'Parent',ax);
    text(textX,textY,'1','fontSize',fntSz,'Color',numClr);
    ax.XTick = [];
    ax.YTick = [];
    tit1 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(2,:)) '}' imgTitsWT{1}];
    tit2 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(1,:)) '}'  imgTitsWT{2}];
    t = title({tit1,tit2},'fontsize',10, 'Interpreter', 'tex');
    ylh = ylabel(phgTitles{phg});
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle',...
        'HorizontalAlignment','right','fontSize',10)
    switch phg
        case 3
            text(letters(1),letters(2),'b','color','black','FontSize',letter_fs)
        case 2
            text(letters(1),letters(2),'c','color','black','FontSize',letter_fs)
        case 1
            text(letters(1),letters(2),'d','color','black','FontSize',letter_fs)
    end
end
% Plot Ms
for phg = 1:numel(plotPhg)
    for bac = 1:numel(relImgsM)
        alldrops = imgs(drops_centers(1,2)-r:drops_centers(8,2)+r,drops_centers(1,1)-r:drops_centers(96,1)+r,:,relImgsM(bac));
        maxR = double(median(maxk(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
        maxY = double(median(maxk(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
        minR = double(median(mink(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
        minY = double(median(mink(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
        crpImg = imgs(drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)-r:...
            drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)+r,...
            drops_centers((plotPhg(phg)-1)*8+1,1)-r:drops_centers((plotPhg(phg)-1)*8+1,1)+r,:,relImgsM(bac));
        crpAdjImg = imadjust(crpImg,[minR,minY,minR; maxR,maxY,maxR]);
        ax = axes('Units','centimeters','position',axesXY+[(gap_x+w)*(bac) (h+gap_y)*(phg-1) 0 0]);
        imagesc(crpAdjImg,'Parent',ax);
        text(textX,textY,num2str(bac+1),'fontSize',fntSz,'Color',numClr);
        ax.XTick = [];
        ax.YTick = [];
        tit1 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(2,:)) '}' imgTitsM{bac}{1}];
        tit2 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(1,:)) '}'  imgTitsM{bac}{2}];
        t = title({tit1,tit2},'fontsize',10, 'Interpreter', 'tex');
    end 
end

%%
anaFolder_mlaA = ['..' f 'script_data' f 'opgGmlaA' f 'mlaA_analysis'];
imgFolder_mlaA = ['..' f 'input_files' f 'opgGmlaA' f '20210325_MG_mlaA_preinduced' f 'images'];
imgsMlaAfn = 'imgsMagenta.mat';
flip = 0;
%%
if ~exist(['..' f 'script_data' f 'opgGmlaA'],'dir')
    mkdir(['..' f 'script_data' f 'opgGmlaA'])
end
if ~exist(['..' f 'script_data' f 'opgGmlaA' f 'mlaA_analysis'],'dir')
    mkdir(['..' f 'script_data' f 'opgGmlaA' f 'mlaA_analysis'])
end
%% Read images:
if ~exist([anaFolder_mlaA f imgsMlaAfn],'file')
    imgsy = dir([imgFolder_mlaA f '*_YFP.CR2']);
    imgsm = dir([imgFolder_mlaA f '*_RFP.CR2']);
    bgY = imread([imgFolder_mlaA f 'BG_YFP.CR2']);
    bgY_sm = imgaussfilt(bgY,10);
    bgR = imread([imgFolder_mlaA f 'BG_RFP.CR2']);
    bgR_sm = imgaussfilt(bgR,10);
    bg = cat(3,bgR_sm(:,:,1),bgY_sm(:,:,2),bgR_sm(:,:,1));
    for im = numel(imgsy):-1:1
        fprintf('Loading img %d\n',im)
        imgy = imread([imgFolder_mlaA f imgsy(im).name]);
        imgm = imread([imgFolder_mlaA f imgsm(im).name]);
        if flip==1
            imgs(:,:,:,im) = flipud(cat(3,imgm(:,:,1),imgy(:,:,2),imgm(:,:,1))-bg);
        else
            imgs(:,:,:,im) = cat(3,imgm(:,:,1),imgy(:,:,2),imgm(:,:,1))-bg;
        end
        imgNms{im} = imgsy(im).name(1:end-8);
    end
    clear imgy imgm bgY bgY_sm bgR bgR_sm bg imgsy imgsm im
    save([anaFolder_mlaA f imgsMlaAfn],'imgs','imgNms')
else
    load([anaFolder_mlaA f imgsMlaAfn]);
end
%% Get drop centers from drop cropping:
if ~exist([anaFolder_mlaA f 'drops_centers.mat'],'file')
    get_drops_centers([imgFolder_mlaA f],[anaFolder_mlaA f],flip)
else 
    load([anaFolder_mlaA f 'drops_centers.mat'])
end

%% Chosen dilution from each phage:
if ~exist([anaFolder_mlaA f 'dilutionVec.mat'],'file')
    dilutionVec = repmat(7,1,12); 
    save([anaFolder_mlaA f 'dilutionVec.mat'],'dilutionVec')
else
    load([anaFolder_mlaA f 'dilutionVec.mat']);
end
%% Plot both colors in one figure, in centimeters: 
phgTitles = {'WT T7'};
r = 180;
imgNmsVecWT = {'Ay_Am'};
imgNmsVecM = {'Ay_102m'}; 
imgTitsWT = {'p-mlaA','p-mlaA'};
imgTitsM = {'p-mlaA','p-mlaA*'}; 
maxMinPxl = 255;
gamma = 1.3;
plotPhg = 2;
[~,relImgsWT] = ismember(imgNmsVecWT,imgNms);
[~,relImgsM] = ismember(imgNmsVecM,imgNms);

axesXY = [margin_x margin_y+gap_mlaAopgG w h];

for phg = 1:numel(plotPhg)
    alldrops = imgs(drops_centers(1,2)-r:drops_centers(8,2)+r,drops_centers(1,1)-r:drops_centers(96,1)+r,:,relImgsWT);
    maxR = double(median(maxk(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
    maxY = double(median(maxk(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
    minR = double(median(mink(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
    minY = double(median(mink(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
    crpImg = imgs(drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)-r:...
        drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)+r,...
        drops_centers((plotPhg(phg)-1)*8+1,1)-r:drops_centers((plotPhg(phg)-1)*8+1,1)+r,:,relImgsWT);
    crpAdjImg = imadjust(crpImg,[minR,minY,minR; maxR,maxY,maxR],[],gamma);
    ax = axes('Units','centimeters','position',axesXY+[0 (h*7/3+gap_y)*(phg-1) 0 0]);
    imagesc(crpAdjImg,'Parent',ax);
    ax.XTick = [];
    ax.YTick = [];
    tit1 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(2,:)) '}' imgTitsWT{1}];
    tit2 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(1,:)) '}'  imgTitsWT{2}];
    t = title({tit1,tit2},'fontsize',10, 'Interpreter', 'tex');
    ylh = ylabel(phgTitles{phg});
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle',...
        'HorizontalAlignment','right','fontSize',10)
    text(letters(1),letters(2),'a','color','black','FontSize',letter_fs)
end
% Plot Ms
for phg = 1:numel(plotPhg)
    for bac = 1:numel(relImgsM)
        alldrops = imgs(drops_centers(1,2)-r:drops_centers(8,2)+r,drops_centers(1,1)-r:drops_centers(96,1)+r,:,relImgsM(bac));
        maxR = double(median(maxk(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
        maxY = double(median(maxk(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
        minR = double(median(mink(reshape(alldrops(:,:,1),[],1),maxMinPxl)))/255;
        minY = double(median(mink(reshape(alldrops(:,:,2),[],1),maxMinPxl)))/255;
        crpImg = imgs(drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)-r:...
            drops_centers((plotPhg(phg)-1)*8+dilutionVec(plotPhg(phg)),2)+r,...
            drops_centers((plotPhg(phg)-1)*8+1,1)-r:drops_centers((plotPhg(phg)-1)*8+1,1)+r,:,relImgsM(bac));
        crpAdjImg = imadjust(crpImg,[minR,minY,minR; maxR,maxY,maxR],[],gamma);
        ax = axes('Units','centimeters','position',axesXY+[(gap_x+w)*(bac) (h+gap_y)*(phg-1) 0 0]);
        imagesc(crpAdjImg,'Parent',ax);
        ax.XTick = [];
        ax.YTick = [];
        tit1 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(2,:)) '}' imgTitsM{1}];
        tit2 = ['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', titClr(1,:)) '}'  imgTitsM{2}];
        t = title({tit1,tit2},...
        'fontsize',10, 'Interpreter', 'tex');
    end 
end

print([figFolder f 'SuppFigure10'],'-dpng','-r300')


