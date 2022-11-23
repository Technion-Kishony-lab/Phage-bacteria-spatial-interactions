%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot a full phenotype matrix based on the cross-infection assay,  %%%
%%% print host-switch example plaques                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Printing Figure 2 infection matrix\n');
f = filesep;
linkageMethod = 'complete'; 
figure_location =  'Figs' ;
pltMap = [0.8 0.8 0.8; hex2rgb({'#D81B60','#1E88E5','#FFC107','#004D40'})];

%% load data
infectStructFile1 = ['..' f 'script_data' f 'infectStruct1.mat'];
infectStructFile2 = ['..' f 'script_data' f 'infectStruct2.mat'];
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 
BacGenotypeTable = ['..' f 'Tables' f 'Bac_genotypes_uni.xlsx']; 
PhgGenotypeTable = ['..' f 'Tables' f 'Phg_genotypes_uni.xlsx']; 
%% Read source data
phenotypes = table2array(readtable(phenotypeTable)); 
genInfoBac = readtable(BacGenotypeTable,'PreserveVariableNames',true);
genInfoPhg = readtable(PhgGenotypeTable,'PreserveVariableNames',true);
BacIsolates = readtable(BacIsolatesTable);
PhgIsolates = readtable(PhgIsolatesTable);

if ~exist('infectStruct1','var')
    load(infectStructFile1,'infectStruct1')
end
if ~exist('infectStruct2','var')
    load(infectStructFile2,'infectStruct2')
end
infectStruct = [infectStruct1; infectStruct2];
clear infectStruct1 infectStruct2

%% Calculate matrix pairwise distances and dendrogram
tmp = unique(sort(reshape(phenotypes,1,[])));
log_param = tmp(2);
phenMatLog = log(phenotypes+log_param); 
% clusters
pd_bac = pdist(phenMatLog,'naneucdist');
pd_phg = pdist(phenMatLog','naneucdist');
% plates color map
%% plot figure
fig = figure(2);clf;
set(gcf,'name','Infection matrix','units','centimeters','Position',[1 2 18 21])

Bheight = 4.5;
pxl = 0.12; 
mat_h = size(phenMatLog,2)*pxl;
mat_w = size(phenMatLog,1)*pxl;
den_w = 1.2; 
den_h = 1.2; 
margin_left = 2; 
margin_bot = 0.7; 
space_x = 0.08; 
space_y = 0.08; 
plt_w = 0.2; 
plt_h = 0.2; 
letter_fs = 10;
arrLen = 0.4; 
clsTh = 0.4;
% For b:
marg_x = 1; 
marg_y = 0.05;
dif_x = 0.18; 
dif_y = 0.18;
img_wh = 1; 
% phage dendrogram:
ax_den_phg = axes('units','centimeters','position',[margin_left margin_bot+Bheight den_w mat_h]) ;
lkgPhg = linkage(pd_phg,linkageMethod);
d_phg = dendrogram_pinks(lkgPhg,inf,'Orientation','left','ColorThreshold',max(lkgPhg(:,3))*clsTh) ;
set(d_phg,'LineWidth',1)

clus_phg = cluster(linkage(pd_phg,linkageMethod),'MaxClust',12);
ylim([0.5 size(phenMatLog,2)+0.5])
set(gca,'Visible','off')
set(gca,'YDir','reverse')
axlimphg = get(gca, 'XLim');                                         
aylimphg = get(gca, 'YLim');  
text(max(axlimphg)+20,(max(aylimphg)-min(aylimphg))/2+5,...
    'Phages','FontSize',12,'FontWeight','bold','Rotation',90);

% bac dendrogram:
ax_den_bac = axes('units','centimeters','position',[den_w+margin_left+space_x+plt_w+space_x...
    margin_bot+Bheight+mat_h+plt_h+space_y*1.5 mat_w den_h]) ;
lkgBac = linkage(pd_bac,linkageMethod);
d_bac = dendrogram_blues(lkgBac,inf,'Orientation','top','ColorThreshold', max(lkgBac(:,3))*clsTh) ;

set(d_bac,'LineWidth',1)

clus_bac = cluster(linkage(pd_bac,linkageMethod),'Cutoff',11,'Criterion','distance');
xlim([0.5 size(phenMatLog,1)+0.5])
set(gca,'Visible','off')
axlimbac = get(gca, 'XLim');                                         
aylimbac = get(gca, 'YLim');  
text((max(axlimbac)-min(axlimbac))/2-15,max(aylimbac),...
    'Bacteria','FontSize',12,'FontWeight','bold');
text(-16.5,42,'a','color','black','FontSize',letter_fs)

% Phg plate origin:
phgPltNum = PhgIsolates.replicateNum;
exp2Phg = strcmp(PhgIsolates.experiment,'Cont');
ax_plt_phg = axes('units','centimeters','position',[den_w+margin_left+space_x margin_bot+Bheight plt_w mat_h]);
imagesc(phgPltNum,'parent',ax_plt_phg);
hold on
plot(ax_plt_phg,ones(numel(find(exp2Phg)),1),find(exp2Phg),'*','Color','w','MarkerSize',2)
yyaxis right
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                       
ylp = get(ylh, 'Position');
h_plt_phg = gca;
set(h_plt_phg,'XTick',[],'YTick',[],'xcolor','none','ycolor','none')
h_plt_phg.YAxis(2).Label.Color=[0 0 0]; h_plt_phg.YAxis(2).Label.Visible='on';
colormap(ax_plt_phg,pltMap);

yyaxis left
set(gca,'XTick',[],'YTick',[],'xcolor','none','ycolor','none')

% Bac plate origin:
bacPltNum = BacIsolates.replicateNum;
exp2Bac = strcmp(BacIsolates.experiment,'Cont');
ax_plt_bac = axes('units','centimeters','position',[den_w+margin_left+space_x+plt_w+space_x margin_bot+Bheight+mat_h+space_x mat_w plt_h ]);
imagesc(bacPltNum','parent',ax_plt_bac);
hold on
plot(ax_plt_bac,find(exp2Bac),ones(numel(find(exp2Bac)),1),'*','Color','w','MarkerSize',2)
h_plt_bac = gca;
set(h_plt_bac,'XTick',[],'YTick',[],'xcolor','none','ycolor','none')
h_plt_bac.XAxis.Label.Color=[0 0 0]; h_plt_bac.XAxis.Label.Visible='on';
colormap(ax_plt_bac,pltMap);

% Infection matrix:
ax_mat = axes('units','centimeters','position',[den_w+margin_left+space_x+plt_w+space_x margin_bot+Bheight mat_w mat_h]) ;
imagesc(phenMatLog','Parent',ax_mat);
xticks([]);
yticks([]);
axis([0.5, size(phenMatLog,1)+0.5, 0.5, size(phenMatLog,2)+0.5])        
colormap(ax_mat,flipud(gray))
hold on
[nanR, nanC] = find(isnan(phenMatLog'));
plot(ax_mat,[nanC'-0.5;nanC'+0.5],[nanR'-0.5;nanR'+0.5],'color',[0.2 0.2 0.2]);
plot(ax_mat,[nanC'-0.5;nanC'+0.5],[nanR'+0.5;nanR'-0.5],'color',[0.2 0.2 0.2]);
caxis([min(min(phenMatLog)),max(max(phenMatLog))]);
% link axes:
linkaxes([ax_den_phg,ax_mat],'y');
linkaxes([ax_den_bac,ax_mat],'x');

% source data:
phenMatLogTrans = phenMatLog';
save(['output_files' f 'Fig2a.mat'],'phenMatLogTrans');

% Add colorbars
barH = mat_h/30;
barW = mat_w/3;
axes_clrBarInf = axes('units','centimeters','position',...
    [margin_left+space_x+den_w margin_bot+barH*5.5 barW barH]);
infBar = linspace(min(min(phenMatLog)),max(max(phenMatLog)),100);
imagesc(infBar,'Parent',axes_clrBarInf)
colormap(axes_clrBarInf,flipud(gray))
set(gca,'XTick',1:24:100,'XTickLabel',round([infBar(1) infBar(25) infBar(49) infBar(73) infBar(97)],1),...
    'YTick',[],'fontsize',8)
xlabel('Log(Infectivity)')

axes_clrBarPlt = axes('units','centimeters','position',[margin_left+space_x+den_w margin_bot+2*barH barW barH]);
imagesc(1:5,'Parent',axes_clrBarPlt)
hold on
text([0.7,1.9:4.9],repmat(0.9,1,5),{'WT','1','2','3','4'})
set(gca,'XTick',[],'YTick',[],'fontsize',8);
xlabel('Independent Coevolution Replicate')
colormap(axes_clrBarPlt,pltMap)


%% Infection examples for infectivity scale bar
min_plq = 30; % ignore small detections 
bacExamples = {'Cont-20-B','Cont-14-B','DL-27-A'};
phgExamples = {'1_UL-S04-A','1_DL-38-A','2_DR-S16-A'};

for exm = numel(bacExamples):-1:1
    b(exm) = BacIsolates.infectionIdx(contains(erase(BacIsolates.SequenceNms,'Sample_Bac_'),bacExamples{exm}));
    p(exm) = PhgIsolates.infectionIdx(contains(erase(PhgIsolates.SequenceNms,'Sample_Phage'),phgExamples{exm}));
    if p(exm)>96
        plate = 2;
        p(exm) = p(exm)-96;
    else
        plate = 1;
    end
    infVal = log(infectStruct(plate).plaques(b(exm)).infectivity(p(exm)) + log_param) ;
    x = barW * find(infVal<infBar,1) / 100;
    axes_InfExamples =  axes('units','centimeters','position',[margin_left+space_x+den_w+x-(barH*1.5/2)...
        margin_bot+barH*6.5 barH*1.5 barH*1.5 ]);
    bwsize =  infectStruct(plate).params.wellSize*2 - 20;
    drops_centers = infectStruct(plate).params.drops_centers;
    curr_plate = infectStruct(plate).plaques(b(exm)).normalized_plate;
    well_img = curr_plate(drops_centers(p(exm),2)-bwsize:drops_centers(p(exm),2)+bwsize,...
        drops_centers(p(exm),1)-bwsize:drops_centers(p(exm),1)+bwsize);
    well_img(1,:) = 255;
    well_img(size(well_img,1),:) = 255;
    well_img(:,1) = 255;
    well_img(:,size(well_img,1)) = 255;
    imagesc(axes_InfExamples,well_img);
    set(gca,'XTick',[],'YTick',[]);
    colormap(axes_InfExamples,'gray');
    caxis([0 255])
    hold on
    init_thresh = infectStruct(plate).plaques(b(exm)).well_threshold(p(exm)); % in case it was updated recently
    BW = imbinarize(well_img,init_thresh);
    [B,~] = bwboundaries(BW);
    for k = 1:length(B)
        boundary = B{k};
        boundary(boundary(:,1)==1,:) = []; % delete frame line
        boundary(boundary(:,2)==1,:) = [];
        boundary(boundary(:,1)==size(BW,1),:) = [];
        boundary(boundary(:,2)==size(BW,1),:) = [];
        if length(boundary)>min_plq
           plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 0.5);
        end
    end     
    plot([bwsize bwsize],[bwsize*2 bwsize*2-40],'k-','lineWidth',1)
    rectangle('position',[1 1 size(well_img,1)-1 size(well_img,2)-1],'edgecolor',[0 0 0])
end

%% B: plot host switch examples over same/diff replicate and WT bacteria
PhgPlqSmpls = {'1_Cont-08-A','1_Cont-07-B','1_T7_WT_1'};
BacPlqSmpls = {'Cont-05-A','UR-S07-A','MG-Y_1'};
% Phages:
[~,A] = ismember(PhgPlqSmpls{1},[strcat('1_',infectStruct(1).FullPhgNms.Var1);...
        strcat('2_',infectStruct(2).FullPhgNms.Var1)]); % host switch phage
[~,B] = ismember(PhgPlqSmpls{2},[strcat('1_',infectStruct(1).FullPhgNms.Var1);...
        strcat('2_',infectStruct(2).FullPhgNms.Var1)]); % host switch phage % host range phage; before was 49;
phgWT =  find(strcmp(PhgPlqSmpls{3},[strcat('1_',infectStruct(1).FullPhgNms.Var1);...
        strcat('2_',infectStruct(2).FullPhgNms.Var1)])); % 
phgInfAll = [phgWT(1),B,A];
% Bacteria:
samePlateInf = find(strcmp(BacPlqSmpls{1},table2cell(infectStruct(1).FullBacNms)));
diffPlateInf = find(strcmp(BacPlqSmpls{2},table2cell(infectStruct(1).FullBacNms)));
bacWTInf = find(strcmp(BacPlqSmpls{3},table2cell(infectStruct(1).FullBacNms)));
bacInfPlates = [bacWTInf(3),samePlateInf,diffPlateInf];
% Figure Params:
reduce_factor = 0.9;
plate = 1;
dc = infectStruct(plate).params.drops_centers;
bwsize = infectStruct(plate).params.wellSize*2-20;

ylabels = {{'WT'; 'Phage'},{'Host'; 'Range'; '(HR)'},{'Host'; 'Switch'; '(HS)'}};
titles = {{'WT'; 'Bac.'},{'Same'; 'Rep.'},{'Different'; 'Rep.'}};
X0 = margin_left+space_x+mat_w*3/4;
Y0 = 0;
pltRectWidth = bwsize*0.3;
% Location in the matrix figure:
HSline(1) = size(phenMatLog,2) - find(cellfun(@(x) contains(PhgPlqSmpls{3},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')),1);
HSline(2) =  size(phenMatLog,2) - find(cellfun(@(x) contains(PhgPlqSmpls{2},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')));
HSline(3) = size(phenMatLog,2) - find(cellfun(@(x) contains(PhgPlqSmpls{1},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')));
HSline(4) = find(cellfun(@(x) contains(BacPlqSmpls{3},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')),1,'last');
HSline(5) = find(cellfun(@(x) contains(BacPlqSmpls{1},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')));
HSline(6) = find(cellfun(@(x) contains(BacPlqSmpls{2},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')));

for bac = 1:numel(bacInfPlates) 
    curr_plate = infectStruct(plate).plaques(bacInfPlates(bac)).normalized_plate;
    for phg = 1:numel(phgInfAll)
        ax = axes('Units','centimeters','Position',...
            [X0+marg_x+(img_wh+dif_x)*(bac-1) Y0+marg_y+(numel(phgInfAll)-phg)*(dif_y+img_wh) img_wh img_wh ]);
        if bac == 1 && phg ==2 % Center only this image
            well_img = curr_plate(dc(phgInfAll(phg),2)-bwsize:dc(phgInfAll(phg),2)+bwsize,...
                dc(phgInfAll(phg),1)-bwsize-40:dc(phgInfAll(phg),1)+bwsize-40);
        else
            well_img = curr_plate(dc(phgInfAll(phg),2)-bwsize:dc(phgInfAll(phg),2)+bwsize,...
            dc(phgInfAll(phg),1)-bwsize:dc(phgInfAll(phg),1)+bwsize);
        end
        imagesc(ax,well_img);
        colormap(ax,'gray');
        set(gca,'XTick',[],'YTick',[]);
        caxis([0 255])
        hold on
        if bac ==1
            text(-320,140,ylabels{phg},'VerticalAlignment','middle','FontSize',8);
            rectangle('Position',[1,1,pltRectWidth,bwsize*2-2],'FaceColor',...
                pltMap(phgPltNum(size(phenMatLog,2)-HSline(phg))+1,:),'EdgeColor','none');
        end
        if phg ==1
            title(titles{bac})
            rectangle('Position',[1,1,bwsize*2-2,pltRectWidth],'FaceColor',pltMap(bacPltNum(HSline(bac+3))+1,:),'EdgeColor','none');
        end
        if bac==1 && phg==1
            text(-330,-200,'b','color','black','FontSize',letter_fs)

        end
        clear ax
    end
end


figPos = get(gcf,'position');
arrowString = {'WT','HR','HS','WT',{'Same Rep.'},{'Diff. Rep.'}};
xArr = [];
yArr = [];
for i = 3:-1:1
    xArr(i,:) = [den_w+margin_left+space_x+plt_w+space_x+mat_w+arrLen den_w+...
        margin_left+space_x+plt_w+space_x+mat_w ]/figPos(3);
    yArr(i,:) = [margin_bot+Bheight+(HSline(i)+0.5)*pxl ...
        margin_bot+Bheight+(HSline(i)+0.5)*pxl]/figPos(4);
    annotation('textarrow', xArr(i,:),yArr(i,:),'string',arrowString{i},...
        'HorizontalAlignment','center','fontsize',8);
end
for i = 6:-1:4
    xArr(i,:) = [den_w+margin_left+space_x+plt_w+space_x+(HSline(i)-0.5)*pxl  den_w+margin_left+space_x+plt_w+space_x+(HSline(i)-0.5)*pxl ]/figPos(3);
    yArr(i,:) = [margin_bot+Bheight-arrLen margin_bot+Bheight]/figPos(4);
    annotation('textarrow', xArr(i,:),yArr(i,:),'string',arrowString{i},...
        'HorizontalAlignment','center','fontsize',8);
end

print([figure_location f 'Figure2'],'-dpng','-r300');
print([figure_location f 'Figure2'],'-depsc2')