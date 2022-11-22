%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot all infection assay images  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
fprintf('Printing Supplmentary Figure 11\n');

f = filesep;
figure_location = ['Figs' f 'SuppFig11'];
if ~exist(['Figs' filesep 'SuppFig11'],'dir')
     mkdir(['Figs' filesep 'SuppFig11'])
end
% load data
infectStructFile1 = ['..' f 'script_data' f 'infectStruct1.mat'];
infectStructFile2 = ['..' f 'script_data' f 'infectStruct2.mat'];
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 

BacIsolates = readtable(BacIsolatesTable);
PhgIsolates = readtable(PhgIsolatesTable);
load(infectStructFile1,'infectStruct1')
load(infectStructFile2,'infectStruct2')
infectStruct = [infectStruct1; infectStruct2];
clear infectStruct1 infectStruct2
%% Infectivity parameters
% figure parameters
w = 0.8;
h = 0.8; 
dx = 0;
dy = 0.18;
W = 18;
H = 22;
margin_x = 0.2;
margin_y = 0;
fs = 7;
lw = 0.5;
reduce_factor = 0.9;
drops_centers1 = infectStruct(1).params.drops_centers;
drops_centers2 = infectStruct(2).params.drops_centers;

jump = 3;
min_plq = 30/(jump^2); % ignore small detections
wsFull = infectStruct(1).params.wellSize; % same for both experiments
bwsizeFull = wsFull*2-20;
ws = round(wsFull/jump); % same for both experiments
bwsize = round(ws*2-(20/jump));
%%
col = 1;
colF = 18;
row0 = ceil(116/18)*3; 21; %23;
row = row0;
figNum = 200;
fprintf('Supp Figure 11 - Page %d\n',figNum-200+1);
fig = figure(figNum);clf;
set(gcf,'units','centimeters','position',[1 1 W H]);
for b = 1:size(BacIsolates,1)
    fprintf('Bac %d\n',b);
    bacInfecNum = BacIsolates.infectionIdx(b);
    curr_plate1 = infectStruct(1).plaques(bacInfecNum).normalized_plate;
    curr_plate2 = infectStruct(2).plaques(bacInfecNum).normalized_plate;  
    for p = 1:size(PhgIsolates,1)
        if p==1 &&  mod(b,3)==1 && b>1 % new page every 4th bac 
            print([figure_location f num2str(figNum) '.png'],'-dpng','-r300');
            close(fig);
            figNum = figNum+1;
            fprintf('Page %d\n',figNum-200+1);
            fig = figure(figNum);
            set(gcf,'units','centimeters','position',[1 1 W H]);
            row = row0;
            col = 1;
        elseif  p==1 &&  mod(b,3)~=1 % new line every bac
            row = row - 1;
            col = 1;
        elseif p>1 && col>colF % new line in same bac
            row = row - 1;
            col = 1;
        end
        if mod(b,3)==1
            ax = axes('Units','centimeters','Position',[margin_x+(col-1)*(w+dx) ...
                margin_y+(row-1)*(h+dy)+2*dy w h]);
        elseif mod(b,3)==2
            ax = axes('Units','centimeters','Position',[margin_x+(col-1)*(w+dx) ...
                margin_y+(row-1)*(h+dy)+dy w h]);
        else 
            ax = axes('Units','centimeters','Position',[margin_x+(col-1)*(w+dx) ...
                margin_y+(row-1)*(h+dy) w h]);
        end

        phgInfecNum = PhgIsolates.infectionIdx(p);
        if phgInfecNum>96
            plt = 2;
            phgInfecNum = phgInfecNum-96;
            drops_centers = drops_centers2;
            curr_plate = curr_plate2;
        else
            plt = 1;
            drops_centers = drops_centers1;
            curr_plate = curr_plate1;
        end

        well_img_full = curr_plate(drops_centers(phgInfecNum,2)-bwsizeFull:drops_centers(phgInfecNum,2)+bwsizeFull,...    
            drops_centers(phgInfecNum,1)-bwsizeFull:drops_centers(phgInfecNum,1)+bwsizeFull);
        well_img = imresize(well_img_full,1/jump,'box');
        well_img(1,:) = 255;
        well_img(size(well_img,1),:) = 255;
        well_img(:,1) = 255;
        well_img(:,size(well_img,1)) = 255;
        imagesc(double(well_img),'Parent',ax);
        colormap gray 
        set(gca,'XTick',[],'YTick',[])
        caxis([0 255])
        hold on
        thresh = infectStruct(plt).plaques(bacInfecNum).well_threshold(phgInfecNum);

        BW = imbinarize(well_img,thresh);
        [B,~] = bwboundaries(BW);
        for k = 1:length(B)
            boundary = B{k};
            boundary(boundary(:,1)==1,:) = []; % delete frame line
            boundary(boundary(:,2)==1,:) = [];
            boundary(boundary(:,1)==size(BW,1),:) = [];
            boundary(boundary(:,2)==size(BW,1),:) = [];
            if length(boundary)>min_plq
               plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', lw);
            end
        end
        rectangle(ax,'Position',[floor(size(well_img,2)/2)-ws floor(size(well_img,1)/2)-ws 2*ws 2*ws],...
            'EdgeColor','b','LineWidth',0.3, 'LineStyle', '--');
        rectangle(ax,'position',[1 1 size(well_img,1) size(well_img,2)],'edgecolor',[1 1 1])
        text(2,size(well_img,2)/jump-3,num2str(p),'fontSize',fs);
        if p==1
            text(3,-17,BacIsolates.bacNms{b},'fontSize',fs+1);
        end
        col = col+1;
        clear well_img B boundary
    end
    clear curr_plate
end

print([figure_location f num2str(figNum) '.png'],'-dpng','-r300');
