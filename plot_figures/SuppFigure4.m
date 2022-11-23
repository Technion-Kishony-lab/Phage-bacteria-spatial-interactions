%%% Supplementary Figure 4- show the final coevolution images together with
%%% locations and names of all samples
f = filesep;
out_imgs_main = ['..' f 'script_data' f 'imgs_main']; % folder for movie images - main experiment
out_imgs_cont = ['..' f 'script_data' f 'imgs_cont']; % folder for movie images - main experiment
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx'];
figure_location =  'Figs' ;

%%
BacIso = readtable(BacIsolatesTable);
PhgIso = readtable(PhgIsolatesTable);
%%
rawIni = dir([out_imgs_main f '*.jpg']);
rawCont = dir([out_imgs_cont f '*.jpg']);
imgNmsIni = natsort({rawIni.name}');
imgNmsCont = natsort({rawCont.name}');
%%
smplImgIni = imread([out_imgs_main f  imgNmsIni{end}]);
smplImgCont = imread([out_imgs_cont f  imgNmsCont{end}]);

%% Plot samples on plate image - initial coevolution
iBac = find(strcmp(BacIso.experiment,'Init'));
iPhg = find(strcmp(PhgIso.experiment,'Init'));
fs = 4;
ms = 2;
regClr = 'c';
ini2contClr = 'm';
left1 = 0;
left2 = 0;
up = 25;
dn = 20;
ang = -25;
figure(111); clf;
set(gcf,'name','Samples','units','centimeters','position',[5 1 20 20])
[coorIni,u,r] = unique([BacIso.coordinates_1(iBac),BacIso.coordinates_2(iBac);...
    PhgIso.coordinates_1(iPhg) , PhgIso.coordinates_2(iPhg)],'rows');
iniBacSq = BacIso.SequenceNms(iBac);
ini2contSmpl = {'S01','S05','S09','S13'};

subplot_tight(1,2,1)
imshow(smplImgIni)
hold on
for cr = 1:size(coorIni,1)
    clear bacSmp multiNm phgSmp
    bacSmp = find(BacIso.coordinates_1==coorIni(cr,1) & BacIso.coordinates_2==coorIni(cr,2) &...
        strcmp(BacIso.experiment,'Init'));
    if contains(BacIso.SequenceNms(bacSmp),ini2contSmpl)
        clr = ini2contClr;
    else
        clr = regClr;
    end
    plot(coorIni(cr,1),coorIni(cr,2),'o','MarkerFaceColor',clr,'MarkerSize',ms,'MarkerEdgeColor','none')
    

    if numel(bacSmp)==1
        text(coorIni(cr,1)-left1,coorIni(cr,2)-up,BacIso.bacNms{bacSmp},'Color',clr,'FontSize',fs,'rotation',ang)
    elseif numel(bacSmp)>1
        multiNm = BacIso.bacNms{bacSmp(1)};
        for sp = 2:numel(bacSmp)
            multiNm = [multiNm ',' erase(BacIso.bacNms{bacSmp(sp)},'Bac')];
        end
        text(coorIni(cr,1)-left2,coorIni(cr,2)-up,multiNm,'Color',clr,'FontSize',fs,'rotation',ang)
    end
    phgSmp = find(PhgIso.coordinates_1==coorIni(cr,1) & PhgIso.coordinates_2==coorIni(cr,2) &...
        strcmp(PhgIso.experiment,'Init'));
    if numel(phgSmp)==1
        text(coorIni(cr,1)-left1,coorIni(cr,2)+dn,PhgIso.phgNms{phgSmp},'Color',clr,'FontSize',fs,'rotation',ang)
    elseif numel(phgSmp)>1
        multiNm = PhgIso.phgNms{phgSmp(1)};
        for sp = 2:numel(phgSmp)
            multiNm = [multiNm ',' erase(PhgIso.phgNms{phgSmp(sp)},'Phg')];
        end
        text(coorIni(cr,1)-left2,coorIni(cr,2)+dn,multiNm,'Color',clr,'FontSize',fs,'rotation',ang)
    end
    
end
text(-80,0,'a','fontsize',12)
% Plot samples on plate image - Cont coevolution
cBac = find(strcmp(BacIso.experiment,'Cont'));
cPhg = find(strcmp(PhgIso.experiment,'Cont'));

coorCont = unique([BacIso.coordinates_1(cBac),BacIso.coordinates_2(cBac);...
    PhgIso.coordinates_1(cPhg) , PhgIso.coordinates_2(cPhg)],'rows');
subplot_tight(1,2,2)

imshow(smplImgCont)
hold on
for cr = 1:size(coorCont,1)
    clear bacSmp multiNm phgSmp
    plot(coorCont(cr,1),coorCont(cr,2),'o','MarkerFaceColor',clr,'MarkerSize',ms,'MarkerEdgeColor','none')
    bacSmp = find(BacIso.coordinates_1==coorCont(cr,1) & BacIso.coordinates_2==coorCont(cr,2) &...
        strcmp(BacIso.experiment,'Cont'));
    if numel(bacSmp)==1
        text(coorCont(cr,1)-left1,coorCont(cr,2)-up,BacIso.bacNms{bacSmp},'Color',clr,'FontSize',fs,'rotation',ang)
    elseif numel(bacSmp)>1
        multiNm = BacIso.bacNms{bacSmp(1)};
        for sp = 2:numel(bacSmp)
            multiNm = [multiNm ',' erase(BacIso.bacNms{bacSmp(sp)},'Bac')];
        end
        text(coorCont(cr,1)-left2,coorCont(cr,2)-up,multiNm,'Color',clr,'FontSize',fs,'rotation',ang)
    end
    phgSmp = find(PhgIso.coordinates_1==coorCont(cr,1) & PhgIso.coordinates_2==coorCont(cr,2) &...
        strcmp(PhgIso.experiment,'Cont'));
    if numel(phgSmp)==1
        text(coorCont(cr,1)-left1,coorCont(cr,2)+dn,PhgIso.phgNms{phgSmp},'Color',clr,'FontSize',fs,'rotation',ang)
    elseif numel(phgSmp)>1
        multiNm = PhgIso.phgNms{phgSmp(1)};
        for sp = 2:numel(phgSmp)
            multiNm = [multiNm ',' erase(PhgIso.phgNms{phgSmp(sp)},'Phg')];
        end
        text(coorCont(cr,1)-left2,coorCont(cr,2)+dn,multiNm,'Color',clr,'FontSize',fs,'rotation',ang)
    end
    
end
text(-80,0,'b','fontsize',12)
print([figure_location f 'SuppFigure4'],'-dpng','-r300');