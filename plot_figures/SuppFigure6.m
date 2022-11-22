%% SuppFigure6: Reproducibility of Infection score (vs. turbidity) 
%% Run after running Figure2.m
f = filesep;
% Read Infectivity and turbidity tables:
turbidTblFile = ['..' f 'Tables' f 'Turbidity.xlsx']; 
turbidTblOrgOrderFile = ['..' f 'Tables' f 'TurbidityOrgOrder.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 
turbidity = table2array(readtable(turbidTblFile)); 
turbidityInvOrgOrder = 1-table2array(readtable(turbidTblOrgOrderFile)); 
infectivity = table2array(readtable(phenotypeTable)); 
linkageMethod = 'complete';
figure_location =  'Figs' ;

%% Remove nans:
[nanR, nanC] = find(isnan(infectivity));
nanInd = sub2ind(size(infectivity),nanR, nanC);
infectivityNoNan = infectivity;
infectivityNoNan(nanInd) = [];
turbidityOrgInvOrderNoNan = turbidityInvOrgOrder;
turbidityOrgInvOrderNoNan(nanInd) = [];

%% Turbidity map (original order)
figure(107);clf;
set(gcf,'name','Infection matrix - turbidity','units','centimeters','Position',[1 2 18 16])
Bheight = 0; 4.5;
pxl = 0.113;
mat_h = size(phenMatLog,2)*pxl;
mat_w = size(phenMatLog,1)*pxl;
margin_left = 0.5; 
margin_bot = 0.7; 
space_x = 0.08; 
space_y = 1;
letter_fs = 8;
letter_arrow = 7.5;
letter_panel = 12;
arrLen = 0.4;
bar_yspace = 1;
bcd_w = (mat_h + bar_yspace/2 - 2*space_y)/3;

lineClrs = [0.9 0 0.9; 
    0 0.9 0;
    0 0 0.9
    0.1 0.1 0.1];
lineNms = {'m','g','b','k'};

ax_mat = axes('units','centimeters','position',[margin_left+space_x margin_bot+Bheight+bar_yspace mat_w mat_h]) ;
imagesc(ax_mat,turbidityInvOrgOrder')
xticks([]);
yticks([]);
axis([0.5, size(turbidityInvOrgOrder,1)+0.5, 0.5, size(turbidityInvOrgOrder,2)+0.5])        
colormap(flipud(gray))
xlabel('Bacteria')
ylabel('Phages')
caxis([min(min(turbidityOrgInvOrderNoNan)),max(max(turbidityOrgInvOrderNoNan))])
% plot X on nans:
hold on
[nanRt, nanCt] = find(isnan(turbidityInvOrgOrder'));
plot(ax_mat,[nanCt'-0.5;nanCt'+0.5],[nanRt'-0.5;nanRt'+0.5],'color',[0.2 0.2 0.2]);
plot(ax_mat,[nanCt'-0.5;nanCt'+0.5],[nanRt'+0.5;nanRt'-0.5],'color',[0.2 0.2 0.2]);

text(-2.5,0,'a','color','black','FontSize',letter_panel)

% Add colorbars
barH = mat_h/30;
barW = mat_w/3;
axes_clrBarInf = axes('units','centimeters','position',[margin_left+space_x margin_bot+Bheight*0.9 barW barH]);
turbBar = linspace(min(min(turbidityOrgInvOrderNoNan)),max(max(turbidityOrgInvOrderNoNan)),100);
imagesc(turbBar,'Parent',axes_clrBarInf)
colormap(axes_clrBarInf,flipud(gray))
set(gca,'XTick',1:49:100,'XTickLabel',sprintfc('%0.1f',[turbBar(1)  turbBar(50) turbBar(100)])...
    ,'YTick',[],'fontsize',8)
title('1-Turbidity')

% B: plot host switch examples over same/diff replicate and WT bacteria
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
BacIsolates = readtable(BacIsolatesTable);
PhgIsolates = readtable(PhgIsolatesTable);

PhgPlqSmpls = {'1_Cont-08-A','1_Cont-07-B','1_T7_WT_1'};
BacPlqSmpls = {'Cont-05-A','UR-S07-A','MG-Y_1'};

% Location in the matrix figure (original order):
HSline(1) = size(turbidity,2) - find(cellfun(@(x) contains(PhgPlqSmpls{3},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')),1);
HSline(2) =  size(turbidity,2) - find(cellfun(@(x) contains(PhgPlqSmpls{2},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')));
HSline(3) = size(turbidity,2) - find(cellfun(@(x) contains(PhgPlqSmpls{1},x), erase(PhgIsolates.SequenceNms,'Sample_Phage')));
HSline(4) = find(cellfun(@(x) contains(BacPlqSmpls{3},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')),1,'last');
HSline(5) = find(cellfun(@(x) contains(BacPlqSmpls{1},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')));
HSline(6) = find(cellfun(@(x) contains(BacPlqSmpls{2},x), erase(BacIsolates.SequenceNms,'Sample_Bac_')));

figPos = get(gcf,'position');
arrowString = {'WT','HR','HS'};
xArr = [];
yArr = [];
for i = 3:-1:1
    xArr(i,:) = [margin_left+space_x+mat_w+arrLen ...
        margin_left+space_x+mat_w ]/figPos(3);
    yArr(i,:) = [margin_bot+Bheight+bar_yspace+(HSline(i)+0.6)*pxl ...
        margin_bot+Bheight+bar_yspace+(HSline(i)+0.6)*pxl]/figPos(4);
    annotation('textarrow', xArr(i,:),yArr(i,:),'string',arrowString{i},...
        'HorizontalAlignment','center','fontsize',letter_arrow);
end

ax_sctB = axes('Units','centimeters','Position',...
            [margin_left+space_x*22+mat_w+arrLen margin_bot+Bheight+bar_yspace+mat_h-bcd_w bcd_w bcd_w]);
scatter(ax_sctB,infectivityNoNan(:),turbidityOrgInvOrderNoNan(:),5,'filled')
hold on
box(ax_sctB,'on')

xlabel('Infectivity Score')
ylabel('1-Turbidity')
xlim([0 0.75])
ylim([0 0.75]);
set(gca,'fontsize',letter_fs);
text(-0.18,0.75,'b','color','black','FontSize',letter_panel)

% Calculate R-square of correlation between WT replicate measurements
ax_C = axes('Units','centimeters','Position',...
            [margin_left+space_x*22+mat_w+arrLen margin_bot+Bheight+bar_yspace+mat_h-2*bcd_w-space_y bcd_w bcd_w]);       
bacWTpos = find(cellfun(@(x) contains(x,'MG'), BacIsolates.SequenceNms));
phgWTpos = find(cellfun(@(x) contains(x,'T7'), PhgIsolates.SequenceNms));

bacWTcpls = [bacWTpos(1),bacWTpos(2);
    bacWTpos(2),bacWTpos(3);
    bacWTpos(3),bacWTpos(1)];
phgWTcpls = [phgWTpos(1),phgWTpos(2);
    phgWTpos(2),phgWTpos(3);
    phgWTpos(3),phgWTpos(4);
    phgWTpos(4),phgWTpos(1)];
cplInd = [1,2,3,1];

for cp = size(bacWTcpls,1):-1:1
    [p,N] = fit(infectivity(bacWTcpls(cp,1),:)',infectivity(bacWTcpls(cp,2),:)','poly1');
    plot(ax_C,infectivity(bacWTcpls(cp,1),:),infectivity(bacWTcpls(cp,2),:),'.','color',lineClrs(cp,:))
    hold on
    plot(p,lineNms{cp})
    Rbac(cp) = round(N.rsquare,2);
    cplNm{cp} = [num2str(cplInd(cp)) 'vs.' num2str(cplInd(cp+1)) ];
    legend('hide')
    leg = [cplNm{cp} ' R^2=' num2str(round(N.rsquare,2))];
    text(0.02,0.77-cp*0.06,leg,'color',lineClrs(cp,:),'fontsize',7);
end
axis square
xlim([0 0.75])
ylim([0 0.75])
xlabel('WT Replicate 1')
ylabel('WT Replicate 2')
text(-0.18,0.75,'c','color','black','FontSize',letter_panel)


ax_D = axes('Units','centimeters','Position',...
        [margin_left+space_x*22+mat_w+arrLen margin_bot+Bheight+bar_yspace+mat_h-3*bcd_w-2*space_y bcd_w bcd_w]);
 
cplInd1 = [1,2,3,4,1];
for cp = size(phgWTcpls,1):-1:1
    [p,N] = fit(infectivity(:,phgWTcpls(cp,1)),infectivity(:,phgWTcpls(cp,2)),'poly1');
    plot(ax_D,infectivity(:,phgWTcpls(cp,1)),infectivity(:,phgWTcpls(cp,2)),'.','color',lineClrs(cp,:))
    hold on
    plot(p,lineNms{cp})
    Rphg(cp) = N.rsquare;    
    cplNm2{cp} = [num2str(cplInd1(cp)) 'vs.' num2str(cplInd1(cp+1)) ];
    legend('hide')
    leg = [cplNm2{cp} ' R^2=' num2str(round(N.rsquare,2))];
    text(0.02,0.77-cp*0.06,leg,'color',lineClrs(cp,:),'fontsize',7);
end
axis square
xlim([0 0.75])
ylim([0 0.75])
xlabel('WT Replicate 1')
ylabel('WT Replicate 2')

text(-0.18,0.75,'d','color','black','FontSize',letter_panel)
print([figure_location f 'SuppFigure6'],'-dpng','-r300');