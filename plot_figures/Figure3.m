%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply Lasso regression to identify the   %%%
%%% most significant mutations               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Printing Figure 3 - Lasso\n');
f = filesep;
fig_data_location = 'output_files';
fig_location = 'Figs';
pwClrMapFile = [fig_data_location f 'pathwayClrMap.mat'];
clrTblFile = ['..' f 'input_files' f 'PathwayTables' f 'colorBlindPalete.xlsx'];
%% parameters
plotLasso = 1;
unify = 1; % unify co-ocuring mutations
unify_same_genome = 1;
include_coverage = 1; % only in bacteria, FIRST RUN COVERAGE! FOR DPNN
include_syn = 0;
aa_names = 1; % take amino-acid mutation names
cv = 5;
nLambda = 200;
%% source tables
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 
BacGenotypeTable = ['..' f 'Tables' f 'Bac_genotypes_uni.xlsx']; 
PhgGenotypeTable = ['..' f 'Tables' f 'Phg_genotypes_uni.xlsx']; 
pw_loc = ['..' f 'input_files' f 'PathwayTables'];

%% Read source data
phenotypesOrg = table2array(readtable(phenotypeTable)); 
genInfoBac = readtable(BacGenotypeTable,'PreserveVariableNames',true);
bacCol1 = find(strcmp(genInfoBac.Properties.VariableNames,'genotype_ 1'));
bacColEnd = find(strcmp(genInfoBac.Properties.VariableNames,['genotype_' num2str( size(phenotypesOrg,1))] ));
genotype_bac2 = table2array(genInfoBac(:,bacCol1:bacColEnd));

genInfoPhg = readtable(PhgGenotypeTable,'PreserveVariableNames',true);
phgCol1 = find(strcmp(genInfoPhg.Properties.VariableNames,'genotype_  1'));
phgColEnd = find(strcmp(genInfoPhg.Properties.VariableNames,['genotype_' num2str( size(phenotypesOrg,2))] ));
genotype_phg2 = table2array(genInfoPhg(:,phgCol1:phgColEnd)); 

BacIsolates = readtable(BacIsolatesTable);
isonms_bac = BacIsolates.SequenceNms;
PhgIsolates = readtable(PhgIsolatesTable);
isonms_phg =  PhgIsolates.SequenceNms;

bacWTinfectiviy = median(reshape(phenotypesOrg(95:97,:),[],1));

%% Build predictors matrix and phenotype vector
if ~include_syn
    genotype_bac1 = genotype_bac2(~strcmp(genInfoBac.mutType,'SNP-Syn'),:);
    genotype_phg1 = genotype_phg2(~strcmp(genInfoPhg.mutType,'SNP-Syn'),:);
    cnvs = find(strcmp(genInfoBac.mutType(~strcmp(genInfoBac.mutType,'SNP-Syn')),'CN'));
    mutnms_bac_cnv = genInfoBac.mutName(strcmp(genInfoBac.mutType,'CN') & ~strcmp(genInfoBac.mutType,'SNP-Syn'));
    mutnms_bac_no_cnv = genInfoBac.mutName(~strcmp(genInfoBac.mutType,'CN') & ~strcmp(genInfoBac.mutType,'SNP-Syn'));
    mutnms_phg = genInfoPhg.mutName(~strcmp(genInfoPhg.mutType,'SNP-Syn'));
else
    genotype_bac1 = genotype_bac2;
    genotype_phg1 = genotype_phg2;
    cnvs = find(strcmp(genInfoBac.mutType,'CN'));
    mutnms_bac_no_cnv = genInfoBac.mutName(setdiff(1:size(genotype_bac1,1),cnvs));
    mutnms_bac_cnv = genInfoBac.mutName(cnvs);
    mutnms_phg = genInfoPhg.mutName;
end

[uggenotype_bac,iugbac,rugbac ] = unique(genotype_bac1','rows','stable');
uggenotype_bac = uggenotype_bac';
[uggenotype_phg,iugphg,rugphg ] = unique(genotype_phg1','rows','stable');
uggenotype_phg = uggenotype_phg';
uphenotypebac = zeros(max(rugbac),size(phenotypesOrg,2));
for rb = 1:max(rugbac)
    uphenotypebac(rb,:) =  mean(phenotypesOrg(rugbac==rb,:),1,'omitnan');
end
uphenotypephg = zeros(max(rugbac),max(rugphg));
for rp = 1:max(rugphg)
    uphenotypephg(:,rp) =  mean(uphenotypebac(:,rugphg==rp),2,'omitnan');
end
genotype_bac = uggenotype_bac;
genotype_phg = uggenotype_phg;
phenotypes = uphenotypephg;
isonms_bac = isonms_bac(iugbac);
isonms_phg = isonms_phg(iugphg);

[ugenotype_bac_no_amp,ubac,rbac] = unique(genotype_bac(setdiff(1:size(genotype_bac,1),cnvs),:),'rows','stable');
[ugenotype_phg,uphg,rphg] = unique(genotype_phg,'rows','stable');

for n = 1:numel(ubac)
    r = find(rbac==n);
    umutnms_bac_no_amp(n) = mutnms_bac_no_cnv(r(1));
    for w = 2:numel(r)
        umutnms_bac_no_amp{n} = [umutnms_bac_no_amp{n} '+' mutnms_bac_no_cnv{r(w)}];
    end
end
ugenotype_bac_amp = double(genotype_bac(cnvs,:));
ugenotype_bac = [ugenotype_bac_no_amp; ugenotype_bac_amp];
umutnms_bac = [umutnms_bac_no_amp, mutnms_bac_cnv'];
for n = 1:numel(uphg)
    r = find(rphg==n);
    umutnms_phg(n) = mutnms_phg(r(1));
    for w = 2:numel(r)
        umutnms_phg{n} = [umutnms_phg{n} '+' mutnms_phg{r(w)}];
    end
end

predictors_all = [reshape(repmat(ugenotype_bac,[size(phenotypes,2) 1]),...
    [size(ugenotype_bac,1) size(phenotypes,1)*size(phenotypes,2)])',...
    repmat((ugenotype_phg'),[size(phenotypes,1) 1])];
predictor_nms_all = [umutnms_bac' ; umutnms_phg'];    

%% remove intermidiate infectivity, No phenotype and No mut lines:
phen_vector_all = reshape(phenotypes',[size(phenotypes,1)*size(phenotypes,2) 1]);
nodata = find(isnan(phen_vector_all)); % find NaN interactions
remove_ii= nodata;
phen_vec_cln = phen_vector_all(setdiff(1:end,remove_ii));
nomut = find(sum(predictors_all,1)==0); % find mutations with no isolates carrying them
predictors_noNaN = predictors_all(setdiff(1:end,remove_ii),:);
pred_cln = predictors_noNaN(:,setdiff(1:end,nomut));
pred_nms_cln = predictor_nms_all(setdiff(1:end,nomut));
%% run Lasso
lassoName = sprintf('lasso_cv%d_nLambda%d_unify%d_uniGenomes%d_includeCoverage%d_includeSyn%d.mat',...
    cv,nLambda,unify,unify_same_genome,include_coverage,include_syn);

if ~isfile([fig_data_location f lassoName])
    fprintf('Number of predictors: %d, Number of phenotypes:%d, running Lasso...\n',numel(pred_nms_cln),size(phen_vec_cln,1));
    [B,FitInfo] = lasso(double(pred_cln),phen_vec_cln,'CV',cv,'NumLambda',nLambda,...
        'PredictorNames',pred_nms_cln);

    save([fig_data_location f lassoName],'B','FitInfo');
else
    fprintf('loading Lasso... Number of phenotypes:%d\n',size(phen_vec_cln,1));
    load([fig_data_location f lassoName]);
end
if plotLasso
    lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');
    figure(120); clf; plot(1:200,FitInfo.MSE); hold on; errorbar(1:200,FitInfo.MSE,FitInfo.SE)
    xline(FitInfo.IndexMinMSE,'-')
    xline(FitInfo.Index1SE,'--')
end
idxLambdaChosen = FitInfo.Index1SE; % group meeting/seminar: FitInfo.Index1SE;
if idxLambdaChosen == nLambda % of no minimum, take 3% of predictors
    for i = 1:nLambda
        if numel(find(B(:,i)~=0)) > size(pred_cln,2)*0.03 && numel(find(B(:,i)~=0)) < size(pred_cln,2)*0.032
            idxLambdaChosen = i;
        end
    end
end
lassoPredInd = find(B(:,idxLambdaChosen)~=0); %FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0);
fprintf('%d predictors remained\n',numel(lassoPredInd))
% Plot sorted coeffs:
[sortedCoeffs,Bi] = sort(B(:,idxLambdaChosen),'descend');
pred_nms = pred_nms_cln(Bi(sortedCoeffs~=0));
% divide to bac and phg
[~,bacMutInd] = ismember(umutnms_bac, pred_nms);
bacMutInd = sort(bacMutInd(bacMutInd>0));
[~,phgMutInd] = ismember(umutnms_phg,pred_nms);
phgMutInd = sort(phgMutInd(phgMutInd>0));

%% Plot Figure:
clrTbl = readtable(clrTblFile);
bacPwsClrMap = hex2rgb([clrTbl.Hex ; clrTbl.Hex]);
phgPwsClrMap = hex2rgb(clrTbl.Hex);
pw_fs = 9;
panel_fs = 10;
mutfs = 7;
barW = 8.6;
barH = 0.28;
space_y = 0.5;
space_x = 0.3;
margin_x = 0.3;
margin_y = 1.2;
tit_fs = 10;
label_fs = 9;
xlim_val_bac = 0.6; % 6.1;
xlim_val_phg = 0.3; %4.1;
text_space = xlim_val_bac/60;
addMultipleManually = 0;
darkGray = [0.2 0.2 0.2];
lightGray =  [0.6 0.6 0.6];
leg_tit_h = 1; 0.55;
leg_phgTit_x = 0.09;
H = max(numel(phgMutInd),numel(bacMutInd))*barH+2*margin_y;
[bac_pw_all_multi,bac_pw] = readPathways(pred_nms(bacMutInd),1,pw_loc,'Uni');
[phg_pw_all_multi,phg_pw] = readPathways(pred_nms(phgMutInd),2,pw_loc,'Uni');

if any(contains(bac_pw,'uncharacterized'))
    bacPwsClrMap(strcmp(bac_pw,'uncharacterized'),:) = lightGray;
end
bacPwsClrMap(strcmp(bac_pw,'multiple'),:) = darkGray;
if any(contains(phg_pw,'uncharacterized'))
    phgPwsClrMap(strcmp(phg_pw,'uncharacterized'),:) = lightGray;
end
if any(strcmp(phg_pw,'multiple'))
    phgPwsClrMap(strcmp(phg_pw,'multiple'),:) =darkGray;
else
    addMultipleManually = 1;
end

sortedCoeffsNo0 = sortedCoeffs(sortedCoeffs~=0);
figure(3); clf;
set(gcf,'units','centimeters','Position',[1 2 18 21]);

% Bacteria:
ax_bar_Bac = axes('units','centimeters','position',...
    [margin_x margin_y+(barH*(numel(phgMutInd)-numel(bacMutInd)))...
    barW barH*numel(bacMutInd)]);
for i = 1:numel(bacMutInd)
    bhB = barh(i,sortedCoeffsNo0(bacMutInd(i)));
    pw = bac_pw_all_multi{i};
    xMutNm = pred_nms{bacMutInd(i)};
    xMutNmSplit = strsplit(xMutNm,'+');
    mutNameText = [];
    if numel(xMutNmSplit)<4
        for sm = 1:numel(xMutNmSplit)
            if size(pw,1)==1
                pwNm = pw;
            else
                pwNm = pw{sm};
            end
            nimcerMutNm = xMutNmSplit{sm};
            mutNameText = [mutNameText '\color[rgb]{' ...
                sprintf('%1.2f,%1.2f,%1.2f', bacPwsClrMap(strcmp(pwNm,bac_pw),:)) '}' nimcerMutNm ,'+'];

        end
    else
        mutNameText = '\color[rgb]{0,0,0}4 mutations or more ';
        bac_pw_all_multi{i} = {''};

    end
    mutNameText = mutNameText(1:end-1); % remove ',' from end
    if size(pw,1)==1 % if it's a single mutation          
        bhB.FaceColor = bacPwsClrMap(strcmp(pw,bac_pw),:);
        if find(strcmp(pw,bac_pw))>10 && ~strcmp(pw,'multiple')
            hatchfill2(bhB,'cross','HatchAngle',45);
        end
    else % if mutations appear together
        bhB.FaceColor = darkGray;
    end

    if sortedCoeffsNo0(bacMutInd(i))<0
        x = text_space;
        align = 'left';
    else
        x = -text_space;
        align = 'right';
    end
    text(x,i,mutNameText,'HorizontalAlignment',align,'fontsize',mutfs)
    hold on
end   
set(gca,'YTick',[],'YTickLabel',[],'FontSize',label_fs,...
    'XTick',-xlim_val_bac:0.2:xlim_val_bac,'XTickLabel',sprintfc('%0.1f',-(xlim_val_bac):0.2:(xlim_val_bac)))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',8);
xlim([-xlim_val_bac,xlim_val_bac]);  
ylim([0 numel(sortedCoeffsNo0(bacMutInd))+1]);
xlabel('Association with Infectivity','FontSize',label_fs)
title('Bacteria','FontSize',tit_fs)
text(-xlim_val_bac,numel(sortedCoeffsNo0(bacMutInd))+2.3,'a','color','black','FontSize',panel_fs);
dx = margin_x/(barW*2+space_x+margin_x*2);
bardX = barW/(barW*2+space_x+margin_x*2);
spaceX = bardX/8;
y_arr = 0.95;
test_ht = 0.02;
annotation('arrow',[dx+bardX/2-spaceX dx+spaceX],[y_arr y_arr])
annotation('textbox',[dx+spaceX+0.01 y_arr+0.005 bardX/2-2*spaceX test_ht],'String',...
    'Resistance','HorizontalAlignment','center','EdgeColor','none')
annotation('arrow',[dx+bardX/2+spaceX dx+bardX-spaceX],[y_arr y_arr])
annotation('textbox',[dx+bardX/2+spaceX-0.01 y_arr+0.005 bardX/2-2*spaceX test_ht],...
    'String','Sensitivity','HorizontalAlignment','center','EdgeColor','none')
% Phages:
ax_bar_Phg = axes('units','centimeters','position',...
    [margin_x+space_x+barW margin_y barW barH*numel(phgMutInd)]);
for i = 1:numel(phgMutInd)
    bhP = barh(i,sortedCoeffsNo0(phgMutInd(i)));
    pw = phg_pw_all_multi{i};
    xMutNm = pred_nms{phgMutInd(i)};
    xMutNmSplit = strsplit(xMutNm,'+');
    mutNameText = [];
    if numel(xMutNmSplit)<4
        for sm = 1:numel(xMutNmSplit)
            if size(pw,1)==1
                pwNm = pw;
            else
                pwNm = pw{sm};
            end
            nimcerMutNm = xMutNmSplit{sm};
            mutNameText = [mutNameText '\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f',...
                phgPwsClrMap(strcmp(pwNm,phg_pw),:)) '}' nimcerMutNm '+'];
        end
    else
    mutNameText = '\color[rgb]{0,0,0}4 mutations or more ';
    phg_pw_all_multi{i} = {''};
    end 
    mutNameText = mutNameText(1:end-1); % remove ',' from end
    if size(pw,1)==1 % if it's a single mutation          
        bhP.FaceColor = phgPwsClrMap(strcmp(pw,phg_pw),:);
    else % if mutations appear together
        bhP.FaceColor = darkGray;
    end
    if sortedCoeffsNo0(phgMutInd(i))<0
        x = text_space;
        align = 'left';
    else
        x = -text_space;
        align = 'right';
    end
    text(x,i,mutNameText,'HorizontalAlignment',align,'fontsize',mutfs)
    hold on
end 
set(gca,'YTick',[],'YTickLabel',...
    [],'FontSize',label_fs,'XTick',-(xlim_val_phg):0.1:(xlim_val_phg),...
    'XTickLabel',sprintfc('%0.1f',-(xlim_val_phg):0.1:(xlim_val_phg)))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',8);
xlim([-xlim_val_phg,xlim_val_phg])
ylim([0 numel(sortedCoeffsNo0(phgMutInd))+1]);
xlabel('Association with Infectivity','FontSize',label_fs)
title('Phages','FontSize',tit_fs)
text(-xlim_val_phg,numel(sortedCoeffsNo0(phgMutInd))+2.3,'b','color','black','FontSize',panel_fs)

% Print Pathway legend:
ax_legend_Bac = axes('units','centimeters','position',...
    [margin_x 0 ...
    barW/2 barH*(numel(phgMutInd)-numel(bacMutInd))]);
hold on
set(gca,'XTick',[],'YTick',[],'color','none','Visible','off');
xlim([0 1]);
ylim([0 1]);
[A,B] = groupcounts(cat(1,bac_pw_all_multi{:}));
[~,freqPwInd] = sort(A,'descend');
freqPw = B(freqPwInd);
bacLegend = unique(cat(1,bac_pw_all_multi{:}));
bacLegend(strcmp(bacLegend,'')) = [];
bacLegend(strcmp(bacLegend,'multiple')) = [];
freqPw(strcmp(freqPw,'')) = [];
freqPw(strcmp(freqPw,'multiple')) = [];
[~,bacOrd] = ismember(freqPw,bacLegend);
bacLegend = bacLegend(bacOrd);
[~,clrCode] = ismember(bacLegend,bac_pw);
text(0,leg_tit_h,'\underline{Bacterial Pathways}','color','k','fontSize',pw_fs,'interp','latex','fontname','Arial')
for pLeg = 1:numel(bacLegend)
    if pLeg<=9
        text(ax_legend_Bac,0.07,leg_tit_h-0.045-0.092*pLeg,bacLegend{pLeg},'color',bacPwsClrMap(clrCode(pLeg),:),...
            'fontSize',pw_fs)
        rectangle(ax_legend_Bac,'Position',[0,leg_tit_h-0.045-0.092*pLeg-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
        if clrCode(pLeg)>10
            plot(ax_legend_Bac,[0 0.05],[leg_tit_h-0.045-0.092*pLeg-0.05 leg_tit_h-0.045-0.092*pLeg+0.05],'-k')
        end
    else
        text(ax_legend_Bac,0.77,leg_tit_h-0.045-0.092*(pLeg-9),bacLegend{pLeg},'color',bacPwsClrMap(clrCode(pLeg),:),'fontSize',pw_fs)
        rectangle(ax_legend_Bac,'Position',[0.70,leg_tit_h-0.045-0.092*(pLeg-9)-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
        if clrCode(pLeg)>10
            plot(ax_legend_Bac,[0.70 0.70+0.05],[leg_tit_h-0.045-0.092*(pLeg-9)-0.05 leg_tit_h-0.045-0.092*(pLeg-9)+0.05],'-k')
        end
    end
end

ax_legend_phg = axes('units','centimeters','position',...
    [margin_x+barW/2+1.2 0 ...
    barW/2 barH*(numel(phgMutInd)-numel(bacMutInd))]);
set(gca,'XTick',[],'YTick',[],'color','none','Visible','off');
xlim([0 1]);
ylim([0 1]);
[A,B] = groupcounts(cat(1,phg_pw_all_multi{:}));
[~,freqPwInd] = sort(A,'descend');
freqPw = B(freqPwInd);
phgLegend = unique(cat(1,phg_pw_all_multi{:}));
phgLegend(strcmp(phgLegend,'')) = [];
phgLegend(strcmp(phgLegend,'multiple')) = [];
freqPw(strcmp(freqPw,'')) = [];
freqPw(strcmp(freqPw,'multiple')) = [];
[~,phgOrd] = ismember(freqPw,phgLegend);
phgLegend = phgLegend(phgOrd);
[~,clrCode] = ismember(phgLegend,phg_pw);
text(leg_phgTit_x,leg_tit_h,'\underline{Phage Pathways}','color','k','fontSize',...
    pw_fs,'interp','latex','fontname','Arial')
for pLeg = 1:numel(phgLegend)
    text(ax_legend_phg,leg_phgTit_x+0.07,leg_tit_h-0.045-0.092*pLeg,phgLegend{pLeg},...
        'color',phgPwsClrMap(clrCode(pLeg),:),'fontSize',pw_fs)
    rectangle(ax_legend_phg,'Position',[leg_phgTit_x,leg_tit_h-0.045-0.092*pLeg-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
end
if addMultipleManually
    text(ax_legend_phg,leg_phgTit_x+0.07,leg_tit_h-0.045-0.092*(pLeg+1),'multiple',...
        'color',darkGray,'fontSize',pw_fs)
    rectangle(ax_legend_phg,'Position',[leg_phgTit_x,leg_tit_h-0.045-0.092*(pLeg+1)-0.05,0.05, 0.1],...
        'FaceColor',darkGray,'EdgeColor','none');
end
ddx = (margin_x+barW)/(barW*2+space_x+margin_x*2);
eps = space_x/(barW*2+space_x+margin_x*2);
annotation('arrow',[ddx+dx+bardX/2-spaceX ddx+dx+spaceX],[y_arr y_arr])
annotation('textbox',[ddx+dx+spaceX+eps y_arr+0.008 bardX/2-2*spaceX test_ht],...
    'String','Lower Infectivity','HorizontalAlignment','center','EdgeColor','none')
annotation('arrow',[ddx+dx+bardX/2+spaceX ddx+dx+bardX-spaceX],[y_arr y_arr])
annotation('textbox',[ddx+dx+bardX/2+spaceX-eps y_arr+0.008 bardX/2-2*spaceX test_ht],...
    'String','Higher Infectivity','HorizontalAlignment','center','EdgeColor','none')

print([fig_location f 'Figure3'],'-dpng','-r300');
save(pwClrMapFile,'bacPwsClrMap','phgPwsClrMap','bac_pw','phg_pw');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print Figure 3 with full mutation names for reviewer comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bac_pw_all_multi_fn,bac_pw] = readPathways(pred_nms(bacMutInd),1,pw_loc,'Uni');
[phg_pw_all_multi_fn,phg_pw] = readPathways(pred_nms(phgMutInd),2,pw_loc,'Uni');
figure(3003); clf;
set(gcf,'units','centimeters','Position',[2 1 barW*2+space_x+margin_x*2 H]);
% Bacteria:
ax_bar_Bac = axes('units','centimeters','position',...
    [margin_x margin_y+(barH*(numel(phgMutInd)-numel(bacMutInd)))...
    barW barH*numel(bacMutInd)]);
for i = 1:numel(bacMutInd)
    bhB = barh(i,sortedCoeffsNo0(bacMutInd(i)));
    pw = bac_pw_all_multi_fn{i};
    xMutNm = pred_nms{bacMutInd(i)};
    xMutNmSplit = strsplit(xMutNm,'+');
    mutNameText = [];
        for sm = 1:numel(xMutNmSplit)
            if size(pw,1)==1
                pwNm = pw;
            else
                pwNm = pw{sm};
            end
            nimcerMutNm = xMutNmSplit{sm};
            mutNameText = [mutNameText '\color[rgb]{' ...
                sprintf('%1.2f,%1.2f,%1.2f', bacPwsClrMap(strcmp(pwNm,bac_pw),:)) '}' nimcerMutNm ,'+'];

        end

    mutNameText = mutNameText(1:end-1); % remove ',' from end
    if size(pw,1)==1 % if it's a single mutation          
        bhB.FaceColor = bacPwsClrMap(strcmp(pw,bac_pw),:);
        if find(strcmp(pw,bac_pw))>10
            hatchfill2(bhB,'cross','HatchAngle',45);
        end
    else % if mutations appear together
        bhB.FaceColor = darkGray;
    end

    if sortedCoeffsNo0(bacMutInd(i))<0
        x = text_space;
        align = 'left';
    else
        x = -text_space;
        align = 'right';
    end
    text(x,i,mutNameText,'HorizontalAlignment',align,'fontsize',mutfs)
    hold on
end   
set(gca,'YTick',[],'YTickLabel',[],'FontSize',label_fs,...
    'XTick',-xlim_val_bac:0.2:xlim_val_bac,'XTickLabel',-(xlim_val_bac):0.2:(xlim_val_bac))
xlim([-xlim_val_bac,xlim_val_bac]);  
ylim([0 numel(sortedCoeffsNo0(bacMutInd))+1]);
xlabel('Association with Infectivity','FontSize',label_fs)
title('Bacteria','FontSize',tit_fs)
text(-xlim_val_bac,numel(sortedCoeffsNo0(bacMutInd))+2.3,'a','color','black','FontSize',panel_fs);
dx = margin_x/(barW*2+space_x+margin_x*2);
bardX = barW/(barW*2+space_x+margin_x*2);
spaceX = bardX/8;
y_arr = 0.95;
test_ht = 0.02;
annotation('arrow',[dx+bardX/2-spaceX dx+spaceX],[y_arr y_arr])
annotation('textbox',[dx+spaceX+0.01 y_arr+0.005 bardX/2-2*spaceX test_ht],'String',...
    'Resistance','HorizontalAlignment','center','EdgeColor','none')
annotation('arrow',[dx+bardX/2+spaceX dx+bardX-spaceX],[y_arr y_arr])
annotation('textbox',[dx+bardX/2+spaceX-0.01 y_arr+0.005 bardX/2-2*spaceX test_ht],...
    'String','Sensitivity','HorizontalAlignment','center','EdgeColor','none')
% Phages:
ax_bar_Phg = axes('units','centimeters','position',...
    [margin_x+space_x+barW margin_y barW barH*numel(phgMutInd)]);
for i = 1:numel(phgMutInd)
    bhP = barh(i,sortedCoeffsNo0(phgMutInd(i)));
    pw = phg_pw_all_multi_fn{i};
    xMutNm = pred_nms{phgMutInd(i)};
    xMutNmSplit = strsplit(xMutNm,'+');
    mutNameText = [];
        for sm = 1:numel(xMutNmSplit)
            if size(pw,1)==1
                pwNm = pw;
            else
                pwNm = pw{sm};
            end
            nimcerMutNm = xMutNmSplit{sm};
            mutNameText = [mutNameText '\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f',...
                phgPwsClrMap(strcmp(pwNm,phg_pw),:)) '}' nimcerMutNm '+'];
        end

    mutNameText = mutNameText(1:end-1); % remove ',' from end
    if size(pw,1)==1 % if it's a single mutation          
        bhP.FaceColor = phgPwsClrMap(strcmp(pw,phg_pw),:);
    else % if mutations appear together
        bhP.FaceColor = darkGray;
    end
    if sortedCoeffsNo0(phgMutInd(i))<0
        x = text_space;
        align = 'left';
    else
        x = -text_space;
        align = 'right';
    end
    text(x,i,mutNameText,'HorizontalAlignment',align,'fontsize',mutfs)
    hold on
end 
set(gca,'YTick',[],'YTickLabel',...
    [],'FontSize',label_fs,'XTick',-(xlim_val_phg):0.1:(xlim_val_phg),...
    'XTickLabel',-(xlim_val_phg):0.1:(xlim_val_phg))
xlim([-xlim_val_phg,xlim_val_phg])
ylim([0 numel(sortedCoeffsNo0(phgMutInd))+1]);
xlabel('Association with Infectivity','FontSize',label_fs)
title('Phages','FontSize',tit_fs)
text(-xlim_val_phg,numel(sortedCoeffsNo0(phgMutInd))+2.3,'b','color','black','FontSize',panel_fs)

% Print Pathway legend:
ax_legend_Bac = axes('units','centimeters','position',...
    [margin_x 0 ...
    barW/2 barH*(numel(phgMutInd)-numel(bacMutInd))]);
hold on
set(gca,'XTick',[],'YTick',[],'color','none','Visible','off');
xlim([0 1]);
ylim([0 1]);
[A,B] = groupcounts(cat(1,bac_pw_all_multi_fn{:}));
[~,freqPwInd] = sort(A,'descend');
freqPw = B(freqPwInd);
bacLegend = unique(cat(1,bac_pw_all_multi_fn{:}));
bacLegend(strcmp(bacLegend,'')) = [];
bacLegend(strcmp(bacLegend,'multiple')) = [];
freqPw(strcmp(freqPw,'')) = [];
freqPw(strcmp(freqPw,'multiple')) = [];
[~,bacOrd] = ismember(freqPw,bacLegend);
bacLegend = bacLegend(bacOrd);
[~,clrCode] = ismember(bacLegend,bac_pw);
text(0,leg_tit_h,'\underline{Bacterial Pathways}','color','k','fontSize',pw_fs,'interp','latex')
for pLeg = 1:numel(bacLegend)
    if pLeg<=9
        text(ax_legend_Bac,0.08,leg_tit_h-0.045-0.092*pLeg,bacLegend{pLeg},'color',bacPwsClrMap(clrCode(pLeg),:),...
            'fontSize',pw_fs)
        rectangle(ax_legend_Bac,'Position',[0,leg_tit_h-0.045-0.092*pLeg-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
        if clrCode(pLeg)>10
            plot(ax_legend_Bac,[0 0.05],[leg_tit_h-0.045-0.092*pLeg-0.05 leg_tit_h-0.045-0.092*pLeg+0.05],'-k')
        end
    else
        text(ax_legend_Bac,0.67,leg_tit_h-0.045-0.092*(pLeg-9),bacLegend{pLeg},'color',bacPwsClrMap(clrCode(pLeg),:),'fontSize',pw_fs)
        rectangle(ax_legend_Bac,'Position',[0.59,leg_tit_h-0.045-0.092*(pLeg-9)-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
        if clrCode(pLeg)>10
            plot(ax_legend_Bac,[0.59 0.59+0.05],[leg_tit_h-0.045-0.092*(pLeg-9)-0.05 leg_tit_h-0.045-0.092*(pLeg-9)+0.05],'-k')
        end
    end
end

ax_legend_phg = axes('units','centimeters','position',...
    [margin_x+barW/2+1.2 0 ...
    barW/2 barH*(numel(phgMutInd)-numel(bacMutInd))]);
set(gca,'XTick',[],'YTick',[],'color','none','Visible','off');
xlim([0 1]);
ylim([0 1]);
[A,B] = groupcounts(cat(1,phg_pw_all_multi_fn{:}));
[~,freqPwInd] = sort(A,'descend');
freqPw = B(freqPwInd);
phgLegend = unique(cat(1,phg_pw_all_multi_fn{:}));
phgLegend(strcmp(phgLegend,'')) = [];
phgLegend(strcmp(phgLegend,'multiple')) = [];
freqPw(strcmp(freqPw,'')) = [];
freqPw(strcmp(freqPw,'multiple')) = [];
[~,phgOrd] = ismember(freqPw,phgLegend);
phgLegend = phgLegend(phgOrd);
[~,clrCode] = ismember(phgLegend,phg_pw);
text(leg_phgTit_x,leg_tit_h,'\underline{Phage Pathways}','color','k','fontSize',pw_fs,'interp','latex')
for pLeg = 1:numel(phgLegend)
    text(ax_legend_phg,leg_phgTit_x+0.1,leg_tit_h-0.045-0.092*pLeg,phgLegend{pLeg},...
        'color',phgPwsClrMap(clrCode(pLeg),:),'fontSize',pw_fs)
    rectangle(ax_legend_phg,'Position',[leg_phgTit_x,leg_tit_h-0.045-0.092*pLeg-0.05,0.05, 0.1],'FaceColor',bacPwsClrMap(clrCode(pLeg),:),...
            'EdgeColor','none');
end

ddx = (margin_x+barW)/(barW*2+space_x+margin_x*2);
eps = space_x/(barW*2+space_x+margin_x*2);
annotation('arrow',[ddx+dx+bardX/2-spaceX ddx+dx+spaceX],[y_arr y_arr])
annotation('textbox',[ddx+dx+spaceX+eps y_arr+0.005 bardX/2-2*spaceX test_ht],'String','Lower Infectivity','HorizontalAlignment','center','EdgeColor','none')
annotation('arrow',[ddx+dx+bardX/2+spaceX ddx+dx+bardX-spaceX],[y_arr y_arr])
annotation('textbox',[ddx+dx+bardX/2+spaceX-eps y_arr+0.005 bardX/2-2*spaceX test_ht],'String','Higher Infectivity','HorizontalAlignment','center','EdgeColor','none')

print([fig_location f 'Figure3_fullMutationNames'],'-dpng','-r300');
print([figure_location f 'Figure3'],'-depsc2')

