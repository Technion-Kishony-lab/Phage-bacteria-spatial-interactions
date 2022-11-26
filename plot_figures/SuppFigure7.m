%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate nestedness and modularity of each replicate with BiMat %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('BiMat...\n');
tic
f = filesep;
% data location:
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
bp_figure_location = 'Figs';
bp_location = 'output_files';
%% Bimat params:
BiMatMatrix = 'original'; % options: 'original' / 'sigMutsStruct';
binThresh = -2.5;  % for Bimat binarization(only original)
removeAllZeros = 1; % remove all 0 rows and columns 
randomMats = 5000;
randomMatsPerPlate = 10000;
NullModel = 'EQUIPROBABLE';
optAlgo = 'AdaptiveBrim';
pltMap = [0.7 0.7 0.7; brewermap(4,'Pastel2')];
pltMap(2,:) = pltMap(2,:)-[0 0 0.1]; % colors in original colormap are too similar
pltMap(3,:) = pltMap(3,:)-[0 0 0.1];

%% load phenotypes and isolate info
BacIsolates = readtable(BacIsolatesTable);
PhgIsolates = readtable(PhgIsolatesTable);
mat = table2array(readtable(phenotypeTable)); 
tmp = unique(sort(reshape(mat,1,[])));
log_param = tmp(2);
phenMatLog = log(mat+log_param); 

%% Arrange matrix
% Remove NaNs, by averaging 
nanMeanArea = 2; 
[a,b] = find(isnan(mat));
for i = 1:numel(a)
    mat(a(i),b(i)) = mean(nanmean(mat(max(a(i)-nanMeanArea,1):...
        min(a(i)+nanMeanArea,size(phenMatLog,1)),max(b(i)-nanMeanArea,1):min(b(i)+nanMeanArea,size(phenMatLog,2)))));
end  

% Binarize :
mat_log = log(mat + log_param);
matrix_bp = zeros(size(mat_log));
matrix_bp(mat_log <= binThresh) = 0;
matrix_bp(mat_log > binThresh) = 1;
% run BiMat: 
bp_name = sprintf('%s%sbp_%s_%s_%d_%s_remove0_%0.1f_%d_Matrix.mat',bp_location,f,...
    NullModel,optAlgo,randomMats, BiMatMatrix, binThresh,removeAllZeros);
bp_name = bp_name(~isspace(bp_name));

%% BiMat, per plate analysis:
bp_name_per_plate = sprintf('%s %s bpPerPlate_%s_%s_%d_%s_%0.1f_remove0_%d_Matrix.mat',bp_location,f,...
    NullModel,optAlgo,randomMatsPerPlate, BiMatMatrix,binThresh,removeAllZeros);
bp_name_per_plate = bp_name_per_plate(~isspace(bp_name_per_plate));
nbins = 20;
if ~exist(bp_name_per_plate,'file')
    for pl = 4:-1:1   
        fprintf('plate %d\n',pl)
        indBac = [find(BacIsolates.replicateNum==pl); find(BacIsolates.replicateNum==0,1)];
        indPhg = [find(PhgIsolates.replicateNum==pl); find(PhgIsolates.replicateNum==0,1)];
        matrix_bp_plate = matrix_bp(indBac,indPhg);
        bpPerPlate(pl).bp = Bipartite(matrix_bp_plate);
        % Print general properties:
        bpPerPlate(pl).bp.printer.PrintGeneralProperties();
        % Optimization algorithm:
        switch optAlgo 
            case 'LeadingEigenvector'
                bpPerPlate(pl).bp.community = LeadingEigenvector(bpPerPlate(pl).bp.matrix);
                bpPerPlate(pl).bp.community.DoKernighanLinTunning = true;
            case 'LPBrim'
                bpPerPlate(pl).bp.community = LPBrim(bpPerPlate(pl).bp.matrix);
            case 'AdaptiveBrim'
                bpPerPlate(pl).bp.community = AdaptiveBrim(bpPerPlate(pl).bp.matrix);
        end
        % Calculate the modularity
        bpPerPlate(pl).bp.community.Detect();
        % Calcualte the nestedness:
        bpPerPlate(pl).bp.nestedness.Detect();
        switch NullModel
            case 'FIXED'
                bpPerPlate(pl).bp.statistics.DoCompleteAnalysis(randomMatsPerPlate, @NullModels.FIXED);
            case 'AVERAGE'
                bpPerPlate(pl).bp.statistics.DoCompleteAnalysis(randomMatsPerPlate, @NullModels.AVERAGE);
            case 'EQUIPROBABLE'
                bpPerPlate(pl).bp.statistics.DoCompleteAnalysis(randomMatsPerPlate, @NullModels.EQUIPROBABLE); 
        end
        bpPerPlate(pl).bp.printer.PrintStructureStatistics();
    end
    save(bp_name_per_plate,'bpPerPlate');
else
    load(bp_name_per_plate)
end
%% Plot BiMat per plate:
sz = 3;
margin = 1;
space = 0.7;
spaceAB = 0.5;
figure(1006); clf;
set(gcf,'units','centimeters','position',[1 1 sz*4+margin*2+space*3 sz*3+space+margin*4+spaceAB])
% Nestedness
for pl = 4:-1:1
    axes('units','centimeters','position', [margin+(pl-1)*(space+sz) ...
        margin+space+sz/2+2*sz+spaceAB sz sz])
    bpPerPlate(pl).bp.plotter.use_labels = 0;
    bpPerPlate(pl).bp.plotter.PlotNestedMatrix();
    title(['Replicate' num2str(pl)]);
    if pl==1
        annotation('textbox',[0.02 0.87 0.05 0.05],'string','a','EdgeColor','none','FontSize',12)
    end
end
for pl = 4:-1:1
    axes('units','centimeters','position', [margin+(pl-1)*(space+sz) ...
        margin+2*sz+spaceAB sz sz/2])
    hstNest(pl).hst = histogram(bpPerPlate(pl).bp.statistics.N_values.random_values,'FaceColor',pltMap(pl+1,:));
    ylim([0 1000]);
    pval_nest = sum(bpPerPlate(pl).bp.statistics.N_values.random_values>=bpPerPlate(pl).bp.nestedness.N)/randomMatsPerPlate;
    pval_antinest = sum(bpPerPlate(pl).bp.statistics.N_values.random_values<=bpPerPlate(pl).bp.nestedness.N)/randomMatsPerPlate;
    xline(bpPerPlate(pl).bp.nestedness.N,'r-','linewidth',2);
    if pval_nest>=1E-4
        tit_nest = sprintf('Nestedness, p= %0.2f',pval_nest);
    else
        tit_nest = sprintf('Nestedness, p<10^{-4}');
    end
    title(tit_nest);
    if pl==1
        annotation('textbox',[0.02 0.66 0.05 0.05],'string','b','EdgeColor','none','FontSize',12)
    end
end
% Modularity
for pl = 4:-1:1
    axes('units','centimeters','position', [margin+(pl-1)*(space+sz) margin+space+sz/2 sz sz])
    bpPerPlate(pl).bp.plotter.use_labels = 0;
    bpPerPlate(pl).bp.plotter.PlotModularMatrix();
    title(['Replicate' num2str(pl)]);
    if pl==1
        annotation('textbox',[0.02 0.41 0.05 0.05],'string','c','EdgeColor','none','FontSize',12)
    end
end
for pl = 4:-1:1
    axes('units','centimeters','position', [margin+(pl-1)*(space+sz) margin sz sz/2])
    hstMod(pl).hst = histogram(bpPerPlate(pl).bp.statistics.Qb_values.random_values,'FaceColor',pltMap(pl+1,:));
    ylim([0 1000]);
    xline(bpPerPlate(pl).bp.community.Qb,'r-','linewidth',2);
    pval_mod = sum(bpPerPlate(pl).bp.statistics.Qb_values.random_values>=bpPerPlate(pl).bp.community.Qb)/randomMatsPerPlate;
    pval_antimod = sum(bpPerPlate(pl).bp.statistics.Qb_values.random_values<=bpPerPlate(pl).bp.community.Qb)/randomMatsPerPlate;
    if pval_mod>=1E-4
        tit_mod = sprintf('Modularity, p= %0.2f',pval_mod);
    else
        tit_mod = sprintf('Modularity, p<10^{-4}');
    end
    title(tit_mod);
    if pl==1
        annotation('textbox',[0.02 0.2 0.05 0.05],'string','d','EdgeColor','none','FontSize',12)
    end
end
print([bp_figure_location f 'SuppFigure7'], '-dpng','-r300')

% Save source data
for pl = 4:-1:1
    strctSupp7b.(['Rep' num2str(pl) '_BinEdges']) = hstNest(pl).hst.BinEdges';
    strctSupp7b.(['Rep' num2str(pl) '_Nest']) = hstNest(pl).hst.Values';
    strctSupp7d.(['Rep' num2str(pl) '_BinEdges']) = hstMod(pl).hst.BinEdges';
    strctSupp7d.(['Rep' num2str(pl) '_Mod']) = hstMod(pl).hst.Values';
end
save(['output_files' f 'sourceData' f 'SuppFig7_bd.mat'],'strctSupp7b','strctSupp7d');
toc