%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Count the weighted number of mutation events per gene,%%%
%%% compared to a random mutation simulation              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Multi-mutated genes \n');
tic
f = filesep;
organisms = {'bacteria','phage'};
figure_location = 'Figs';
refs_folder = ['..' f 'input_files' f 'sequencing'];
PhgGeneList = ['..' f 'input_files' f 'PathwayTables' f 'phage_genes.xlsx'];
pathwayFolder =  ['..' f 'input_files' f 'PathwayTables'];
pwClrMapFile = ['output_files' f 'pathwayClrMap.mat'];
multiMut = [];
%% Params:
wghtFactor = 1000;
binSizes = [1,1];
simulNums = [500,5000];
letter_fs = 14;
sigP = 0.05; % P-val for significance
%% source tables
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
BacGenotypeTable = ['..' f 'Tables' f 'Bac_genotypes_org.xlsx']; 
PhgGenotypeTable = ['..' f 'Tables' f 'Phg_genotypes_org.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 

%% load plate information
phenotypes = table2array(readtable(phenotypeTable)); 

isolatesTable(1).isolates = readtable(BacIsolatesTable);
isolatesTable(2).isolates = readtable(PhgIsolatesTable);

mutationInfo(1).mutinfo = readtable(BacGenotypeTable,'PreserveVariableNames',true);
bacCol1 = find(strcmp(mutationInfo(1).mutinfo.Properties.VariableNames,'genotype_ 1'));
bacColEnd = find(strcmp(mutationInfo(1).mutinfo.Properties.VariableNames,['genotype_' num2str( size(phenotypes,1))] ));
genotype(1).gen = table2array(mutationInfo(1).mutinfo(:,bacCol1:bacColEnd));

mutationInfo(2).mutinfo = readtable(PhgGenotypeTable,'PreserveVariableNames',true);
phgCol1 = find(strcmp(mutationInfo(2).mutinfo.Properties.VariableNames,'genotype_  1'));
phgColEnd = find(strcmp(mutationInfo(2).mutinfo.Properties.VariableNames,['genotype_' num2str( size(phenotypes,2))] ));
genotype(2).gen = table2array(mutationInfo(2).mutinfo(:,phgCol1:phgColEnd));

%% load pathway colomap from regression
if exist(pwClrMapFile,'file')
    load(pwClrMapFile)
end
%%
for org = 2:-1:1
    organism = organisms{org};
    binSize = binSizes(org);
    simulNum = simulNums(org);
    nms = isolatesTable(org).isolates.SequenceNms;
    simulFile = ['output_files' f 'simulMultiMut_' organism '.mat'];  
    switch organism
        case 'bacteria'
            ref_folder = [refs_folder f 'bac' f 'genome.mat'];
        case 'phage'
            ref_folder = [refs_folder f 'phg' f 'genome.mat'];
            phage_genes = readtable(PhgGeneList);
    end
    load(ref_folder,'pan_g');
    indVec = [cellfun(@(x) min(x), {pan_g.CDS.indices}); cellfun(@(x) max(x),...
        {pan_g.CDS.indices}) ]'; % gene min and max indices

    %%
    geneLen = zeros(numel(pan_g.CDS),1);
    for gn = numel(pan_g.CDS):-1:1
        geneLen(gn) = max(pan_g.CDS(gn).indices)-min(pan_g.CDS(gn).indices);
    end
    if org == 2  % compensate for long edges
        geneLen(1) = max(pan_g.CDS(1).indices);
        geneLen(numel(pan_g.CDS)) = size(pan_g.Sequence,2)-min(pan_g.CDS(numel(pan_g.CDS)).indices);
    end
    meanGeneLen = mean(geneLen);
    gen = genotype(org).gen;
    mutinfo = mutationInfo(org).mutinfo;
    mutGenes = struct('genes',[],'wgts',[],'mutPltIdx',[]);
    plates = zeros(size(gen,1),1);
    wgts = wghtFactor./...
        max(mutinfo.genomePos2-mutinfo.genomePos1,wghtFactor);

    for mutIdx = size(gen,1):-1:1
        isos = find(gen(mutIdx,:));
        if ~isempty(isos) && ~mutinfo.preexisting_loci(mutIdx)==1
            plates(mutIdx) = numel(unique(isolatesTable(org).isolates.replicateNum(isos)));
            mutGenes(mutIdx).mutPltIdx = repmat(mutIdx,plates(mutIdx),1);
        else
            plates(mutIdx) = 0;
        end
       
        % Find first and last mutated genes:        
        geneNums = assignMutGene(pan_g,[mutinfo.genomePos1(mutIdx) ... 
            mutinfo.genomePos2(mutIdx)],mutinfo.mutType{mutIdx}); 
        g1 = geneNums.gene1; 
        g2 = geneNums.gene2; 
        
        mutGenes(mutIdx).genes = (g1:g2)';
        mutGenes(mutIdx).wgts = (repmat(wgts(mutIdx),g2-g1+1,1))*plates(mutIdx);
    end
    allWgts = cat(1,mutGenes.wgts);
    allGns =  cat(1,mutGenes.genes);
    allMutPltIdx = cat(1,mutGenes.mutPltIdx);
    uGns = unique(allGns);
    geneName = [];
    wgtPrGn = zeros(numel(uGns),1);
    for gn = numel(uGns):-1:1
        if org ==2
            geneName{gn} = erase(erase(erase((phage_genes.genes{uGns(gn)}),'gene '),' possible'),'/helicase [14');
            geneSizeFact = meanGeneLen/geneLen(uGns(gn));
            wgtPrGn(gn) = sum(allWgts(allGns==uGns(gn)));
        else
            geneName{gn} = pan_g.CDS(uGns(gn)).gene;
            geneSizeFact = meanGeneLen/sum(geneLen(strcmp(geneName{gn}, {pan_g.CDS.gene}))); % sum because of yibX gene #85
            wgtPrGn(gn) = sum(allWgts(allGns==uGns(gn))); 
        end
    end 
    geneCount = ceil(max(wgtPrGn));
    binEdges = 0:binSize:ceil(max(wgtPrGn));
    geneCountHistWeighted = histcounts(wgtPrGn,binEdges);

    %% Simulate random mutation events on the genome simulNum times
    if ~exist(simulFile,'file')
        fprintf('simulating mutations...\n');
        Values = zeros(simulNum,numel(binEdges)-1);
        for sim = 1:simulNum
            geneCountSimulWeighted = zeros(size(pan_g.CDS,2),1);
            if mod(sim,50)==0
                fprintf('%d\n',sim)
            end
            [allMutgenesSimul,allSimulWeights] = simulateMutations(mutinfo,...
                org,pan_g,allMutPltIdx,wghtFactor);   
            [~,~,iw] = unique(allMutgenesSimul,'stable');  
            for i = max(iw):-1:1
                geneCountSimulWeighted(i) = sum(allSimulWeights(iw==i));
            end
            
            Values(sim,:) = histcounts(geneCountSimulWeighted,binEdges); 
        end
        save(simulFile,'Values');
    else
        fprintf('loading simulation results...');
        load(simulFile,'Values');
        fprintf('done.\n');
    end

    %% find limit of significance 
    lngPltCnt = max(geneCount);
    geneCntSum = zeros(size(geneCountHistWeighted));
    ValuesSum = zeros(simulNum,size(geneCountHistWeighted,2));
    ValuesMean = mean(Values,1);
    ValuesErr = std(Values,1);
    ValuesSumMean = zeros(size(Values,2),1);
    for i = 1:numel(geneCountHistWeighted)
        geneCntSum(i) = sum(geneCountHistWeighted(i:end));
        ValuesSum(:,i) = sum(Values(:,i:end),2);
        ValuesSumMean(i) = mean(ValuesSum(:,i));
    end
    ValuesSumCI = prctile(ValuesSum,[2.5,97.5]); % confidence intervals for error bars
    sigWght = (find(sum(ValuesSum >= geneCntSum)/simulNum<sigP,1)-1)*binSize;
    sigGnInd = find(wgtPrGn>=sigWght);
    sigGenes = geneName(sigGnInd);

    multiMut(org).ValuesSum = ValuesSum;
    multiMut(org).ValuesSumMean = ValuesSumMean;
    multiMut(org).ValuesSumCI = ValuesSumCI;
    multiMut(org).sigWght = sigWght;
    multiMut(org).sigGnInd = sigGnInd;
    multiMut(org).sigGenes = sigGenes;
    multiMut(org).wgtPrGn = wgtPrGn;   
    multiMut(org).sigGnCnt = wgtPrGn(sigGnInd);
    multiMut(org).lngPltCnt = lngPltCnt;
    multiMut(org).geneCntSum = geneCntSum;
end    

%% Plot commulative  Weight per gene
% Figure params:
sz = 6.5;
margin_x = 1.5;
margin_y = 1.5;
space_x = 1.5;
space_y = 1.5;
gn_fs = 10;
pw_fs = 11;
simBarClr = [0.7 0.7 0.7];
orgTit = {'Bacteria','Phages'};

figure(102); clf;
set(gcf,'units','centimeters','position',[0.5 0.5 margin_x*2+space_x+2*sz margin_y*2+space_y+2*sz]);
for org = 1:2
    ax_weights = axes('units','centimeters','position',[margin_x+(space_x+sz)*(org-1)...
        margin_y+space_y+sz sz sz]);
    bar(ax_weights,binSizes(org):binSizes(org):(multiMut(org).lngPltCnt-binSizes(org)),...
        multiMut(org).geneCntSum(2:end),'FaceColor',simBarClr)
    hold on
    errorbar(ax_weights,binSizes(org):binSizes(org):(multiMut(org).lngPltCnt-binSizes(org)),...
        multiMut(org).ValuesSumMean(2:end),...
        multiMut(org).ValuesSumMean(2:end)-reshape(multiMut(org).ValuesSumCI(1,2:end),[],1),...
        reshape(multiMut(org).ValuesSumCI(2,2:end),[],1)-multiMut(org).ValuesSumMean(2:end),...
        '.','MarkerSize',12,'color','red')
    switch org
        case 1
            set(gca,'xtick',[1,3,5,7],'xticklabels',{'\geq1','\geq3','\geq5','\geq7'})
            xlim1 = 8.5;
            ylim1 = 90;
            ylim2 = max(multiMut(org).sigGnCnt)+1;
            ltrX = [0.025,0.025];
            ltrY = [0.94,0.48];
            letters = {'a','c'};
        case 2
            set(gca,'xtick',[1,10,20],'xticklabels',{'\geq1','\geq10','\geq20'})
            xlim1 = 26;
            ylim1 = 50;
            ylim2 = max(multiMut(org).sigGnCnt)+1;
            ltrX = [0.475,0.475];
            ltrY = [0.94,0.48];
            letters = {'b','d'};
    end
    hold on
    xl = xline(multiMut(org).sigWght,'k--');
    l = legend({'Data','Simulation'});
    ylim([0 ylim1]);
    xlim([0 xlim1]);
    title(orgTit{org});
    xlabel('Weighted number of mutations per gene')
    ylabel('Number of genes')
    annotation('textbox',[ltrX(1),ltrY(1),0.01,0.01],'String',letters{1},'EdgeColor','none','FontSize',letter_fs);

    % Bar plot of all significant genes:
    ax_bar = axes('units','centimeters','position',[margin_x+(space_x+sz)*(org-1) margin_y sz sz]);
    [sortCnt,sortI] = sort(multiMut(org).sigGnCnt,'descend');
    switch org
        case 1
            pwFile = [pathwayFolder f 'PathwayBacOrg.xlsx'];
        case 2 
            pwFile = [pathwayFolder f 'PathwayPhgOrg.xlsx'];
    end
    pwClrMap = zeros(numel(multiMut(org).sigGenes),3);
    pwStrings = [];
    if exist(pwFile,'file') && exist(pwClrMapFile,'file')
        load(pwClrMapFile);
        allPW = readtable(pwFile);
        for gn = numel(multiMut(org).sigGenes(sortI)):-1:1
            pws = allPW.pathway(cellfun(@(x) ...
                contains(x,[multiMut(org).sigGenes{sortI(gn)} ':']),allPW.mutation));
            if numel(pws)>1
                pw = pws{find(~strcmp(pws,'multiple'),1)};
            elseif ~isempty(pws)
                pw = pws{1};
            else
                pw = [];
            end
            switch org
                case 1
                    if ~isempty(pw) && any(strcmp(pw,bac_pw)) && ~strcmp(pw,'multiple')
                        pwStrings{gn} = pw;
                        pwClrMap(gn,:) = bacPwsClrMap(strcmp(pw,bac_pw),:);
                    elseif contains(sigGenes{sortI(gn)},'waa') % because waaA is not wrriten explicitly in the pathway list
                        pwStrings{gn} = 'LPS';
                        pwClrMap(gn,:) = bacPwsClrMap(strcmp(pwStrings{gn},bac_pw),:);
                    end
                case 2
                    if ~isempty(pw) && any(strcmp(pw,phg_pw)) && ~strcmp(pw,'multiple')
                        pwStrings{gn} = pw;
                        pwClrMap(gn,:) = phgPwsClrMap(strcmp(pw,phg_pw),:);
                    end            
            end
        end
    end
    for br = 1:numel(sortCnt)
        barClr = pwClrMap(br,:);
        bar(br,sortCnt(br),'facecolor',barClr,'Parent',ax_bar);    
        hold on
    end
    upwStrings = unique(pwStrings(cellfun(@(x) ~isempty(x),pwStrings)));
    if org==1
        for pw = 1:numel(upwStrings)
        text(numel(sortCnt)-5,sortCnt(1)+1-pw,upwStrings{pw},'Color',...
            bacPwsClrMap(strcmp(upwStrings{pw},bac_pw),:),'FontSize',pw_fs)
        end
    else
        for pw = 1:numel(upwStrings)
        text(numel(sortCnt),sortCnt(1)+1-4*pw,upwStrings{pw},'Color',...
            phgPwsClrMap(strcmp(upwStrings{pw},phg_pw),:),'FontSize',pw_fs)
        end
    end

    if org==1
        set(gca,'Xtick',1:numel(sortCnt),'XTickLabel',multiMut(org).sigGenes(sortI),'XTickLabelRotation',60,'FontSize',gn_fs)
    else
        set(gca,'Xtick',1:numel(sortCnt),'XTickLabel',strcat('gp',multiMut(org).sigGenes(sortI)),'XTickLabelRotation',60,'FontSize',gn_fs)
    end
    ylabel('Weighted number of mutations')
    ylim([0 ylim2]);
    annotation('textbox',[ltrX(2),ltrY(2),0.01,0.01],'String',letters{2},'EdgeColor','none','FontSize',letter_fs);

    print([figure_location f 'SuppFigure9'],'-dpng','-r300'); 
end
toc