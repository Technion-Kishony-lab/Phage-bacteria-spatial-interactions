%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform dN/dS analysis, maintain transition/transvertion ratio %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Supp fig 8 - dN/dS\n');
tic
f = filesep;
simulNum = 10000;
pangBacFile = ['..' f 'input_files' f 'sequencing' f 'bac' f 'genome.mat'];
pangPhgFile = ['..' f 'input_files' f 'sequencing' f 'phg' f 'genome.mat'];
BacGenotypeTable = ['..' f 'Tables' f 'Bac_genotypes_org.xlsx']; 
PhgGenotypeTable = ['..' f 'Tables' f 'Phg_genotypes_org.xlsx']; 
BacIsolatesTable = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsolatesTable = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
phenotypeTable = ['..' f 'Tables' f 'Phenotypes.xlsx']; 

%% Read source data
saveFig = 'Figs';
saveFile = ['output_files' f 'dNdS'];

%% load plate information
phenotypes = table2array(readtable(phenotypeTable)); 

genInfoBac = readtable(BacGenotypeTable,'PreserveVariableNames',true);
genInfoPhg = readtable(PhgGenotypeTable,'PreserveVariableNames',true);

bacCol1 = find(strcmp(genInfoBac.Properties.VariableNames,'genotype_ 1'));
bacColEnd = find(strcmp(genInfoBac.Properties.VariableNames,['genotype_' num2str( size(phenotypes,1))] ));
genotype_bac2 = table2array(genInfoBac(:,bacCol1:bacColEnd));

phgCol1 = find(strcmp(genInfoPhg.Properties.VariableNames,'genotype_  1'));
phgColEnd = find(strcmp(genInfoPhg.Properties.VariableNames,['genotype_' num2str( size(phenotypes,2))] ));
genotype_phg2 = table2array(genInfoPhg(:,phgCol1:phgColEnd)); 

genotype(1).gen = table2array(genInfoBac(:,bacCol1:bacColEnd));
genotype(2).gen = table2array(genInfoPhg(:,phgCol1:phgColEnd)); 
organisms = {'Bacteria','Phage'};

%% Print mutation count:
for org = 1:2
    switch org
        case 1
            genInfo = genInfoBac;
        case 2
            genInfo = genInfoPhg;
    end
    fprintf('%s:\n',organisms{org})
    fprintf('SNPs: %d, nonsynonemous: %d, synonemous: %d, intergenic: %d\n',...
        sum(cellfun(@(x) strcmp(x(1:2),'SN'),genInfo.mutType)),...
        sum(strcmp(genInfo.mutType,'SNP-NonSyn')),...
        sum(strcmp(genInfo.mutType,'SNP-Syn')),...
        sum(strcmp(genInfo.mutType,'SNP-intergenic')));    
    fprintf('CNVs: %d\n',sum(strcmp(genInfo.mutType,'CN')));
    fprintf('InDels: %d\n',sum(cellfun(@(x) (strcmp(x(1:2),'IN') ||...
        strcmp(x(1:2),'MO') || strcmp(x(1:2),'DE')),genInfo.mutType)))
end

%% NS analysis (simulate random SNPs, maintain transition/transversion ratio)
transitionPool = {'AT','GC','TA','CG'};
transversionPool = {'AG','AC','TG','TC','CA','CT','GA','GT'};
if ~exist([saveFile '.mat'],'file')
    dNdS = zeros(simulNum,2);
    for org = 1:2
        if org == 1
            panG = load(pangBacFile);
            genInfo = genInfoBac;
            isolates = readtable(BacIsolatesTable);
        else
            panG = load(pangPhgFile);
            genInfo = genInfoPhg;
            isolates = readtable(PhgIsolatesTable);
        end
        
        indVec = cellfun(@(x) [min(x) max(x)], {panG.pan_g.CDS.indices}, 'UniformOutput', false)';
        ingeneSNPsInd = find((strcmp(genInfo.mutType,'SNP-NonSyn') | strcmp(genInfo.mutType,'SNP-Syn'))...
            & ~genInfo.preexisting_loci);
        ingeneSNPs = genInfo.Call(ingeneSNPsInd);
        gen = genInfo(ingeneSNPsInd,find(strcmp(genInfo.Properties.VariableNames,'genotype_ 1')):end);     
        transitions = ismember(ingeneSNPs,transitionPool);
        transversions = ismember(ingeneSNPs,transversionPool);
        % find number of replicates:
        transitionsInd = find(transitions);
        tranversionsInd = find(transversions);
        platesTransitions = zeros(sum(transitions),4);
        platesTransversions = zeros(sum(transversions),4);
        for snp = 1:numel(transitionsInd)
            isos = find(table2array(gen(transitionsInd(snp),:)));
            platesTransitions(transitionsInd(snp),unique(isolates.replicateNum(isos))) = 1;
        end
        for snp = 1:numel(tranversionsInd)
            isos = find(table2array(gen(tranversionsInd(snp),:)));
            platesTransversions(tranversionsInd(snp),unique(isolates.replicateNum(isos))) = 1;
        end
        transitionsCnt = sum(sum(platesTransitions));
        transversionsCnt = sum(sum(platesTransversions));
        intraSeq = cellfun(@(z) panG.pan_g.Sequence(z(1):z(2)),indVec,'UniformOutput',false);
        allLngth = cellfun(@(x) x(2)-x(1)+1, indVec);
        for g = 1:numel(panG.pan_g.CDS)
            if length(panG.pan_g.CDS(g).indices)==4
                sortedIndices = sort(panG.pan_g.CDS(g).indices);
                intraSeq{g} = [panG.pan_g.Sequence(sortedIndices(1): sortedIndices(2)),...
                   panG.pan_g.Sequence(sortedIndices(3):sortedIndices(4))];
                allLngth(g) = sortedIndices(4)-sortedIndices(3)+...
                    sortedIndices(2)-sortedIndices(1)+2;
            elseif length(panG.pan_g.CDS(g).indices)==6
                sortedIndices = sort(panG.pan_g.CDS(g).indices);
                intraSeq{g} = [panG.pan_g.Sequence(sortedIndices(1): sortedIndices(2)),...
                   panG.pan_g.Sequence(sortedIndices(3):sortedIndices(4)),...
                   panG.pan_g.Sequence(sortedIndices(5):sortedIndices(6))];
                allLngth(g) = sortedIndices(6)-sortedIndices(5)+...
                    sortedIndices(4)-sortedIndices(3)+...
                    sortedIndices(2)-sortedIndices(1)+3;
            end
        end  
        syn = zeros(simulNum,1);
        nonsyn = zeros(simulNum,1);
        for sim = 1:simulNum
            if mod(sim,100)==0
                fprintf('%d\n',sim)
            end
            syn(sim) = 0;
            nonsyn(sim) = 0;
            for t = 1:transitionsCnt+transversionsCnt
                pos = randi(sum(allLngth));
                gene = find(pos<=cumsum(allLngth),1);
                ingenepos = pos - sum(allLngth(1:gene-1));
                nt_pos = mod(ingenepos,3);
                if nt_pos == 0
                    nt_pos = 3;
                end
                refCodon = upper(strcat(intraSeq{gene}(ingenepos-nt_pos+1),intraSeq{gene}(ingenepos-nt_pos+2),...
                    intraSeq{gene}(ingenepos-nt_pos+3)));
                if strcmp(panG.pan_g.CDS(gene).location(1:4),'comp')  
                    refCodon = seqrcomplement(refCodon);
                end 
                altcodon = refCodon;
                if t <= transitionsCnt
                    [~,from] = ismember(refCodon(nt_pos),cellfun(@(x) x(1),transitionPool));
                    altcodon(nt_pos) = transitionPool{from}(2);
                else
                    [~,from] = ismember(refCodon(nt_pos),cellfun(@(x) x(1),transversionPool));
                    to = rand(1);
                    if to<0.5
                        altcodon(nt_pos) = transversionPool{from}(2);
                    else
                        altcodon(nt_pos) = transversionPool{from+1}(2);
                    end

                end
                if strcmp(nt2aa(refCodon,'AlternativeStartCodons', false),nt2aa(altcodon,'AlternativeStartCodons', false))
                    syn(sim) = syn(sim)+1;
                else
                    nonsyn(sim) = nonsyn(sim)+1;
                end

            end
            dNdS(sim,org) = nonsyn(sim)./syn(sim);
        end
        clear panG
    end
    save(saveFile,'dNdS'); 
else
    load(saveFile,'dNdS')     
end

%% Plot:
figure(103); clf;
set(gcf,'Units','centimeters','Position',[0.5 0.5 20 9])
realdNdS = zeros(1,2);
clr = [0.5 0.5 0.5];
nbins = [30,12];
xlims = [0 13;
    0 8];
for org = 1:2
    axes('Units','centimeters','Position',[2+(org-1)*9 1 7 7])
    switch org
        case 1
            realdNdS(org) = sum(strcmp(genInfoBac.mutType,'SNP-NonSyn') & ~genInfoBac.preexisting_loci)./...
                sum(strcmp(genInfoBac.mutType,'SNP-Syn') & ~genInfoBac.preexisting_loci);
        case 2
            realdNdS(org) = sum(strcmp(genInfoPhg.mutType,'SNP-NonSyn') & ~genInfoPhg.preexisting_loci)./...
                sum(strcmp(genInfoPhg.mutType,'SNP-Syn') & ~genInfoPhg.preexisting_loci);
    end
    histogram(dNdS(:,org),nbins(org),'FaceColor',clr); hold on
    xline(realdNdS(org),'color','red','linewidth',2)
    ylim([0 3500])
    xlim(xlims(org,:))
    legend({'Simulations','Data'})
    xlabel('No. Nonsynonymous / No. synonymous SNPs')
    ylabel('No. simulations')
    tit = sprintf('%s, p = %0.4f',organisms{org},numel(find(dNdS(:,org)>=realdNdS(org)))/simulNum);
    title(tit)
end
print([saveFig f 'SuppFigure8'],'-dpng','-r300')
toc
