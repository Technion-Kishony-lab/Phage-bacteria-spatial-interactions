%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main phage-bacteria co-evolution data analysis pipeline  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stand in: 'D:\Research\OneDrive - Technion\Ecoli-Phage-CoEvolution\scripts'
fprintf('Running main analysis\n');
tic
f = filesep;
addpath(genpath('.'))
if ~exist('output_files','dir')
    mkdir('output_files')
end
if ~exist('Figs','dir')
    mkdir('Figs')
end
if ~exist(['..' f 'Tables'],'dir')
    mkdir(['..' f 'Tables'])
end
if ~exist(['..' f 'script_data'],'dir')
    mkdir(['..' f 'script_data'])
end
%% Set definitions and parameters
prefix_Phage = 'Sample_'; % prefix for phage sample names
prefix_Bac = 'Sample_Bac_'; % prefix for bacterial sample names
WT_bac = 'Sample_Bac_MG-Y_1'; % bacterial WT sequencing sample name
WT_phg = 'Sample_Phage_T7_WT_1'; % phage WT sequencing sample name
bacRef = 'U00096.3';  % NCBI numebr of reference bacteria
phgRef = 'NC_001604.1'; % NCBI numebr of reference phage
polymorph = 0.7; % minimal percentage for a chosen allel in the case of polymorphism
minCov = [5,15]; % minimal total coverage in the case of polymorphism. below that place WT allele. [bac,phg];
hetroCutoff = 1; % number of 'n's that make an isolate a "heterozygote" and removed from analysis
linkage_method = 'complete'; % mathod for linkage based on phenotype pairwise distances
plateSize = 1000; % diameter of each plate [pxl] in the full size images, initial and continual coevolution 
plateSizeSP = 1100; % diameter of each plate [pxl] in the full size images, no-phage control
pltNms = {'UL','UR','DL','DR'}; % plate replicate names (Up-Left, Up-Right etc.)
fixPltTh = 0; % set to 1 to manually adjust plates saturation thresholds in the infection assay 
fixWellTh = 0; % set to 1 to manually adjust single wells saturation thresholds in the infection assay 
%% Input locations of analyzed sequencing data
if exist(['..' f 'input_files'],'dir')
    refBacLoc = ['..' f 'input_files' f 'sequencing' f 'bac' f 'genome.mat'];
    refPhgLoc = ['..' f 'input_files' f 'sequencing' f 'phg' f 'genome.mat'];
    seq_DPs_bac = ['..' f 'input_files' f 'sequencing' f 'bac' f 'DP'];
    Sequencing_analysis_output = ['..' f 'input_files' f 'sequencing' f 'Sequencing_analysis_results']; 
    seq_ana_bac = ['..' f 'input_files' f 'sequencing' f 'bac' f 'Ana_Bacteria_1'];
    seq_ana_phg = ['..' f 'input_files' f 'sequencing' f 'phg' f 'Ana_Phage'];
    breseqResBac = ['..' f 'input_files' f 'sequencing' f 'bac' f 'Bac_isolates_output'];
    breseqResPhg = ['..' f 'input_files' f 'sequencing' f 'phg' f 'Phg_isolates_output'];
    phage_preexisting_loci_File = ['..' f 'input_files' f 'sequencing' f 'phage_preexisting_loci.xlsx']; 
    not_clean_phg_xls = ['..' f 'input_files' f 'crossInfection' f 'check_phg_uniformity.xlsx'];
else
    refBacLoc = ['..' f 'script_data' f 'sequencing' f 'bac' f 'genome.mat'];
    refPhgLoc = ['..' f 'script_data' f 'sequencing' f 'phg' f 'genome.mat'];
    seq_DPs_bac = ['..' f 'script_data' f 'sequencing' f 'bac' f 'DP'];
    Sequencing_analysis_output = ['..' f 'script_data' f 'sequencing' f 'Sequencing_analysis_results']; 
    seq_ana_bac = ['..' f 'script_data' f 'sequencing' f 'bac' f 'Ana_Bacteria_1'];
    seq_ana_phg = ['..' f 'script_data' f 'sequencing' f 'phg' f 'Ana_Phage'];
    breseqResBac = ['..' f 'script_data' f 'sequencing' f 'bac' f 'Bac_isolates_output'];
    breseqResPhg = ['..' f 'script_data' f 'sequencing' f 'phg' f 'Phg_isolates_output'];
    phage_preexisting_loci_File = ['..' f 'script_data' f 'sequencing' f 'phage_preexisting_loci.xlsx']; 
    not_clean_phg_xls = ['..' f 'script_data' f 'crossInfection' f 'check_phg_uniformity.xlsx'];   
end

%% Coevolution raw images and analysis files (sampling coordinates, plate mask files)
raw_imgs_main = ['..' f 'input_files' f 'coevolution' f '20191206_plate_raw_imgs']; % coevolution images - first experiment (Fig 1)
raw_imgs_cont = ['..' f 'input_files' f 'coevolution' f '20191206_plate_raw_imgs' f '20191218_continue']; % coevolution images - continual experiment (Supp Fig 1)
raw_imgs_noPhg = ['..' f 'input_files' f 'coevolution' f '220904_MGY_SP'];
main_crop_file = ['..' f 'input_files' f 'coevolution' f 'crop_main.mat']; % contains the left upper corner of the left upper plate and the bottom right corner of the bottom right plate
cont_crop_file = ['..' f 'input_files' f 'coevolution' f 'crop_cont.mat']; % contains the left upper corner of the left upper plate and the bottom right corner of the bottom right plate - continual experiment
noPhg_crop_file = ['..' f 'input_files' f 'coevolution' f 'crop_SP.mat'];
masks_main = ['..' f 'input_files' f 'coevolution' f 'masks_main.mat']; % masks for plate locations in the coevolution images.
masks_cont = ['..' f 'input_files' f 'coevolution' f 'masks_cont.mat']; % masks for plate locations in the continual coevolution images.
noPhg_mask_file = ['..' f 'input_files' f 'coevolution' f 'masks_SP.mat'];

SamplingInfo_file = ['..' f 'input_files' f 'coevolution' f 'SamplingInfo.xlsx']; % Coordinates of all coevolution samples
%% Cross infection input files:
% Raw images of infection assay (total 2 experiments):
infection_imgs_p1 = ['..' f 'input_files' f 'crossInfection' f 'coevo_cross_infection' f '20200205_phage1_interactions' f];
infection_imgs_p2 = ['..' f 'input_files' f 'crossInfection' f 'coevo_cross_infection' f '20200202_phage2_interactions' f];
infection_name_order_file = ['..' f 'input_files' f 'crossInfection' f 'bac+phage_name_order_infection.xlsx'];
% Saved adjusted thresholds for plaque detection:
well_thresh_xls1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'single_well_threshold_20211117.xlsx'];
well_thresh_xls2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'single_well_threshold_20211117.xlsx'];
plate_thresh_xls1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'Adjusted_threshold.xls'];
plate_thresh_xls2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'Adjusted_threshold.xls'];
plateWT1_thresh_xls = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'WT_threshold.xls'];
plateWT2_thresh_xls = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'WT_threshold.xls'];
% map of inoculated phages:
inoc_plate_file1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'dilution_Ph1_20200205.xlsx'];
inoc_plate_file2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'dilution_phage2_20200202.xlsx'];
% images of plate lids (contain wells that should be removed due to merged drops)
lids_folder_p1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'lids_merged'];
lids_folder_p2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'lids_merged'];
remove_merged_1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1' f 'remove_merged_Ph1.xls'];
remove_merged_2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2' f 'remove_merged_Ph2.xls'];
% mat files with adjusted infection images will be saved here:
save_folder_p1 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp1'];
save_folder_p2 = ['..' f 'input_files' f 'crossInfection' f 'infectionExp2'];
%% save locations
% Coevolution imaging:
sngPltMFile = ['..' f 'script_data' f 'sngPltM.mat']; % Initial coevolution
sngPltCFile = ['..' f 'script_data' f 'sngPltC.mat']; % Continual coevolution
sngPltSPFile = ['..' f 'script_data' f 'sngPltSP.mat']; % No-Phage Control 

sngPltMeanMFile = ['..' f 'script_data' f 'sngPltMeanM.mat'];
sngPltMeanCFile = ['..' f 'script_data' f 'sngPltMeanC.mat'];
sngPltMeanSPFile = ['..' f 'script_data' f 'sngPltMeanSP.mat'];
movieParamsFile = ['..' f 'script_data' f 'movieParams.mat'];
% Cross-Infection saving folder
infectStructFile1 = ['..' f 'script_data' f 'infectStruct1.mat'];
infectStructFile2 = ['..' f 'script_data' f 'infectStruct2.mat'];
% Folder for cropped and reduced images for time lapse movies
out_imgs_main = ['..' f 'script_data' f 'imgs_main']; % folder for movie images - main experiment
out_imgs_cont = ['..' f 'script_data' f 'imgs_cont']; % folder for movie images - continual experiment
out_imgs_noPhg  = ['..' f 'script_data' f 'imgs_SP'];
%% Genetic pathways of mutated gene lists:
if exist(['..' f 'input_files'],'dir')
    PathwayBacUni = ['..' f 'input_files' f 'PathwayTables' f 'PathwayBacUni.xlsx'];
    PathwayBacOrg = ['..' f 'input_files' f 'PathwayTables' f 'PathwayBacOrg.xlsx'];
    PathwayPhgUni = ['..' f 'input_files' f 'PathwayTables' f 'PathwayPhgUni.xlsx'];
    PathwayPhgOrg = ['..' f 'input_files' f 'PathwayTables' f 'PathwayPhgOrg.xlsx'];
    phage_readable_genes = ['..' f 'input_files' f 'PathwayTables' f 'phage_genes.xlsx'];
else
    PathwayBacUni = ['..' f 'script_data' f 'PathwayTables' f 'PathwayBacUni.xlsx'];
    PathwayBacOrg = ['..' f 'script_data' f 'PathwayTables' f 'PathwayBacOrg.xlsx'];
    PathwayPhgUni = ['..' f 'script_data' f 'PathwayTables' f 'PathwayPhgUni.xlsx'];
    PathwayPhgOrg = ['..' f 'script_data' f 'PathwayTables' f 'PathwayPhgOrg.xlsx'];
    phage_readable_genes = [];
end
%% Tables saving locations:
BacIsoTblFile = ['..' f 'Tables' f 'Bac_isolates.xlsx']; 
PhgIsoTblFile = ['..' f 'Tables' f 'Phg_isolates.xlsx']; 
fullPhenTblFile = ['..' f 'Tables' f 'FullPhenotypes.xlsx']; 
phenTblFile = ['..' f 'Tables' f 'Phenotypes.xlsx']; 
turbidTblFile = ['..' f 'Tables' f 'Turbidity.xlsx']; 
turbidTblOrgOrderFile = ['..' f 'Tables' f 'TurbidityOrgOrder.xlsx']; 
BacGenTblOrgFile = ['..' f 'Tables' f 'Bac_genotypes_org.xlsx'];
BacGenTblUniFile = ['..' f 'Tables' f 'Bac_genotypes_uni.xlsx'];
PhgGenTblOrgFile = ['..' f 'Tables' f 'Phg_genotypes_org.xlsx'];
PhgGenTblUniFile = ['..' f 'Tables' f 'Phg_genotypes_uni.xlsx'];


%% Start Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Time lapse movies %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proccess coevolution images into a movie, save cropped images averaged over small areas for peak and dominance analysis
chnl = 1; % BW images are taken from the Red channel. 
xyWin = 10; % average pixel intensity on this window size
tWin = 4; %  intensity is averaged on this sliding time window (tWin X 10min)
A = [130, 130, 10000]; % Logarithmic frame rate: parameter in log function [initial, continual, no-phage]
N = [300, 340, 201]; % number of images in movie [initial, continual, no-phage]
ti = [1, 20, 1]; % initial image index in movie [initial, continual, no-phage]
MinPeakDistance = 20; % for findPeaks 
MinPeakProminence = 4; % for findPeaks
save(movieParamsFile,'N','A','ti','chnl','xyWin','tWin','plateSize','plateSizeSP');
if ~exist(sngPltMeanMFile,'file')
    [sngPltM,sngPltMeanM] = buildImgsStruct(raw_imgs_main,main_crop_file,masks_main,sngPltMFile,...
        sngPltMeanMFile,plateSize,chnl,xyWin,tWin,MinPeakDistance,MinPeakProminence);
end
if ~exist(sngPltMeanCFile,'file')
    [sngPltC,sngPltMeanC] = buildImgsStruct(raw_imgs_cont,cont_crop_file,masks_cont,sngPltCFile,...
        sngPltMeanCFile,plateSize,chnl,xyWin,tWin,MinPeakDistance,MinPeakProminence);
end
% create images for time-lapse movies
if ~exist(out_imgs_main,'dir')
    mkdir(out_imgs_main)
end
bar_pos_main = [plateSize-500 (plateSize*2)-20]; % position of time bar in movies
bgCntr(1,:) = [449,487,480,517];
bgCntr(2,:) = [474,503,1504,1540];
bgCntr(3,:) = [1454,1495,459,497];
bgCntr(4,:) = [1495,1530,1472,1508];
SuppMovieImages(raw_imgs_main,main_crop_file,out_imgs_main,plateSize,A(1),...
    N(1),ti(1),bar_pos_main,bgCntr,chnl,'main')
if ~exist(out_imgs_cont,'dir')
    mkdir(out_imgs_cont)
end
bar_pos_cont = [plateSize-560 (plateSize*2)-20];
SuppMovieImages(raw_imgs_cont,cont_crop_file,out_imgs_cont,plateSize,A(2),...
    N(2),ti(2),bar_pos_cont,bgCntr,chnl,'pre')  %% CHANGE PRE TO CONT

% Save timelapse images and analyze peaks in the no-phage control experiment
buildImgsStruct_SP(raw_imgs_noPhg,noPhg_crop_file,noPhg_mask_file,sngPltSPFile,...
        sngPltMeanSPFile,movieParamsFile,out_imgs_noPhg,MinPeakDistance,MinPeakProminence);
%% Read NCBI reference genome files
fprintf('Reading genbank reference files...\n');
gbB = load(refBacLoc);
refGB_bac = gbB.pan_g;
gbP = load(refPhgLoc);
refGB_phg = gbP.pan_g;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Cross-infection %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and analyze cross-infection data
if ~exist(infectStructFile1,'file') || ~exist(infectStructFile2,'file') || fixPltTh || fixWellTh 
    WT_pos = [24,86]; % Wells with the wildtype T7 phage
    wellSize = 90; % 1/2 well length/width in pxl
    min_plq = 30; % minimal total pixel number in a region to be considered as a plaque (and not noise).
    % Cross-infection Exp. 1
    if ~exist(infectStructFile1,'file') 
        fprintf('Building infection matrix data 1...\n');
        infectStruct1 = buildInfectionMat(infection_imgs_p1,save_folder_p1,...
            inoc_plate_file1,infection_name_order_file,WT_pos,lids_folder_p1,remove_merged_1,...
            plate_thresh_xls1,plateWT1_thresh_xls,well_thresh_xls1,...
            wellSize,min_plq,1,infectStructFile1,fixPltTh,fixWellTh);
        save(infectStructFile1,'infectStruct1','-v7.3');
    else
        fprintf('Loading infection matrix data 1...\n');
        load(infectStructFile1,'infectStruct1')
    end
    % Cross-infection Exp. 2
    if ~exist(infectStructFile2,'file') 
        fprintf('Building infection matrix data 2...\n');
        infectStruct2 = buildInfectionMat(infection_imgs_p2,save_folder_p2,...
            inoc_plate_file2,infection_name_order_file,WT_pos,lids_folder_p2,remove_merged_2,...
            plate_thresh_xls2,plateWT2_thresh_xls,well_thresh_xls2,...
            wellSize,min_plq,2,infectStructFile2,fixPltTh,fixWellTh);
        save(infectStructFile2,'infectStruct2','-v7.3');
    else
        fprintf('Loading infection matrix data 2...\n');
        load(infectStructFile2,'infectStruct2')
    end

    clear infection_imgs_p1 infection_imgs_p2 inoc_plate_file1 inoc_plate_file2 ...
        lids_folder_p1 lids_folder_p2 save_folder_p1 save_folder_p2 remove_merged_1 remove_merged_2 
else
    fprintf('Loading infection matrices...\n')
    load(infectStructFile1,'infectStruct1');
    load(infectStructFile2,'infectStruct2');

end
infectivity1 = [infectStruct1.plaques.infectivity];
isMerged1 = [infectStruct1.plaques.takeWell];
infectivity1(isnan(isMerged1)) = NaN;
infectivity2 = [infectStruct2.plaques.infectivity];   
isMerged2 = [infectStruct2.plaques.takeWell];
infectivity2(isnan(isMerged2)) = NaN;

fullInfectionMat = [infectivity1; infectivity2]';
fullInfectionTble = table(fullInfectionMat);
writetable(fullInfectionTble,fullPhenTblFile)
fullInfectionMat = table2array(readtable(fullPhenTblFile));

% Same for turbidity:
turbidity1 = [infectStruct1.plaques.turbidity];
isMerged1 = [infectStruct1.plaques.takeWell];
turbidity1(isnan(isMerged1)) = NaN;
turbidity2 = [infectStruct2.plaques.turbidity];   
isMerged2 = [infectStruct2.plaques.takeWell];
turbidity2(isnan(isMerged2)) = NaN;
fullTurbidityMat = [turbidity1; turbidity2]';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Mutations %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save all isolate names according to their sequencing name and order:
if ~exist(['..' f 'Tables' f 'seqNmsBac.xlsx'],'file')
    seqNmsBacTmp = load([seq_ana_bac f 'AnaIsolatesNames.mat']);
    seqNmsBac = table(seqNmsBacTmp.AnaIsolatesNames.RawSubDirName,'VariableNames',{'Nms'});
    writetable(seqNmsBac,['..' f 'Tables' f 'seqNmsBac.xlsx']);
 else
    seqNmsBac = readtable(['..' f 'Tables' f 'seqNmsBac.xlsx']);   
end 
if ~exist(['..' f 'Tables' f 'seqNmsPhg.xlsx'],'file')
    seqNmsPhgTmp = load([seq_ana_phg f 'AnaIsolatesNames.mat']);
    seqNmsPhg = table(seqNmsPhgTmp.AnaIsolatesNames.RawSubDirName,'VariableNames',{'Nms'});
    writetable(seqNmsPhg,['..' f 'Tables' f 'seqNmsPhg.xlsx']);
else
    seqNmsPhg = readtable(['..' f 'Tables' f 'seqNmsPhg.xlsx']);
end
%% Collect indel data by parsing breseq annotated.gd output files
breseqOutBac = [Sequencing_analysis_output f 'indels_bac.xlsx'];
breseqOutPhg = [Sequencing_analysis_output f 'indels_phg.xlsx'];
% Bacteria:
if ~exist(breseqOutBac,'file')
    indelBac = parse_breseq_output(breseqResBac,breseqOutBac,'bac');
else
    indelBac = readtable(breseqOutBac,'PreserveVariableNames',true);
end
% Phage:
if ~exist(breseqOutPhg,'file')
    indelPhg = parse_breseq_output(breseqResPhg,breseqOutPhg,'phg');
else
    indelPhg = readtable(breseqOutPhg,'PreserveVariableNames',true);
end
%% read manual indel file:
indelManOutBac = [Sequencing_analysis_output f 'indels_man_bac.xlsx'];
indelManOutPhg = [Sequencing_analysis_output f 'indels_man_phg.xlsx'];
% Bacteria:
opts = detectImportOptions(indelManOutBac);
opts = setvartype(opts, 'isolates', 'char');
opts.PreserveVariableNames = true;
indelManBac = readtable(indelManOutBac,opts);
% Phage:
opts = detectImportOptions(indelManOutPhg);
opts = setvartype(opts, 'isolates', 'char');
opts.PreserveVariableNames = true;
indelManPhg = readtable(indelManOutPhg,'PreserveVariableNames',true);
%% Collect CNV data (amplifications) from the manually written cnv.xlsx files
cnvOutBac = [Sequencing_analysis_output f 'cnv_bac.xlsx'];
cnvOutPhg = [Sequencing_analysis_output f 'cnv_phg.xlsx'];
% Bacteria:
opts = detectImportOptions(cnvOutBac);
opts = setvartype(opts, 'isolates', 'char');
opts.PreserveVariableNames = true;
cnvBac = readtable(cnvOutBac,opts);
% Phage:
opts = detectImportOptions(cnvOutPhg);
opts = setvartype(opts, 'isolates', 'char');
opts.PreserveVariableNames = true;
cnvPhg = readtable(cnvOutPhg,'PreserveVariableNames',true);
%% Collect SNPs data from BSAT mutmat1.csv output file
% Bacteria:
bsatOutBac = [Sequencing_analysis_output f 'snps_bac.csv'];
bsatOutPhg = [Sequencing_analysis_output f 'snps_phg.csv'];
mutmatcsv = ['TreeNum1' f 'MutMat' f 'MutMat1.csv'];
if ~exist(bsatOutBac,'file')
    copyfile([seq_ana_bac f mutmatcsv] ,bsatOutBac)
end
opts = detectImportOptions(bsatOutBac);
opts.DataLines = [2 Inf];
opts.VariableNamesLine = 1;
opts.Delimiter = '\t';
opts.PreserveVariableNames = false;
snpBac = readtable(bsatOutBac,opts);
% Phage:
if ~exist(bsatOutPhg,'file')
    copyfile([seq_ana_phg f mutmatcsv] ,bsatOutPhg)
end
opts = detectImportOptions(bsatOutPhg);
opts.DataLines = [2 Inf];
opts.VariableNamesLine = 1;
opts.Delimiter = '\t';
opts.PreserveVariableNames = false;
snpPhg = readtable(bsatOutPhg,opts);

%% Make a readable gene list for phage T7 (original gene names contain long numbers)
if ~exist(phage_readable_genes,'file')
    pg = readablePhageGeneList(phage_readable_genes);
else
    pg = readtable(phage_readable_genes);
end

%% Build a genotype matrix of all mutation types for each organism:
bac_MGV = load([seq_ana_bac f 'MutGenVCF.mat']);
% Bacteria:
[mm_uni_bac,mm_org_bac,geneCNVuniPos_bac] = buildMutMatrix(refGB_bac,snpBac,...
    indelBac,indelManBac,cnvBac,find(contains(seqNmsBac.Nms,WT_bac)),...
    bac_MGV,minCov(1),polymorph);
% Phage:
phg_MGV = load([seq_ana_phg f 'MutGenVCF.mat']);
wts = find(contains(seqNmsPhg.Nms,WT_phg(1:end-2)));
[mm_uni_phg,mm_org_phg,geneCNVuniPos_phg]= buildMutMatrix(refGB_phg,snpPhg,...
    indelPhg,indelManPhg,cnvPhg,find(contains(seqNmsPhg.Nms,WT_phg)),...
    phg_MGV,minCov(2),polymorph,wts);
%% Collect mutation information into tables:
snpInfoBac = buildSNPtbl(refGB_bac,bac_MGV,snpBac);
snpInfoPhg = buildSNPtbl(refGB_phg,phg_MGV,snpPhg,pg);
indelInfoBac = buildINDELtbl(refGB_bac,[indelBac ;indelManBac]);
indelInfoPhg = buildINDELtbl(refGB_phg,[indelPhg ;indelManPhg],pg);
cnvOrgInfoBac = buildCNVtbl(refGB_bac,cnvBac,0);
cnvOrgInfoPhg = buildCNVtbl(refGB_phg,cnvPhg,0);
cnvBacUni = table(geneCNVuniPos_bac(:,1),geneCNVuniPos_bac(:,2),'VariableNames',{'start','end'});
cnvUniInfoBac = buildCNVtbl(refGB_bac,cnvBacUni,1);
% Collect:
mutOrgInfoBac = [snpInfoBac; indelInfoBac; cnvOrgInfoBac];
mutUniInfoBac = [snpInfoBac; indelInfoBac; cnvUniInfoBac];
mutOrgInfoPhg = [snpInfoPhg; indelInfoPhg; cnvOrgInfoPhg];
mutUniInfoPhg = mutOrgInfoPhg; % No CNVs in the phage

%% Mutation matrix of mutations not present in the ancestor:
mm_uni_bac_diff = double((mm_uni_bac - mm_uni_bac(:,strcmp(seqNmsBac.Nms,WT_bac)))~=0);
mm_uni_phg_diff = double((mm_uni_phg - mm_uni_phg(:,strcmp(seqNmsPhg.Nms,WT_phg)))~=0);
mm_org_bac_diff = double((mm_org_bac - mm_org_bac(:,strcmp(seqNmsBac.Nms,WT_bac)))~=0);
mm_org_phg_diff = double((mm_org_phg - mm_org_phg(:,strcmp(seqNmsPhg.Nms,WT_phg)))~=0);

%% Add copy number information to the unified CNV mut matrix (only Bac), according to normalized DPs
mm_uni_bac_diff_wcnv = mm_uni_bac_diff;
minClusterSize = 2;
pos_margin = 120;
dp_mat_location = ['..' f 'script_data'];
numCluster = 4;
jump = 10;
dpnn_bac = build_dpnn(1,jump,seq_DPs_bac,seqNmsBac,...
    numCluster,minClusterSize,pos_margin,dp_mat_location,...
    mutUniInfoBac);
amps = find(strcmp(mutUniInfoBac.mutType,'CNV'));
for amp = 1:numel(amps)
    isoBac = find(mm_uni_bac_diff(amps(amp),:));
    for iso = 1:numel(isoBac)
        ampCov = dpnn_bac(round(mutUniInfoBac.genomePos1(amps(amp))/jump):...
            round(mutUniInfoBac.genomePos2(amps(amp))/jump),isoBac(iso));
        ampCov(ampCov==Inf) = NaN;
        mm_uni_bac_diff_wcnv(amps(amp),isoBac(iso)) = ceil(nanmean(ampCov)) - 1;
    end
end

%% Read isolate names of the cross-infection assay
infect_nms = [];
[~,~,infect_nms.bac] = xlsread(infection_name_order_file,1);
[~,~,infect_nms.phg1] = xlsread(infection_name_order_file,2);
[~,~,infect_nms.phg2] = xlsread(infection_name_order_file,3);
infect_nms.phg = [strcat('Sample_Phage1_',infect_nms.phg1); strcat('Sample_Phage2_',infect_nms.phg2) ];
phg_wt_pos = [24,86,[24,86]+96]; % positions of WT phage
infect_nms.phg(phg_wt_pos) = {WT_phg};

%% Exclude "bad" isolates from analysis and define indices to match betweem infection and sequencing isolate order:
[~,imBac] = ismember(strcat(prefix_Bac,infect_nms.bac), seqNmsBac.Nms);
% Bacterial isolate was CFP
isCFPbac = [22,23,24,79,80]; 
imBac(isCFPbac) = 0;
[~,imPhg] = ismember(infect_nms.phg, seqNmsPhg.Nms);
% No inoculated phage
imPhg([infectStruct1.inoculation(:);infectStruct2.inoculation(:)]==-1) = 0;
% Phage infecting CFP bac isolates
isCFPphg = [22,23,79,80,118,119,175,176];
imPhg(isCFPphg) = 0;
% Isolates with low sequencing quality
lq_bac = textscan(fopen([seq_ana_bac f 'low_quality_isolates.txt'],'r'),'%s');
lq_bac_pos = find(ismember(strcat(prefix_Bac,infect_nms.bac), lq_bac{1}));
imBac(lq_bac_pos) = 0;
lq_phg = textscan(fopen([seq_ana_phg f 'low_quality_isolates.txt'],'r'),'%s');
lq_phg_pos = ismember(infect_nms.phg, lq_phg{1});
imPhg(lq_phg_pos) = 0;
% Phages isolated from non-clean plaques
not_clean_phg = xlsread(not_clean_phg_xls,1);
imPhg(not_clean_phg(not_clean_phg(:,3)==1,1)+(not_clean_phg(not_clean_phg(:,3)==1,2)-1)*96) = 0;
noPhgApplied = 72+96;
imPhg(noPhgApplied) = 0;
% Isolates with heterogeneity in sequencing:
[~,isoP] = ind2sub(size(mm_uni_phg),find(arrayfun(@(x) strcmp(x, 'n'),mm_uni_phg)));
[Ap,Bp] = groupcounts(isoP);
hetroPhg = Bp(Ap>=hetroCutoff);
[~,isoB] = ind2sub(size(mm_uni_bac),find(arrayfun(@(x) strcmp(x, 'n'),mm_uni_bac)));
[Ab,Bb] = groupcounts(isoB);
hetroBac = Bb(Ab>=hetroCutoff);
imBac(ismember(imBac,hetroBac)) = 0;
imPhg(ismember(imPhg,hetroPhg)) = 0;
clear Ab Bb Ap Bp isoP isoB 

pBac = find(imBac~=0); % bac phenotype
pPhg = find(imPhg~=0); % phage phenotype
nBac = length(pBac); 
nPhg = length(pPhg);
gBac = imBac(imBac~=0); % bac genotype
gPhg = imPhg(imPhg~=0); % bac genotype

%% Apply indices on full infection matrix, perform log on the killing score for the figure
tmp = unique(sort(reshape(fullInfectionMat(pBac,pPhg),1,[])));
log_param = tmp(2);
infectionMatrixLog = log(fullInfectionMat(pBac,pPhg)+log_param); 
%% For the infection matrix: calculate pairwise distances and dendrogram order 
pd_bac = pdist(infectionMatrixLog,'naneucdist');
clus_bac = cluster(linkage(pd_bac,linkage_method),'MaxClust',12);
[~,~,orderBac] = dendrogram(linkage(pd_bac,linkage_method),inf,'Orientation','top','ColorThreshold',11) ;
pd_phg = pdist(infectionMatrixLog','naneucdist');
clus_phg = cluster(linkage(pd_phg,linkage_method),'MaxClust',12);
[~,~,orderPhg] = dendrogram(linkage(pd_phg,linkage_method),inf,'Orientation','left','ColorThreshold',11);
close all

%% For comparison to the standard turbidity:
turbidityMat = fullTurbidityMat(pBac,pPhg);
pd_bac_turb = pdist(turbidityMat,'naneucdist');
pd_phg_turb = pdist(turbidityMat','naneucdist');
clus_bac = cluster(linkage(pd_bac_turb,linkage_method),'MaxClust',12);
[~,~,orderBacTurbid] = dendrogram(linkage(pd_bac_turb,linkage_method),inf,'Orientation','top','ColorThreshold',11) ;
clus_phg = cluster(linkage(pd_phg_turb,linkage_method),'MaxClust',12);
[~,~,orderPhgTurbid] = dendrogram(linkage(pd_phg_turb,linkage_method),inf,'Orientation','left','ColorThreshold',11);
close all
%% Order mutations by genome position:
% Unified CNV
[~,geneOrderBacUni] = sort(mutUniInfoBac.genomePos1);
[~,geneOrderPhgUni] = sort(mutUniInfoPhg.genomePos1);
hasIsoMutsBacUni = any(mm_uni_bac_diff(geneOrderBacUni,gBac(orderBac))');
geneOrderBacUni = geneOrderBacUni(hasIsoMutsBacUni);
hasIsoMutsPhgUni = any(mm_uni_phg_diff(geneOrderPhgUni,gPhg(orderPhg))');
geneOrderPhgUni = geneOrderPhgUni(hasIsoMutsPhgUni);
% Original CNV
[~,geneOrderBacOrg] = sort(mutOrgInfoBac.genomePos1);
[~,geneOrderPhgOrg] = sort(mutOrgInfoPhg.genomePos1);
hasIsoMutsBacOrg = any(mm_org_bac_diff(geneOrderBacOrg,gBac(orderBac))');
geneOrderBacOrg = geneOrderBacOrg(hasIsoMutsBacOrg);
hasIsoMutsPhgOrg = any(mm_org_phg_diff(geneOrderPhgOrg,gPhg(orderPhg))');
geneOrderPhgOrg = geneOrderPhgOrg(hasIsoMutsPhgOrg);

%%
%%%%%%%%%%%%%%%%%%%
%%% Save Tables %%%
%%%%%%%%%%%%%%%%%%%
%% Phenotype:
% Infectivity:
phenotype = table(fullInfectionMat(pBac(orderBac),pPhg(orderPhg)));
writetable(phenotype,phenTblFile,'WriteVariableNames', 0);
% Turbidity:
turbidityTbl = table(fullTurbidityMat(pBac(orderBacTurbid),pPhg(orderPhgTurbid)));
writetable(turbidityTbl,turbidTblFile,'WriteVariableNames', 0);
turbidityTblOriginalOrder = table(fullTurbidityMat(pBac(orderBac),pPhg(orderPhg)));
writetable(turbidityTblOriginalOrder,turbidTblOrgOrderFile,'WriteVariableNames', 0);

%% Write isolate information table
SamplingInfo = readtable(SamplingInfo_file); % Table was constructed with organizeSamplesOnPlates.m
% NOTE: coordinates are only relative to the ***cropped*** images files
% (without spaces between plates)
% Bacteria:
BacIsoTbl = buildBacIsoTbl(SamplingInfo,pltNms,gBac,pBac,orderBac,WT_bac,seqNmsBac);
writetable(BacIsoTbl,BacIsoTblFile)    
% Phages:
PhgIsoTbl = buildPhgIsoTbl(SamplingInfo,pltNms,gPhg,pPhg,orderPhg,WT_phg,seqNmsPhg,BacIsoTbl);
writetable(PhgIsoTbl,PhgIsoTblFile)    

%% Write Original Bac genotypes:
% Bacteria:
genotypes_bac_org = table(mm_org_bac_diff(geneOrderBacOrg,gBac(orderBac)),'VariableNames',{'genotype'});
preexisting_loci = table(zeros(size(genotypes_bac_org)),'VariableNames',{'preexisting_loci'});
BacGenTblOrg = [mutOrgInfoBac(geneOrderBacOrg,:),preexisting_loci, genotypes_bac_org];
writetable(BacGenTblOrg,BacGenTblOrgFile)

genotypes_bac_uni = table(mm_uni_bac_diff_wcnv(geneOrderBacUni,gBac(orderBac)),'VariableNames',{'genotype'});
preexisting_loci = table(zeros(size(genotypes_bac_uni)),'VariableNames',{'preexisting_loci'});
BacGenTblUni = [mutUniInfoBac(geneOrderBacUni,:),preexisting_loci, genotypes_bac_uni];
writetable(BacGenTblUni,BacGenTblUniFile)

% Phage:
loci = readtable(phage_preexisting_loci_File);
genotypes_phg_org = table(mm_org_phg_diff(geneOrderPhgOrg,gPhg(orderPhg)),'VariableNames',{'genotype'});
preexisting_loci = table(ismember(mutOrgInfoPhg.genomePos1(geneOrderPhgOrg,:),...
    loci.phage_preexisting_loci),'VariableNames',{'preexisting_loci'});
PhgGenTblOrg = [mutOrgInfoPhg(geneOrderPhgOrg,:), preexisting_loci, genotypes_phg_org];
writetable(PhgGenTblOrg,PhgGenTblOrgFile)

genotypes_phg_uni = table(mm_uni_phg_diff(geneOrderPhgUni,gPhg(orderPhg)),'VariableNames',{'genotype'});
preexisting_loci = table(ismember(mutUniInfoPhg.genomePos1(geneOrderPhgUni,:),...
    loci.phage_preexisting_loci),'VariableNames',{'preexisting_loci'});
PhgGenTblUni = [mutUniInfoPhg(geneOrderPhgUni,:) ,preexisting_loci, genotypes_phg_uni];
writetable(PhgGenTblUni,PhgGenTblUniFile)

%% Mutation: load relevant pathways
if ~exist(PathwayBacUni,'file')
    fprintf('No bac uni pathways file available, writing mutations to xls\n');
    pwtbl = table(BacGenTblUni.mutName,cell(size(BacGenTblUni.mutName)),'VariableNames',{'mutation','pathway'});
    writetable(pwtbl,PathwayBacUni)
end
if ~exist(PathwayBacOrg,'file')
    fprintf('No bac org pathways file available, writing mutations to xls\n');
    pwtbl = table(BacGenTblOrg.mutName,cell(size(BacGenTblOrg.mutName)),'VariableNames',{'mutation','pathway'});
    writetable(pwtbl,PathwayBacOrg)
end
if ~exist(PathwayPhgUni,'file')
    fprintf('No phage uni pathways file available, writing mutations to xls\n');
    pwtbl = table(PhgGenTblUni.mutName,cell(size(PhgGenTblUni.mutName)),'VariableNames',{'mutation','pathway'});
    writetable(pwtbl,PathwayPhgUni)
end
if ~exist(PathwayPhgOrg,'file')
    fprintf('No phage org pathways file available, writing mutations to xls\n');
    pwtbl = table(PhgGenTblOrg.mutName,cell(size(PhgGenTblOrg.mutName)),'VariableNames',{'mutation','pathway'});
    writetable(pwtbl,PathwayPhgOrg)
end
toc