% Main analysis script for "Multistep diversification in spatiotemporal
% bacterial-phage coevolution"
% Produces all manuscript figures and tables
% Run from the 'scripts' folder
tic
fprintf('Analyzing E.coli - T7 co-evolution...\n')
analysis % main analysis
clear
close all
fprintf('Printing main figures...\n')
Figure1_SuppFigure1_3 % Coevolution experiment analysis and peak analysis + Intensity peak distribution in coevolution and no-phage control.
Figure2 % Infection matrix
Figure3 % Lasso regression  
Figure4 % Two color plaque assay with mlaA and opgG mutants 
fprintf('Printing supplementary figures...\n')
SuppFigure4 % Maps of the coevolution sampled sites.
SuppFigure5 % Examples of phage isolation via two-color plaque assay
SuppFigure6 % Cross infection measurement via turbidity score
SuppFigure7 % BiMat analysis for nestedness and modularity 
SuppFigure8 % dN/dS 
SuppFigure9 % Multi-mutated genes simulation
SuppFigure10 %  mlaA and opgG mutants: dye swap
%%
close all
clear
SuppFigure11 % Images of all phage-bacteria interactions
toc