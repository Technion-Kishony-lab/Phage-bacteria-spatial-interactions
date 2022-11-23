
# Multistep diversification in spatiotemporal bacterial-phage coevolution
# Data analysis pipeline  

System requirements:
All scripts were coded with MATLAB R2020a on a Windows 10 Enterprise 64-bit Operating System. 
MATLAB Statistics and Machine Learning Toolbox is required.

Instructions: 
1. Download all files from https://github.com/Technion-Kishony-lab/Phage-bacteria-spatial-interactions.git and place in a "scripts" folder.
2. Extract "Additional_scripts".
3. Download all raw data files from Zenodo:
	10.5281/zenodo.7347986
About 60GB of storage space is required.
4. Place all compressed files, except for the "script_data" file, under a folder named "input_files".
5. extract all folders.
Also, extract all images folders (into their current location):
input_files\coevolution\20191206_plate_raw_imgs
input_files\coevolution\20191206_plate_raw_imgs\20191218_continue
input_files\coevolution\220904_MGY_SP
input_files\crossInfection\coevo_cross_infection\20200205_phage1_interactions
input_files\crossInfection\coevo_cross_infection\20200202_phage2_interactions

6. extract "script_data" (optional)
Place "script_data", "input_files", and "scripts" folders under the same folder (the three folders should be at the same level) 

The raw data folders contain:
* Raw images of the coevolution experiment 
* Raw images of the cross-infection assay and additional technical data.
* Raw and analyzed sequencing analysis files 
* Illustration of the experimental setup (for Fig. 1a and Supplementary Fig. 1a) 
* Isolate sampling coordinates, cropping coordinate of the full coevolution images.
* Genetic pathways information
* Raw data of the mlaA and opgG mutants assay.
* Example images of a two-color plaque assay on a petri dish.
 
7. Run the analysis and plot figures:
In order to analyze the data and produce figures and tables, stand inside the scripts folder and run the script
main.m

The main script will organize the data into tables (in 'Tables' folder) and mat files (in 'script_data' folder).
These files are later read by the script as part of the analysis and the generation of figures. 
The script will generate and save all main and supplementary figures under the 'scripts\Figs' folder.
Additional calculation results (median number of peaks, BiMat analysis, dN/dS ratio, the Lasso model and the multi-mutated gene simulation results)
are saved under scripts\output_files

The script can run without the script_data folder and generate all mat files from scratch. 
However, there are several time consuming steps in the full analysis:
a. Generation of cropped images for the coevolution time-lapse movies.
b. Conversion of coevolution images into matrices and peaks analysis. 
c. Conversion of cross-infection images and experimental data into structs.
d. Cross-infection plate binarization: for infectivity score calculation, plaque regions are detected according to a threshold. 
e. Conversion of the mlaA and opgG assay images into matrices
After the first run, these steps will be skipped.

Runtime for the initial run (without script_data): ~3 hours. 
Runtime after the first run (or with script_data): ~30min, with most of the time spent on plotting Supplementary Figure 11. 

