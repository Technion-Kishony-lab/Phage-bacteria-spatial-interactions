
# Multistep diversification in spatiotemporal bacterial-phage coevolution: Data analysis pipeline  

System requirements:
All scripts were coded with MATLAB R2020a on a Windows 10 Enterprise 64-bit Operating System. 
MATLAB Statistics and Machine Learning Toolbox is required.

Instructions: 
1. Download all files from https://github.com/Technion-Kishony-lab/Phage-bacteria-spatial-interactions.git and place in a "scripts" folder. 
2. Extract "Additional_scripts".
	Git repositiry content should be directly under "scripts". (for example, scripts\main.m , scripts\Additional_scripts\...)
3. Download all raw data files from Zenodo:
	10.5281/zenodo.7347986
	About 60GB of storage space is required. 
	
	The raw data folders contain:
	* Raw images of the coevolution experiment 
	* Raw images of the cross-infection assay and additional technical data.
	* Raw and analyzed sequencing analysis files 
	* Illustration of the experimental setup (for Fig. 1a and Supplementary Fig. 1a) 
	* Isolate sampling coordinates, cropping coordinate of the full coevolution images.
	* Genetic pathways information
	* Raw data of the mlaA and opgG mutants assay.
	* Example images of a two-color plaque assay on a petri dish.

	Follow the next steps to organize the data:
	
4. Place all compressed files, except for the "script_data" file, under a new folder named "input_files".
5. Make a folder named "coevolution" under "input_files".
	Place the following compressed files under "coevolution" and extract as follows:
	a. "coevolution_additionalFiles.zip" - extract files directly into "coevolution"
	b. "220904_MGY_SP.zip" - extract images into the folder "220904_MGY_SP"
	c. "20191206_plate_raw_imgs" - extract images into the folder "20191206_plate_raw_imgs" 
	d. Place the compressed file "20191218_continue.zip" under the folder "20191206_plate_raw_imgs" (see c.)
	   extract images in "20191218_continue.zip" into the folder "20191218_continue" (final hierarchy: input_files\coevolution\20191206_plate_raw_imgs\20191218_continue\*.CR2)

6. Extract all other folders in "input_files".
	Also, extract the remaining images folders:
	input_files\crossInfection\coevo_cross_infection\20200205_phage1_interactions
	input_files\crossInfection\coevo_cross_infection\20200202_phage2_interactions

7. Extract "script_data" (optional)
Place "script_data", "input_files", and "scripts" folders under the same folder (the three folders should be at the same level) 

8. Run the analysis and plot figures:
In order to analyze the data and produce figures and tables, stand inside the "scripts" folder and run the script
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

