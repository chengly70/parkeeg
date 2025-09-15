# parkeeg
Code to calc and display results for EEG Motor Cortex Parkinson's Paper
Relies on publicly available data from UCSD at: https://openneuro.org/datasets/ds002778/versions/1.0.5  
All 46 bdf files from this openneuro site should downloaded and put into a sub-directory named "dataUCSD"
All computer code is in MATLAB

Many other figures generated for this project that were too cumbersome to put in the paper are here: https://drive.google.com/drive/folders/1QdOytuiETdtfEhSMKIuVl8vTBzvteLXZ?usp=sharing

--------- Functions to process EEG data ------------

getData.m -- code to load and process raw EEG data, relies on eeg_read_pdf.m file by Gleb Tcheslavski, via https://www.mathworks.com/matlabcentral/fileexchange/13070-eeg-bdf-reader 
This script creates dataFiltrd.mat that is used in sub-sequent methods. 

calcBandpssData.m -- code to load and process bandpass filtering & envelope extraction of EEG data, final form of EEG data used. Save results in dataEEG_Ct[Mn/k2/k3].mat where Mn=mean, k2=C3 (left prim.motor.cortex), k3=C4 (right prim.motor.cortex). 
Relies on getContn.m (our fcn to downsample in diff time windows) and bandpass.m (must install `Filter Analyzer' or `Signal Analyzer' MATLAB toolbox). 

getContn.m -- our downsampling code, used in previous function

get_sideSeverity.m -- code to get physical impairment side of each Parkinson's patient, data from participants.tsv file at openneuro.org website above. Saves struct whSide in dSeveritySide.mat

------- Scripts/Functions for Main Results in the Paper ----

ACF results (Fig 2): from script getAuto_threshTime_bySide.m, calculates individual and population ACF, as well as time ACF cross below prescribed threshold (0.1), saved in struct variable CharTim

relies on DFA_fun.m file by Martin Magris, via https://www.mathworks.com/matlabcentral/fileexchange/67889-detrended-fluctuation-analysis-dfa

------- MAT files --------

dataFiltrd.mat -- extracted EEG data, some subjects have SMALL snippets at end removed because they appear to be artifacts
dSeveritySide.mat -- denotes the side of physical impairment in Parkinson's patients, created from get_sideSeverity.m script
