# parkeeg
Code to calc and display results for EEG Motor Cortex Parkinson's Paper
Relies on publicly available data from UCSD at: https://openneuro.org/datasets/ds002778/versions/1.0.5  
All 46 bdf files from this openneuron site were downloaded and put into a sub-directory named "dataUCSD"

getData.m -- MATLAB code to load and process raw EEG data, relies on eeg_read_pdf.m file by Gleb Tcheslavski, gleb@vt.edu via https://www.mathworks.com/matlabcentral/fileexchange/13070-eeg-bdf-reader 
This script creates dataFiltrd.mat that is used in sub-sequent methods. 
