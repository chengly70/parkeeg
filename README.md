# parkeeg
Code to calc and display results for EEG Motor Cortex Parkinson's Paper
Relies on publicly available data from UCSD at: https://openneuro.org/datasets/ds002778/versions/1.0.5  
All 46 bdf files from this openneuro site were downloaded and put into a sub-directory named "dataUCSD"

Many other figures generated for this project that were too cumbersome to put in the paper are here: https://drive.google.com/drive/folders/1QdOytuiETdtfEhSMKIuVl8vTBzvteLXZ?usp=sharing

getData.m -- MATLAB code to load and process raw EEG data, relies on eeg_read_pdf.m file by Gleb Tcheslavski, via https://www.mathworks.com/matlabcentral/fileexchange/13070-eeg-bdf-reader 
This script creates dataFiltrd.mat that is used in sub-sequent methods. 
