% script to process UCSD (Aron Lab) EEG data. The only processing we do is
% cut-off a small piece at the end where there are artifacts: all electrodes appear to shift together
% result is saved in ctData (also calc rawData for internal testing)

cd dataUCSD/
allF=dir('*.bdf'); %get all .bdf file names
fileHealth=dir('sub-hc*.bdf');          
fileParkOffMed=dir('sub-pd*off*.bdf');  
fileParkOnMed=dir('sub-pd*on*.bdf');
cd ..

preF='dataUCSD/';

nH=length(fileHealth); %number of healthy subjects
nP=length(fileParkOffMed); % #Parkinsons; same as using *On* Medication

%[data,numChan,labels,txt,fs,gain,prefiltering,ChanDim]=eeg_read_bdf(flName,'all','n');

%hard-code cutoffs
cutfof=[0; 94612; 0; 0; 0; 0; 0; 0; 97500; 95500; 0; 0; 96600; 88890; 95560; 96300; ...
    0; 0; 0; 95312; 0; 0; 0; 94000; 96000; 95510; 70000; 0; 0; 94500; 0; ...
    0; 95800; 0; 0; 95500; 95100; 0; 90200; 0; 97000; 0; 0; 104100; 0; 94300]; %0 denotes NO cut-off

% -- OUTPUTS ----
kID=[5;8;23;26];   %IDs of channels/electrodes to save; 8 & 23 are C3 and C4 (left & right motor cortex)
prefilts=cell(nP*2+nH,1);
gains=zeros(nP*2+nH,1);
nTim=struct('hlt',[],'prkOff',[],'prkOn',[]); %length of (raw) data
nTct=struct('hlt',[],'prkOff',[],'prkOn',[]); % length of cut data (remove ends on some)
rawData=struct('hlt',[],'prkOff',[],'prkOn',[]); %raw
ctData=struct('hlt',[],'prkOff',[],'prkOn',[]);    %just cutoff at end (on some)

rawData.hlt=cell(nH,1); nTim.hlt=zeros(nH,1); nTct.hlt=zeros(nH,1); 
ctData.hlt=cell(nH,1);
for j=1:nH
    [data,~,labels,~,fs,gain,prefiltering]=eeg_read_bdf([preF,fileHealth(j).name],'all','n');
    nTim.hlt(j,1)=size(data,2); gains(j,1)=gain;  prefilts{j,1}=prefiltering; 
    rawData.hlt{j,1}=data(kID,:);
    if(cutfof(j)~=0)
        ctData.hlt{j,1}=data(kID,1:cutfof(j));
    else
        ctData.hlt{j,1}=rawData.hlt{j,1};
    end
    nTct.hlt(j,1)=size(ctData.hlt{j},2);
end

rawData.prkOff=cell(nP,1);  rawData.prkOn=cell(nP,1);
ctData.prkOff=cell(nP,1);    ctData.prkOn=cell(nP,1);
nTim.prkOff=zeros(nP,1); nTct.prkOff=zeros(nP,1); 
nTim.prkOn=zeros(nP,1); nTct.prkOn=zeros(nP,1); 
for j=1:nP
    [data,~,labels,~,fs,gain,prefiltering]=eeg_read_bdf([preF,fileParkOffMed(j).name],'all','n');
    nTim.prkOff(j,1)=size(data,2); gains(j+nH,1)=gain;  prefilts{j+nH,1}=prefiltering;
    rawData.prkOff{j,1}=data(kID,:);
    if(cutfof(j+nH)~=0)
        ctData.prkOff{j,1}=data(kID,1:cutfof(j+nH));
    else
        ctData.prkOff{j,1}=rawData.prkOff{j,1};
    end
    nTct.prkOff(j,1)=size(ctData.prkOff{j},2);
    %-- repeat for ParkOn ---
    [data,~,labels,~,fs,gain,prefiltering]=eeg_read_bdf([preF,fileParkOnMed(j).name],'all','n');
    nTim.prkOn(j,1)=size(data,2); gains(j+nH+nP,1)=gain;  prefilts{j+nH+nP,1}=prefiltering;
    rawData.prkOn{j,1}=data(kID,:);
    if(cutfof(j+nH+nP)~=0)
        ctData.prkOn{j,1}=data(kID,1:cutfof(j+nH+nP));
    else
        ctData.prkOn{j,1}=rawData.prkOn{j,1};
    end
    nTct.prkOn(j,1)=size(ctData.prkOn{j},2);
end

dt=1/fs;

% save relevant results
save dataFiltrd labels fs gains prefilts kID dt rawData nTim ctData nTct 