% calc resulting data after bandpassing, coarsegraining, & envelope, but using loading dataFiltrd.mat
% Uses fcns bandpass.m and getContn.m
% !!!! saving result in  __ go to lines 128-133 to comment in/out to set
% corresopnding file name. AND change line 16!!!!!! 

clear

load dataFiltrd.mat

nH=size(ctData.hlt,1);
nP=size(ctData.prkOff,1); %same as using *On*

deltaTVec=[dt; .004; (.005:.0025:.1)']; % maybe downsample, sampling freq is large
lenDT=length(deltaTVec);

% which electrodes test, or how process data?
k=[2 3];%[1 2 3 4]; %1=FC1, 2=C3, 3=C4 , 4=FC2

%which data use?
%whichD=rawData;   Lt=nTim;
whichD=ctData;   Lt=nTct;

lenF=5; %total ways to BP, 5 freq bands

% --- outputs ---
eegD=struct('hlt',[],'prkOff',[],'prkOn',[]); %time-average d2 with Sam's myYuleWalker.m

eegD.hlt=cell(nH,lenDT,lenF);
eegD.prkOff=cell(nP,lenDT,lenF);
eegD.prkOn=cell(nP,lenDT,lenF);

tic

for fInd=1:lenF

    switch fInd
        case 1
            pssLow=1; %in Hz
            pssHigh=4;
        case 2
            pssLow=4;
            pssHigh=8;
        case 3
            pssLow=8;
            pssHigh=13;
        case 4
            pssLow=13;
            pssHigh=30;
        case 5
            pssLow=30;
            pssHigh=50;
    end
    wPss=[pssLow pssHigh];
    nmTrim= round(1.25/dt)+1; %number of bins to remove start&end b/c bandpass transients

    for dtInd=1:lenDT

        for j=1:nH
            if(isscalar(k))
                tmpDat=whichD.hlt{j}(k,:);
            else %take simple avg of all channels
                tmpDat=mean(whichD.hlt{j}(k,:)); %take avg over channels first
            end

            %do bandpassing
            tmpDat=bandpass(tmpDat-mean(tmpDat),wPss,fs);
            % extract envelope
            tmpDat=abs(hilbert(tmpDat));

            tmpDat=tmpDat(nmTrim+1:end-nmTrim+1); %trim transients start&end

            %down-sample to match window (dtInd)
            windows=[0 dt*(length(tmpDat)-1)];
            contM=[(windows(1):dt:windows(2))' tmpDat'];
            tmpDat = getContn(contM,windows,deltaTVec(dtInd)); %down-sample here

            eegD.hlt{j,dtInd,fInd}=tmpDat;
        end

        % Parkinson's Off of meds
        for j=1:nP
            if(isscalar(k))
                tmpDat=whichD.prkOff{j}(k,:);
            else %take simple avg of all channels
                tmpDat=mean(whichD.prkOff{j}(k,:)); %take avg over channels first
            end
            
            %band pass
            tmpDat=bandpass(tmpDat-mean(tmpDat),wPss,fs);
            % extract envelope
            tmpDat=abs(hilbert(tmpDat));

            tmpDat=tmpDat(nmTrim+1:end-nmTrim+1); %trim transients start&end

            %down-sample
            windows=[0 dt*(length(tmpDat)-1)];
            contM=[(windows(1):dt:windows(2))' tmpDat'];
            tmpDat = getContn(contM,windows,deltaTVec(dtInd)); %down-sample here

            eegD.prkOff{j,dtInd,fInd}=tmpDat;
        end

        % Parkinson's ON medication
        for j=1:nP
            if(isscalar(k))
                tmpDat=whichD.prkOn{j}(k,:);
            else %take simple avg of all channels
                tmpDat=mean(whichD.prkOn{j}(k,:)); %take avg over channels first
            end
            
            %band pass
            tmpDat=bandpass(tmpDat-mean(tmpDat),wPss,fs);
            % extract envelope
            tmpDat=abs(hilbert(tmpDat));

            tmpDat=tmpDat(nmTrim+1:end-nmTrim+1); %trim transients start&end

            %down-sample
            windows=[0 dt*(length(tmpDat)-1)];
            contM=[(windows(1):dt:windows(2))' tmpDat'];
            tmpDat = getContn(contM,windows,deltaTVec(dtInd)); %down-sample here

            eegD.prkOn{j,dtInd,fInd}=tmpDat;
        end
    end

    %save dataEEG_Ctk2 eegD deltaTVec
    %save dataEEG_Ctk3 eegD deltaTVec
    save dataEEG_CtMn eegD deltaTVec
    
end

toc
