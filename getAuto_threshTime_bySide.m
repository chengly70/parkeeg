% Autocorrelation function from processed data 
% Also, set threshold & calc times first go below as a measure of 'char timescale'
% dataEEG_[].mat from dataFiltrd.mat, see calcBandpssData.m: apply bandpass filtering & extracting envelope

%%

Thresh=0.1; %0.1 is good, so are 0.12 & 0.08 (roughly the same). 

ccD=[0 0 0; 1 0 0; 0 0 1];

fl_showInd=1;
fl_normACF=1;

load('dataEEG_Ctk2.mat','eegD')
eegD_r=eegD;
load('dataEEG_Ctk3.mat','eegD')
eegD_l=eegD;
load dataEEG_CtMn.mat %eegD by default is the mean

nH=size(eegD.hlt,1);
nP=size(eegD.prkOff,1);
lenDT=length(deltaTVec);
lenFD=5; 

load dSeveritySide.mat %whSide, 2=right impaired => use C3 (k=2), 1=left impaired => use C4 (k=3)
for j=1:nP
    switch whSide.prkOff(j) 
        case 1
            for dt_i=1:lenDT
                for f_i=1:lenFD
                    eegD.prkOff{j,dt_i,f_i}=eegD_l.prkOff{j,dt_i,f_i};
                    if(j==14) %14 is only subj with On diff Off, and happens to occur in case 1
                        eegD.prkOn{j,dt_i,f_i}=eegD_r.prkOn{j,dt_i,f_i};
                    else
                        eegD.prkOn{j,dt_i,f_i}=eegD_l.prkOn{j,dt_i,f_i};
                    end
                end
            end
        case 2
            for dt_i=1:lenDT
                for f_i=1:lenFD
                    eegD.prkOff{j,dt_i,f_i}=eegD_r.prkOff{j,dt_i,f_i};
                    eegD.prkOn{j,dt_i,f_i}=eegD_r.prkOn{j,dt_i,f_i};
                end
            end
    end
end


dt=deltaTVec(1);

%param for autocorr &
nmLags=10000;

% --- outputs ---
t_Aut=(0:nmLags)'*dt;    %time vector for acf
acfAll=struct('hlt',[],'prkOff',[],'prkOn',[]);  %autocov function
idxt=struct('hlt',[],'prkOff',[],'prkOn',[]);  %index of time-scale
CharTim=struct('hlt',[],'prkOff',[],'prkOn',[]);  %index of time-scale
% pVals checking CharTim --- Are time below thresh statSignif diff?
pWR=zeros(3,5);
efszWR=zeros(3,5);
pTs=zeros(3,5);
efszTs=zeros(3,5);
pA=zeros(3,5);
efszA=zeros(3,5);

CharTim.hlt=NaN(nH,5);
CharTim.prkOff=NaN(nP,5);
CharTim.prkOn=NaN(nP,5);
for fInd=1:lenFD

    acfAll.hlt=zeros(nH,nmLags+1);
    idxt.hlt=NaN(nH,1); %NOT saving indxt over fInd
    if(fl_showInd)
        h1=figure('Renderer', 'Painters'); hold on;
        h2=figure('Renderer', 'Painters'); hold on;
        h3=figure('Renderer', 'Painters'); hold on;
    end
    for j=1:nH
        tmpDat=eegD.hlt{j,1,fInd};
        mn=mean(tmpDat);

        acf=autocorr(tmpDat-mn,'NumLags',nmLags,'NumStd',1.96);
        acfAll.hlt(j,:)=acf;
        
        if(fl_showInd)
            figure(h1)
            plot(t_Aut,acf,'color',.7*ones(1,3))
        end
    end
    if(fl_showInd)
        figure(h1)
        plot(t_Aut,mean(acfAll.hlt),'k','LineWidth',2)
        plot(t_Aut,Thresh*ones(size(t_Aut)),'k--','LineWidth',.5)
        set(gca,'XScale','log')
        set(gca,'XLim',[t_Aut(1) t_Aut(end)])
    end

    acfAll.prkOff=zeros(nP,nmLags+1);
    idxt.prkOff=NaN(nP,1);
    for j=1:nP
        tmpDat=eegD.prkOff{j,1,fInd};
        mn=mean(tmpDat);

        acf=autocorr(tmpDat-mn,'NumLags',nmLags,'NumStd',1.96);
        acfAll.prkOff(j,:)=acf; %is a row vector
        if(fl_showInd)
            figure(h2)
            plot(t_Aut,acf,'color',.7*ones(1,3))
        end
    end
    if(fl_showInd)
        figure(h2)
        plot(t_Aut,sum(acfAll.prkOff)./(nP),'r','LineWidth',2)
        plot(t_Aut,Thresh*ones(size(t_Aut)),'k--','LineWidth',.5)
        set(gca,'XScale','log')
        set(gca,'XLim',[t_Aut(1) t_Aut(end)])
    end

    % ---- finally for Park-On med
    acfAll.prkOn=zeros(nP,nmLags+1);
    idxt.prkOn=NaN(nP,1);
    for j=1:nP
        tmpDat=eegD.prkOn{j,1,fInd};
        mn=mean(tmpDat);

        acf=autocorr(tmpDat-mn,'NumLags',nmLags,'NumStd',1.96);
        acfAll.prkOn(j,:)=acf; %is a row vector
        if(fl_showInd)
            figure(h3)
            plot(t_Aut,acf,'color',.7*ones(1,3))
        end
    end
    if(fl_showInd)
        figure(h3)
        plot(t_Aut,sum(acfAll.prkOn)./(nP),'b','LineWidth',2)
        plot(t_Aut,Thresh*ones(size(t_Aut)),'k--','LineWidth',.5)
        set(gca,'XScale','log')
        set(gca,'XLim',[t_Aut(1) t_Aut(end)])
    end

    % get mean ACF over all subjects; FITS for individ ACF are going to be bad, NEED aABC perhaps!
    mnAcf_h=mean(acfAll.hlt)';
    mnAcf_poff=sum(acfAll.prkOff)/(nP); mnAcf_poff=mnAcf_poff';
    mnAcf_pon=sum(acfAll.prkOn)/(nP); mnAcf_pon=mnAcf_pon';

    for j=1:nH
        tmp=find(acfAll.hlt(j,:)<Thresh,1);
        if(isscalar(tmp))
            idxt.hlt(j,1)=tmp;
            CharTim.hlt(j,fInd)=t_Aut(tmp); 
        end
    end
    for j=1:nP
        tmp=find(acfAll.prkOff(j,:)<Thresh,1);
        if(isscalar(tmp))
            idxt.prkOff(j,1)=tmp;
            CharTim.prkOff(j,fInd)=t_Aut(tmp);
        end
        tmp=find(acfAll.prkOn(j,:)<Thresh,1);
        if(isscalar(tmp))
            idxt.prkOn(j,1)=tmp;
            CharTim.prkOn(j,fInd)=t_Aut(tmp);
        end
    end

    %treating 'NaN'; ONLY occurs in subj 14 (ParkOn) in lowGamma, else
    %everything is ok
    hlt=CharTim.hlt(:,fInd); 
    prkOff=CharTim.prkOff(:,fInd);
    prkOn=CharTim.prkOn(:,fInd); prkOn(isnan(CharTim.prkOn(:,fInd)))=11.74; %at thres of 0.107 
    [pVals,Effsz]=getPvEs(hlt, prkOff, prkOn);
    pWR(:,fInd)=pVals.W; %store them all for each FreqBand
    pTs(:,fInd)=pVals.T;
    pA(:,fInd)=pVals.A;
    efszWR(:,fInd)=Effsz.W;
    efszTs(:,fInd)=Effsz.T;
    efszA(:,fInd)=Effsz.A;

end %ending the for fInd=1:lenFD(5) WAY at beginning!

figure
hold on
plot(1:5,pWR(1,:),'k-')
plot(1:5,pWR(2,:),'k--')
plot(1:5,pTs(1,:),'m-')
plot(1:5,pTs(2,:),'m--')
set(gca,'FontSize',18)
xlabel('Freq Band')
ylabel('p-values')
legend('WRST Hpoff','WRST Hpon','t- Hpoff','t- Hpon')
set(gca,'YScale','log')
figure
hold on
plot(1:5,pWR(3,:),'k-')
plot(1:5,pTs(3,:),'m-')
set(gca,'FontSize',18)
xlabel('Freq Band')
ylabel('p-values')
legend('WRST OffvOn','t- OffvOn')
set(gca,'YScale','log')

%% show box plots, table of p-vals & Effect Sizes
cc=[zeros(1,3); 1 0 0; 0 0 1];
dComp={'H vs Off';'H vs On';'Off vs On'}; %for showing table
h1=figure;
hold on
h2=figure;
hold on
for j=1:5 %use 4 (ignore lowGamma) or 5 (all)
    %throwing out the few times with 'NaN'
    hlt=CharTim.hlt(:,j); hlt=hlt(~isnan(CharTim.hlt(:,j)));
    prkOff=CharTim.prkOff(:,j); prkOff=prkOff(~isnan(CharTim.prkOff(:,j)));
    prkOn=CharTim.prkOn(:,j); prkOn=prkOn(~isnan(CharTim.prkOn(:,j)));
    [pVals,Effsz]=getPvEs(hlt, prkOff, prkOn);
    
    figure(h1)
    boxchart(((j-1)*3+1)*ones(length(hlt),1),hlt,'MarkerStyle','none','BoxFaceColor',cc(1,:),'WhiskerLineColor',cc(1,:))
    boxchart(((j-1)*3+2)*ones(length(prkOff),1),prkOff,'MarkerStyle','none','BoxFaceColor',cc(2,:),'WhiskerLineColor',cc(2,:))
    boxchart(((j-1)*3+3)*ones(length(prkOn),1),prkOn,'MarkerStyle','none','BoxFaceColor',cc(3,:),'WhiskerLineColor',cc(3,:))
    % plot( (j-1)*3+1, mean(hlt), '.','color',cc(1,:), 'MarkerSize', 24) %plot means on boxplots
    % plot( (j-1)*3+2, mean(prkOff), '.','color',cc(2,:), 'MarkerSize', 24)
    % plot( (j-1)*3+3, mean(prkOn), '.','color',cc(3,:), 'MarkerSize', 24)
    figure(h2)
    plot( j, mean(hlt), '.','color',cc(1,:), 'MarkerSize', 24)
    plot( j, mean(prkOff), '.','color',cc(2,:), 'MarkerSize', 24)
    plot( j, mean(prkOn), '.','color',cc(3,:), 'MarkerSize', 24)

        
    switch j
        case 1
            frqBand='delta [1,3]';
        case 2
            frqBand='theta [3,8]';
        case 3
            frqBand='alpha [8,13]';
        case 4
            frqBand='beta [13,30]';
        case 5
            frqBand='gamma [30,50]';
    end
    disp(['FrBand: ',frqBand])
    T_pvals=table(dComp,pTs(:,j),pWR(:,j),pA(:,j))
    T_EffSiz=table(dComp,efszTs(:,j),efszWR(:,j),efszA(:,j))
    pause

end
set(gca,'YScale','log')
set(gca,'FontSize',18)
ylabel('Time (s)')
 if(j==4)
     set(gca,'YLim',[0 3])
 elseif(j==5)
     set(gca,'YLim',[0 4.5])
 end