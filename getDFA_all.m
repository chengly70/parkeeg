% calc DFA from processed data (dataFiltrd.mat); on bandpass & freq bands
% USING spliting (2 segments depending on subj, electrode(s), etc ) to get the best DFA fit
% relies on DFA_fun.m created by Martin Magris. UNcomment code to see figures

% --- !!! IF change which side/electrode, MUST update lines 192-194 (save mat file name) !!!! ---
%spltdfa=readcell('1dfaSplitsC3.xlsx');
%spltdfa=readcell('1dfaSplitsC4.xlsx');
spltdfa=readcell('1dfaSplitsMn.xlsx');

% --- !!! IF change which side/electrode, MUST update lines 192-194 (save mat file name) !!!! ---
%load dataEEG_Ctk2.mat
%load dataEEG_Ctk3.mat
load dataEEG_CtMn.mat

dt=deltaTVec(1);
nH=size(eegD.hlt,1);
nP=size(eegD.prkOff,1);

dtInd=1; %between 1 & 41, see deltaTVec

nStart=140;

% --- outputs ---
Fcdf=struct('hlt',[],'prkOff',[],'prkOn',[]); %the CDF (mean) of lin regr to fluct
Fcdf.hlt=cell(nH,5);
Fcdf.prkOff=cell(nP,5);
Fcdf.prkOn=cell(nP,5);
alps=struct('hlt',[],'prkOff',[],'prkOn',[]); %the CDF (mean) of lin regr to fluct
alps.hlt=zeros(nH,5);
alps.prkOff=zeros(nP,5);
alps.prkOn=zeros(nP,5);
xSclsHlt=cell(nH,5);
xSclsPrkOf=cell(nP,5);
xSclsPrkOn=cell(nP,5);
AvHlt=cell(nH,5);
AvPrkOf=cell(nP,5);
AvPrkOn=cell(nP,5);

for bpi=1:5

    hlt_d=eegD.hlt(:,dtInd,bpi);
    prkOff_d=eegD.prkOff(:,dtInd,bpi);
    prkOn_d=eegD.prkOn(:,dtInd,bpi);

    for j=1:nH

        tmpDat=hlt_d{j};
        % do it for all xScls
        xSclsHlt{j,bpi}.all=(nStart : 10 : round(length(tmpDat)/15))';
        [A,F]=DFA_fun(tmpDat,xSclsHlt{j,bpi}.all);
        Fcdf.hlt{j,bpi}=F; %save whole CDF
        if(spltdfa{j+1,bpi+1}==0)
            AvHlt{j,bpi}(1,:)=A;
            % figure
            % plot(log(xSclsHlt{j,bpi}.all),log(Fcdf.hlt{j,bpi}))
            % hold on
            % plot(log(xSclsHlt{j,bpi}.all),AvHlt{j,bpi}(1,1)*log(xSclsHlt{j,bpi}.all)+AvHlt{j,bpi}(1,2),'k')
        elseif(~ischar(spltdfa{j+1,bpi+1}) && isscalar(spltdfa{j+1,bpi+1})) %just 1 dividing pt
            xall=xSclsHlt{j,bpi}.all;
            dvdx=exp(spltdfa{j+1,bpi+1});
            [~,id]=min(abs(xall-dvdx)); %find closest index
            A=DFA_fun(tmpDat,xall(id+1:end));
            AvHlt{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id));
            AvHlt{j,bpi}(2,:)=A;
            xSclsHlt{j,bpi}.frst=xall(id+1:end); xSclsHlt{j,bpi}.secd=xall(1:id);
            % figure
            % plot(log(xSclsHlt{j,bpi}.all),log(Fcdf.hlt{j,bpi}))
            % hold on
            % plot(log(xSclsHlt{j,bpi}.frst),AvHlt{j,bpi}(1,1)*log(xSclsHlt{j,bpi}.frst)+AvHlt{j,bpi}(1,2),'k')
            % plot(log(xSclsHlt{j,bpi}.secd),AvHlt{j,bpi}(2,1)*log(xSclsHlt{j,bpi}.secd)+AvHlt{j,bpi}(2,2),'k')
        elseif(ischar(spltdfa{j+1,bpi+1})) %format is '<x, >y'
            xall=xSclsHlt{j,bpi}.all;
            leftBnd=str2double(extractBetween(spltdfa{j+1,bpi+1},'<',','));
            rightBnd=str2double(extractAfter(spltdfa{j+1,bpi+1},'>'));
            dvdx=exp(leftBnd);
            [~,id1]=min(abs(xall-dvdx)); %find closest index
            [~,id2]=min(abs(xall-exp(rightBnd))); %find closest index
            A=DFA_fun(tmpDat,xall(id2+1:end));
            AvHlt{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id1));
            AvHlt{j,bpi}(2,:)=A;
            xSclsHlt{j,bpi}.frst=xall(id2+1:end); xSclsHlt{j,bpi}.secd=xall(1:id1);
            %         figure
            %         plot(log(xSclsHlt{j,bpi}.all),log(Fcdf.hlt{j,bpi}))
            %         hold on
            %         plot(log(xSclsHlt{j,bpi}.frst),AvHlt{j,bpi}(1,1)*log(xSclsHlt{j,bpi}.frst)+AvHlt{j,bpi}(1,2),'k')
            %         plot(log(xSclsHlt{j,bpi}.secd),AvHlt{j,bpi}(2,1)*log(xSclsHlt{j,bpi}.secd)+AvHlt{j,bpi}(2,2),'k')
        end
        alps.hlt(j,bpi)=AvHlt{j,bpi}(1,1); %tail alpha only

    end

    %% repeat for Park Off
    for j=1:nP

        tmpDat=prkOff_d{j};
        % do it for all xScls
        xSclsPrkOf{j,bpi}.all=(nStart : 10 : round(length(tmpDat)/15))';
        [A,F]=DFA_fun(tmpDat,xSclsPrkOf{j,bpi}.all);
        Fcdf.prkOff{j,bpi}=F; %save whole CDF
        if(spltdfa{nH+j+1,bpi+1}==0)
            AvPrkOf{j,bpi}(1,:)=A;
            % figure
            % plot(log(xSclsPrkOf{j,bpi}.all),log(Fcdf.prkOff{j,bpi}))
            % hold on
            % plot(log(xSclsPrkOf{j,bpi}.all),AvPrkOf{j,bpi}(1,1)*log(xSclsPrkOf{j,bpi}.all)+AvPrkOf{j,bpi}(1,2),'k')
        elseif(~ischar(spltdfa{nH+j+1,bpi+1}) && isscalar(spltdfa{nH+j+1,bpi+1})) %just 1 dividing pt
            xall=xSclsPrkOf{j,bpi}.all;
            dvdx=exp(spltdfa{nH+j+1,bpi+1});
            [~,id]=min(abs(xall-dvdx)); %find closest index
            A=DFA_fun(tmpDat,xall(id+1:end));
            AvPrkOf{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id));
            AvPrkOf{j,bpi}(2,:)=A;
            xSclsPrkOf{j,bpi}.frst=xall(id+1:end);  xSclsPrkOf{j,bpi}.secd=xall(1:id);
            % figure
            % plot(log(xSclsPrkOf{j,bpi}.all),log(Fcdf.prkOff{j,bpi}))
            % hold on
            % plot(log(xSclsPrkOf{j,bpi}.frst),AvPrkOf{j,bpi}(1,1)*log(xSclsPrkOf{j,bpi}.frst)+AvPrkOf{j,bpi}(1,2),'k')
            % plot(log(xSclsPrkOf{j,bpi}.secd),AvPrkOf{j,bpi}(2,1)*log(xSclsPrkOf{j,bpi}.secd)+AvPrkOf{j,bpi}(2,2),'k')
        elseif(ischar(spltdfa{nH+j+1,bpi+1})) %format is '<x, >y'
            xall=xSclsPrkOf{j,bpi}.all;
            leftBnd=str2double(extractBetween(spltdfa{nH+j+1,bpi+1},'<',','));
            rightBnd=str2double(extractAfter(spltdfa{nH+j+1,bpi+1},'>'));
            dvdx=exp(leftBnd);
            [~,id1]=min(abs(xall-dvdx)); %find closest index
            [~,id2]=min(abs(xall-exp(rightBnd))); %find closest index
            A=DFA_fun(tmpDat,xall(id2+1:end));
            AvPrkOf{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id1));
            AvPrkOf{j,bpi}(2,:)=A;
            xSclsPrkOf{j,bpi}.frst=xall(id2+1:end); xSclsPrkOf{j,bpi}.secd=xall(1:id1);
            %         figure
            %         plot(log(xSclsPrkOf{j,bpi}.all),log(Fcdf.prkOff{j,bpi}))
            %         hold on
            %         plot(log(xSclsPrkOf{j,bpi}.frst),AvPrkOf{j,bpi}(1,1)*log(xSclsPrkOf{j,bpi}.frst)+AvPrkOf{j,bpi}(1,2),'k')
            %         plot(log(xSclsPrkOf{j,bpi}.secd),AvPrkOf{j,bpi}(2,1)*log(xSclsPrkOf{j,bpi}.secd)+AvPrkOf{j,bpi}(2,2),'k')
        end
        alps.prkOff(j,bpi)=AvPrkOf{j,bpi}(1,1); %tail alpha only
    end

    %% repeat for Park On
    for j=1:nP

        tmpDat=prkOn_d{j};
        % do it for all xScls
        xSclsPrkOn{j,bpi}.all=(nStart : 10 : round(length(tmpDat)/15))';
        [A,F]=DFA_fun(tmpDat,xSclsPrkOn{j,bpi}.all);
        Fcdf.prkOn{j,bpi}=F; %save whole CDF
        if(spltdfa{nH+nP+j+1,bpi+1}==0)
            AvPrkOn{j,bpi}(1,:)=A;
            % figure
            % plot(log(xSclsPrkOn{j,bpi}.all),log(Fcdf.prkOn{j,bpi}))
            % hold on
            % plot(log(xSclsPrkOn{j,bpi}.all),AvPrkOn{j,bpi}(1,1)*log(xSclsPrkOn{j,bpi}.all)+AvPrkOn{j,bpi}(1,2),'k')
        elseif(~ischar(spltdfa{nH+nP+j+1,bpi+1}) && isscalar(spltdfa{nH+nP+j+1,bpi+1})) %just 1 dividing pt
            xall=xSclsPrkOn{j,bpi}.all;
            dvdx=exp(spltdfa{nH+nP+j+1,bpi+1});
            [~,id]=min(abs(xall-dvdx)); %find closest index
            A=DFA_fun(tmpDat,xall(id+1:end));
            AvPrkOn{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id));
            AvPrkOn{j,bpi}(2,:)=A;
            xSclsPrkOn{j,bpi}.frst=xall(id+1:end);  xSclsPrkOn{j,bpi}.secd=xall(1:id);
            % figure
            % plot(log(xSclsPrkOn{j,bpi}.all),log(Fcdf.prkOn{j,bpi}))
            % hold on
            % plot(log(xSclsPrkOn{j,bpi}.frst),AvPrkOn{j,bpi}(1,1)*log(xSclsPrkOn{j,bpi}.frst)+AvPrkOn{j,bpi}(1,2),'k')
            % plot(log(xSclsPrkOn{j,bpi}.secd),AvPrkOn{j,bpi}(2,1)*log(xSclsPrkOn{j,bpi}.secd)+AvPrkOn{j,bpi}(2,2),'k')
        elseif(ischar(spltdfa{nH+nP+j+1,bpi+1})) %format is '<x, >y'
            xall=xSclsPrkOn{j,bpi}.all;
            leftBnd=str2double(extractBetween(spltdfa{nH+nP+j+1,bpi+1},'<',','));
            rightBnd=str2double(extractAfter(spltdfa{nH+nP+j+1,bpi+1},'>'));
            dvdx=exp(leftBnd);
            [~,id1]=min(abs(xall-dvdx)); %find closest index
            [~,id2]=min(abs(xall-exp(rightBnd))); %find closest index
            A=DFA_fun(tmpDat,xall(id2+1:end));
            AvPrkOn{j,bpi}(1,:)=A;
            A=DFA_fun(tmpDat,xall(1:id1));
            AvPrkOn{j,bpi}(2,:)=A;
            xSclsPrkOn{j,bpi}.frst=xall(id2+1:end); xSclsPrkOn{j,bpi}.secd=xall(1:id1);
            %         figure
            %         plot(log(xSclsPrkOn{j,bpi}.all),log(Fcdf.prkOn{j,bpi}))
            %         hold on
            %         plot(log(xSclsPrkOn{j,bpi}.frst),AvPrkOn{j,bpi}(1,1)*log(xSclsPrkOn{j,bpi}.frst)+AvPrkOn{j,bpi}(1,2),'k')
            %         plot(log(xSclsPrkOn{j,bpi}.secd),AvPrkOn{j,bpi}(2,1)*log(xSclsPrkOn{j,bpi}.secd)+AvPrkOn{j,bpi}(2,2),'k')
        end
        alps.prkOn(j,bpi)=AvPrkOn{j,bpi}(1,1); %tail alpha only
    end

    %save dDFA_all_Ctk2 Fcdf xScls* Av* alps
    %save dDFA_all_Ctk3 Fcdf xScls* Av* alps
    save dDFA_all_CtMn Fcdf xScls* Av* alps
end