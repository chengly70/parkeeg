function signl = getContn(contM,windows,deltaT)
% INPUT: contM= Lt x 2 matrix (1stCol=time vect dt, 2ndCol=signal),
% windows= N x 2 matrix (each for is time interval (w(i,1) , w(i,2))
% binning in windows with deltaT time bin
% i.e., downsampling

timV=contM(:,1); 
cont=contM(:,2);
dt=mean(diff(timV)); %!!! ASSUMING timV is equally spaced!! If not, change to timV(2)-timV(1) assuming first 2 pts suffice for dt

cont =cont'; %make it a row vector (Sam)

signl=[];
nmSteps=round(deltaT/dt);

for i=1:size(windows,1)
    idTw = (timV>=windows(i,1) & timV<windows(i,2));
    if( deltaT <= dt || (deltaT-dt<1e-6) ) %use smallest possible time bin, deltaT=dt
        signl=[signl cont(idTw)];
    else %downsample b/c deltaT>dt, use avg
        tmp1=cont(idTw);
        nmDT=floor(length(tmp1)/nmSteps);
        tmp2=zeros(1,nmDT);
        if(nmSteps*nmDT == length(tmp1)) %deltaT divides window length
            tmp2=mean(reshape(tmp1,nmSteps,nmDT)); %take average
        else %do all but last window (do last window separately)
            tmp2(1:nmDT-1)=mean(reshape(tmp1(1:nmSteps*(nmDT-1)),nmSteps,nmDT-1));
            tmp2(nmDT) = mean(tmp1(nmSteps*(nmDT-1)+1:end));
        end
        signl=[signl tmp2];
    end
end