%script to get which side Parkinson's severe (0=both, 1=left, 2=right)

cd dataUCSD/
allF=dir('*.bdf'); %get all .bdf file names
fileHealth=dir('sub-hc*.bdf');          
fileParkOffMed=dir('sub-pd*off*.bdf');  
fileParkOnMed=dir('sub-pd*on*.bdf');
cd ..

preF='dataUCSD/';
nH=length(fileHealth); %number of healthy subjects
nP=length(fileParkOffMed); % #Parkinsons; same as using *On* Medication

% --- OUTPUTS ---
whSide=struct('hlt',[],'prkOff',[],'prkOn',[]);

%get data in table
partTab=readtable('participants.xlsx');

cl_nam=1; cl_sdOff=8; cl_sdOn=9; %columns for various data pieces

%manually enter row number for health [CL CHECKED that 
% partTab{j,cl_nam} matches with fileHealth(j).name AND filePark[Off/ON]Med(j).name]
rown_Hlt=[10; 17; 1; 19; 20; 23; 24; 27; 2; 28; 29; 30; 31; 4; 7; 8];
rown_Park=[12; 11; 13; 14; 15; 16; 18; 21; 22; 25; 26; 3; 5; 6; 9];

whSide.hlt=zeros(nH,1);
for j=1:nH
    if( strcmp(partTab{rown_Hlt(j),cl_sdOff}, 'n/a')==1 )
        tmp=0;
    elseif( strcmp(partTab{rown_Hlt(j),cl_sdOff}, 'L')==1 )
        tmp=1;
    elseif( strcmp(partTab{rown_Hlt(j),cl_sdOff}, 'R')==1 )
        tmp=2;
    end
    whSide.hlt(j)=tmp;
end

whSide.prkOff=zeros(nP,1);
whSide.prkOn=zeros(nP,1);
for j=1:nP
    if( strcmp(partTab{rown_Park(j),cl_sdOff}, 'n/a')==1 )
        tmp=0;
    elseif( strcmp(partTab{rown_Park(j),cl_sdOff}, 'L')==1 )
        tmp=1;
    elseif( strcmp(partTab{rown_Park(j),cl_sdOff}, 'R')==1 )
        tmp=2;
    end
    whSide.prkOff(j)=tmp;

    if( strcmp(partTab{rown_Park(j),cl_sdOn}, 'n/a')==1 )
        tmp=0;
    elseif( strcmp(partTab{rown_Park(j),cl_sdOn}, 'L')==1 )
        tmp=1;
    elseif( strcmp(partTab{rown_Park(j),cl_sdOn}, 'R')==1 )
        tmp=2;
    end
    whSide.prkOn(j)=tmp;
end

save('dSeveritySide.mat','whSide')

% % to  manually check spreadsheet, it matches!
% for j=1:nH
%     whSide.hlt(j)
%     rown_Hlt(j)
%     pause
% end
% for j=1:nP
%     [whSide.prkOff(j) whSide.prkOn(j)]
%     rown_Park(j)
%     pause
% end