%Define Case Data (From UW):
% clear;
%NOTE: Before Running a New Case for the RTS System, Make Sure The Files Below Are Updated
%Accordingly:
%UWCase.xlsx
%If you are using this on another network other than the RTS System, then
%the following network-dependent files must also be updated:
%NetworkData.m

WB = 'UWCase3_Dec2013.xlsx';

%Get Network Toplogy and Generator Data:
NetworkData
q=find(gendata(:,1)==107);
gendata(q(3),1)=108;
q=find(gendata(:,1)==207);
gendata(q(2:3),1)=208;
q=find(gendata(:,1)==307);
gendata(q(2:3),1)=308;


%Get Demand Data:
% DemData=xlsread('UWCase.xlsx','Demand');
% DemData=xlsread('UWCase2.xlsx','Demand');
DemData=xlsread(WB,'Demand');
% DemData=xlsread('UWCase3_Dec2013_R1.xlsx','Demand');

hours=DemData(2:length(DemData(:,1)),1);
% Pd(1:73,1:36) = 0;
for i=1:length(hours)
    Pd(:,i)=DemData(i+1,:)/Sb;
end
Pd(1,:)=[];
Pd=[DemData(1,2:length(DemData(1,:))).' Pd];

%Get Generator Parmeter Data:
% GenPar=xlsread('UWCase.xlsx','GenList');
% GenPar=xlsread('UWCase2.xlsx','GenList');
GenPar=xlsread(WB,'GenList');
% GenPar=xlsread('UWCase3_Dec2013_R1.xlsx','GenList');

BusNum=GenPar(:,1);
GenBus=GenPar(:,2);
PMax=GenPar(:,3)/Sb;
CPL(:,1)=BusNum;
CPL(:,2)=GenPar(:,4);
RUp=GenPar(:,5)/Sb;
RDwn=GenPar(:,6)/Sb;
PMin=GenPar(:,7)/Sb;
P0=GenPar(:,8)/Sb;

CPL(:,3:5)=GenPar(:,9:11);
CPL(:,6:8)=GenPar(:,12:14);

% NOTE: This block turns a linear cost curve into a quadratic form. It
% requires additional functions from Jenny.
% [Qcoeff, Lcoeff, Ccoeff]=GetQCC(CPL);
% 
% CstQ=[Qcoeff*100^2];
% CstL=[Lcoeff*100]; 
% CstN=[Ccoeff];
% C_PL=CPL(:,3:5);
% Pgmax=[CPL(:,6)/Sb (CPL(:,6)+CPL(:,7))/Sb (CPL(:,6)+CPL(:,7)+CPL(:,8))/Sb];


%Get Generator Commitment Data:
% GenCom=xlsread('UWCase.xlsx','GenCom');
% GenCom=xlsread('UWCase2.xlsx','GenCom');
GenCom=xlsread(WB,'GenCom');
% GenCom=xlsread('UWCase3_Dec2013_R1.xlsx','GenCom');

Pgmin(1:96,1:length(hours)) = 0;
PgminZ = Pgmin;
GenMaxZ = Pgmin;
GenMax = Pgmin;
for t=1:length(hours)
    for i=1:length(GenBus)
        if(GenCom(t,i+1)==0)
            Pgmin(i,t)=0;
            PgminZ(i,t)=PMin(i);
            GenMaxZ(i,t)=PMax(i);
            GenMax(i,t)=0;
        else
            Pgmin(i,t)=PMin(i);
            PgminZ(i,t)=PMin(i);
            GenMax(i,t)=PMax(i);
            GenMaxZ(i,t)=PMax(i);
        end
    end
end

%Get Storage Data:
% StorData=xlsread('UWCase.xlsx','Stor');
% StorData=xlsread('UWCase2.xlsx','Stor');
StorData=xlsread(WB,'Stor');
% StorData=xlsread('UWCase3_Dec2013_R1.xlsx','Stor');

StorID=StorData(1,2:length(StorData(1,:))).';
effC=StorData(2,2);
effD=StorData(3,2);
rcMax=1/Sb*StorData(4,2:length(StorData(1,:))).';
rdMax=1/Sb*StorData(5,2:length(StorData(1,:))).';
Bmax=1/Sb*StorData(6,2:length(StorData(1,:))).';

Stor_B(1:8,1:length(hours)) = 0;
for i=1:length(hours)
    Stor_B(:,i)=StorData(i+6,2:length(StorID)+1)/Sb;
end

%Get Storage Data:
% StorData=xlsread('UWCase.xlsx','rStor');
% StorData=xlsread('UWCase2.xlsx','rStor');
StorData=xlsread(WB,'rStor');
% StorData=xlsread('UWCase3_Dec2013_R1.xlsx','rStor');

rStor_L1(1:8,1:length(hours)) = 0;
for i=1:length(hours)
    rStor_L1(:,i)=StorData(i+1,2:length(StorID)+1);
end

%Get Wind Data:
% WindData=xlsread('UWCase.xlsx','Wind');
% WindData=xlsread('UWCase2.xlsx','Wind');
WindData=xlsread(WB,'Wind');
% WindData=xlsread('UWCase3_Dec2013_R1.xlsx','Wind');

WindID=WindData(1,2:length(WindData(1,:))).';
WindMax(1:19,1:length(hours)) = 0;
for i=1:length(hours)
    WindMax(:,i)=WindData(i+1,2:length(WindID)+1)/Sb;
end

%Get Wind Curtailment Data:
% WindDataC=xlsread('UWCase.xlsx','WindCur');
% WindDataC=xlsread('UWCase2.xlsx','WindCur');
WindDataC=xlsread(WB,'WindCur');
% WindDataC=xlsread('UWCase3_Dec2013_R1.xlsx','WindCur');

WindCur(1:19,1:length(hours)) = 0;
for i=1:length(hours)
    WindCur(:,i)=WindDataC(i+1,2:length(WindID)+1);
end

WindDataB=xlsread(WB,'WindBnds');
% WindDataB=xlsread('UWCase3_Dec2013_R1.xlsx','WindBnds');

WindBndU(1:19,1:length(hours)) = 0;
WindBndL = WindBndU;
for i=1:length(hours)
    WindBndU(:,i)=WindDataB((i-1)*length(WindID)+1:i*length(WindID),1)/Sb;
    WindBndL(:,i)=WindDataB((i-1)*length(WindID)+1:i*length(WindID),2)/Sb;
end

