% Convert ARPA-E to MATPOWER
% 3-14-14
% Jonas Kersulis

% This script relies on an excel file and two MATLAB scripts developed by
% Jennifer Felder:
%  - UWCase3_Dec2013.xlsx contains all case-specific data
%  - CaseData.m extracts case data from the above Excel file
%  - NetworkData.m is called by CaseData.m and imports network data
% The output of this script is an MPC struct compatible with MATPOWER's
% power flow algorithms called 'caseRTS96'.

% Import all relevant Excel data into MATLAB variables. This can take
% upwards of 20s:
CaseData % This script calls another script called "NetworkData"

% Clear variables containing unnecessary data:
clear Bmax BusNum DemData GenMaxZ GenPar RUp RDwn StorData ...
      StorID Stor_B WB effC effD hours i q rStor_L1 rcMax rdMax t

% Demand and generator data depend on the hour:
hour = 7;

caseRTS96.version = '2';

caseRTS96.baseMVA = 100;

nodes = size(nodedata,1);

% Bus Data has 13 columns:
bus_i = nodedata(:,1);
type = nodedata(:,2);

Wind(1:nodes,1) = 0;
for i = 1:size(WindID,1)
    for j = 1:size(bus_i,1)
        if WindID(i) == bus_i(j)
            index = find(WindData(1,:) == WindID(i));
           Wind(j) = sum(WindData(hour + 1,index));
        end
    end
end

Pd = Pd(:,hour + 1).*Sb - Wind; % Reduce load by wind
Qd = nodedata(:,4).*Sb;
Gs = nodedata(:,5);
Bs = nodedata(:,6);
area = nodedata(:,7);
Vm(1:nodes,1) = 1;
Va(1:nodes,1) = 0;
baseKV = nodedata(:,8);
zone = nodedata(:,9);
Vmax(1:nodes,1) = 1.1;
Vmin(1:nodes,1) = 0.9;

% Use same slack bus as Jenny:
%type(13) = 2;
%type(21) = 3;

% Assemble case data:
caseRTS96.bus = [bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin];

% Gen Data has 21 columns:
gens = size(gendata,1);

bus = gendata(:,1);
Pg = gendata(:,3).*Sb;
Qg = gendata(:,4).*Sb;
Qmax = gendata(:,5).*Sb;
Qmin = gendata(:,6).*Sb;
Vg = gendata(:,7);
mBase(1:gens,1) = Sb;
status = [GenMax(1:14,hour); 0; GenMax(15:47,hour); 0; ...
          GenMax(48:80,hour); 0; GenMax(81:96,hour)]; % Insert Sync Cond
Pmax = Sb.*[PMax(1:14); 0.003; 
            PMax(15:47); 0.003;
            PMax(48:80); 0.003;
            PMax(81:96)];
Pmin = Sb.*[PMin(1:14); 0; PMin(15:47); 0; PMax(48:80); 0; PMin(81:96)];
Pc1(1:gens,1) = 0;
Pc2 = Pc1;
Qc1min = Pc1;
Qc1max = Pc1;
Qc2min = Pc1;
Qc2max = Pc1;
ramp_agc = Pc1;
ramp_10 = Pc1;
ramp_30 = Pc1;
ramp_q = Pc1;
apf = Pc1;

caseRTS96.gen = [bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 ...
                Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ...
                ramp_q apf];

% Branch Data has 13 columns:
branches = size(branchdata,1);

fbus = branchdata(:,1);
tbus = branchdata(:,2);
r = branchdata(:,7);
x = branchdata(:,8);
b = branchdata(:,9);
rateA = 0.8.*Sb.*branchdata(:,10);  % Reduce to 80 %
rateB = 0.8.*Sb.*branchdata(:,11);  % Reduce to 80 %
rateC = 0.8.*Sb.*branchdata(:,12);  % Reduce to 80 %
ratio(1:branches,1) = 0;
angle = ratio;
status(1:branches,1) = 1;
angmin(1:branches,1) = -360;
angmax = -angmin;

caseRTS96.branch = [fbus tbus r x b rateA rateB rateC ratio angle ...
                    status angmin angmax];

% Gen Cost Data has 5 columns:
model(1:gens,1) = 1;
startup = CPL(:,2); 
startup = [startup(1:14); 0; 
           startup(15:47); 0; 
           startup(48:80); 0; 
           startup(81:96)];
shutdown(1:gens,1) = 0;
n(1:gens,1) = 3;
x1 = CPL(:,6); 
x1 = [x1(1:14); 0.001;
      x1(15:47); 0.001;
      x1(48:80); 0.001;
      x1(81:96)];
y1 = CPL(:,3); 
y1 = [y1(1:14); 1; 
      y1(15:47); 2; 
      y1(48:80); 3; 
      y1(81:96)];
x2 = 2.*x1;
x3 = 3.*x1;
% x2 = 2.*CPL(:,7); 
% x2 = [x2(1:14); 0; 
%       x2(15:47); 0; 
%       x2(48:80); 0; 
%       x2(81:96)];
y2 = CPL(:,4); 
y2 = [y2(1:14); 1; 
      y2(15:47); 2; 
      y2(48:80); 3; 
      y2(81:96)];
% x3 = 3.*CPL(:,8); 
% x3 = [x3(1:14); 0; 
%       x3(15:47); 0; 
%       x3(48:80); 0; 
%       x3(81:96)];
y3 = CPL(:,5); 
y3 = [y3(1:14); 1; 
      y3(15:47); 2; 
      y3(48:80); 3; 
      y3(81:96)];

% Finish the case data:
caseRTS96.gencost = [model startup shutdown n x1 y1 x2 y2 x3 y3];

%% Next, use Jon's power flow:
% Vmag = caseRTS96.bus(:,8);
% Vang = caseRTS96.bus(:,9);
% Pinj = caseRTS96.bus(:,8);
% 
% [Vmag,Vang,Pcalc,Qcalc,maxMis,converged] = powerflow(Vmag,Vang,Pinj,Qinj,ty,G,B,nbus);
% 
% %(keyboard allows you to access the command line within the function. It 
% % enables you to see the output of the powerflow problem.)
% keyboard 