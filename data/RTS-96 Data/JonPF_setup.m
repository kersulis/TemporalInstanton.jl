function JonPF_setup
% bus_data = [num, type, P(pu), Q(pu), Vmag(pu), Vang(rad), Bshunt]
bus_data = [1  3  0  0  1  0  0;
            2  2  0  0  1  0  1;
            3  1  4  0  1  0  0;
            4  2 -5 -1  1  0  0];
        
% line_data = [from, to, R, X, B]
line_data = [1  2  0  0.20  0;
             1  4  0  0.05  0;
             2  3  0  0.08  0;
             2  4  0  0.10  0];

bus  = bus_data(:,1);
ty   = bus_data(:,2);                  
Pinj = bus_data(:,3);
Qinj = bus_data(:,4);
Vmag = bus_data(:,5);
Vang = bus_data(:,6);
bs   = 1i*bus_data(:,7);
nbus = length(ty);

p = line_data;
yy = 1./(p(:,3)+1i*p(:,4));
b = 1i*p(:,5)/2;
Y=sparse([p(:,1);p(:,2);p(:,1);p(:,2);p(:,1);p(:,2);bus],...  
         [p(:,2);p(:,1);p(:,1);p(:,2);p(:,1);p(:,2);bus],...
         [-yy;-yy;yy;yy;b;b;bs]);
G = real(Y);
B = imag(Y);

% Solve the power flow problem (note: you can ignore maxMis with ~)
[Vmag,Vang,Pcalc,Qcalc,maxMis,converged] = powerflow(Vmag,Vang,Pinj,Qinj,ty,G,B,nbus);

%(keyboard allows you to access the command line within the function. It 
% enables you to see the output of the powerflow problem.)
keyboard 
end