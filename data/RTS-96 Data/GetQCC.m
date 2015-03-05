function [Qcoeff Lcoeff Ccoeff]=GetQCC(CPL)

CPL(:,7)=CPL(:,7)+CPL(:,6);
CPL(:,8)=CPL(:,7)+CPL(:,8);

% Temp=CPL(12:14,:);
% CPL(12:14,:)=[];
% CPL=[CPL; Temp];

for i=1:length(CPL(:,1))
    P{i}=1:1:CPL(i,8);
end

for k=1:length(P)
    Cst=zeros(length(P{k}),1);
    Pow=P{k};
    for i=1:length(P{k})
        if(Pow(i)<=CPL(k,6))
            Cst(i,1)=CPL(k,2)+CPL(k,3)*Pow(i);
        elseif(Pow(i)>CPL(k,6)&&Pow(i)<=CPL(k,7))
            Cst(i,1)=CPL(k,2)+CPL(k,3)*CPL(k,6)+CPL(k,4)*(Pow(i)-CPL(k,6));
        elseif(Pow(i)>CPL(k,7))
            Cst(i,1)=CPL(k,2)+CPL(k,4)*CPL(k,7)+CPL(k,5)*(Pow(i)-CPL(k,7));
        end
    end
    Cst2{k}=Cst;
end

% plot(P{1},Cst2{1},'r');
% hold on
% plot(P{3},Cst2{3},'b');
% plot(P{9},Cst2{9},'g');

Qcoeff=zeros(length(P),1);
Lcoeff=zeros(length(P),1);
Ccoeff=zeros(length(P),1);
rse=zeros(length(P),1);

for i=1:length(P)

    Pow=P{i};
    Cst=Cst2{i};
    [pfit, gof] = fit(Pow.', Cst,  'poly2' );
    Qcoeff(i,1)=pfit.p1;
    Lcoeff(i,1)=pfit.p2;
    Ccoeff(i,1)=pfit.p3;
    rse(i,1)=gof.rsquare;
end
% plot(P{1},Cst2{1},'b')
% hold on
% plot(p1fit,'r')


end
