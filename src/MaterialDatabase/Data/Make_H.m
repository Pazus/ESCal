function H = Make_H()
load('H.mat');
H.Mat = 'H';
H.M = 1.00797;
H.Z = 1;
H.Density = 0.071; %g/cm^3
H.Density = H.Density*10^-21/H.M*6.022*10^23;
H.NvTPP = 1;
H.NvSGS = 1;
H.Eg = 0;
H.Ep = NaN; %free-electron plasmon energy
H.Ef = NaN; %Fermi energy

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
H.Elastic.x = zeros(numel(E0),1);
H.Elastic.l_el = zeros(numel(E0),1);
H.Elastic.l_tr = zeros(numel(E0),1);
H.Elastic.x = E0;
H.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(H.Z,E0);
toc
H.DECS.x = data(1).x;
for i = 1:numel(E0)
    H.Elastic.l_el(i) = 1/data(i).sigma_el/H.Density;
    H.Elastic.l_tr(i) = 1/data(i).sigma_tr1/H.Density;
    H.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

% H.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(H.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' H.Mat num2str(H.DIIMFP.E0(i)) '.diimfp']);
%     H.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%         H.DIIMFP.x = WernerData(:,1);
%     end
% end

H.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000];
H.XPS.Shells = {'1S1/2'};
H.XPS.EB = 14;

for i = 1:numel(H.XPS.EB)
    H.XPS.PCS(:,i) = Q(:,i+3*(i-1))*10^-7;
    H.XPS.betta(:,i) = Q(:,i+3*(i-1)+1);
    H.XPS.gamma(:,i) = Q(:,i+3*(i-1)+2);
    H.XPS.delta(:,i) = Q(:,i+3*(i-1)+3);
end



