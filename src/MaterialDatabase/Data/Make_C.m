function C = Make_C()
load('C.mat');
C.Mat = 'C';
C.M = 12.011; %g/mol
C.Z = 6; 
C.Density = 2.25; %g/cm^-3
C.Density = C.Density/C.M*6.022*10^23*10^-21;
C.NvTPP = 4;
C.NvSGS = 4;
C.Eg = 0;
C.Ep      = 24.93;
C.Ef      = 20.4;

% C.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(C.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' C.Mat num2str(C.DIIMFP.E0(i)) '.diimfp']);
%     C.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%          C.DIIMFP.x = WernerData(:,1);
%     end
% end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
C.Elastic.x = zeros(numel(E0),1);
C.Elastic.l_el = zeros(numel(E0),1);
C.Elastic.l_tr = zeros(numel(E0),1);
C.Elastic.x = E0;
C.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(C.Z,E0);
toc
C.DECS.x = data(1).x;
for i = 1:numel(E0)
    C.Elastic.l_el(i) = 1/data(i).sigma_el/C.Density;
    C.Elastic.l_tr(i) = 1/data(i).sigma_tr1/C.Density;
    C.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

C.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';};
C.XPS.EB = [288; 16; 11;];
C.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];

for i = 1:numel(C.XPS.EB)
    C.XPS.PCS(:,i) = C_gr(:,i+3*(i-1))*10^-7;
    C.XPS.betta(:,i) = C_gr(:,i+3*(i-1)+1);
    C.XPS.gamma(:,i) = C_gr(:,i+3*(i-1)+2);
    C.XPS.delta(:,i) = C_gr(:,i+3*(i-1)+3);
end
