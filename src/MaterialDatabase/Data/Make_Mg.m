function Mg = Make_Mg
%% Mg
load ('Mg.mat');

Mg.Mat     = 'Mg';
Mg.M       = 24.3050;
Mg.Z       = 12 ;
Mg.Density = 1.738; %g/cm^3
Mg.Density = Mg.Density*10^-21/Mg.M*6.022*10^23;
% Mg.Density = 43.0621 ;
Mg.NvTPP   = 2.0;
Mg.NvSGS   = 2.0;
Mg.Eg      = 0;
Mg.Ep      = 10.89;
Mg.Ef      = 7.1;

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Mg.Elastic.x = zeros(numel(E0),1);
Mg.Elastic.l_el = zeros(numel(E0),1);
Mg.Elastic.l_tr = zeros(numel(E0),1);
Mg.Elastic.x = E0;
Mg.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Mg.Z,E0);
toc
Mg.DECS.x = data(1).x;
for i = 1:numel(E0)
    Mg.Elastic.l_el(i) = 1/data(i).sigma_el/Mg.Density;
    Mg.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Mg.Density;
    Mg.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Mg.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';};
Mg.XPS.EB = [1308;92;54;54;7;];
Mg.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
for i = 1:numel(Mg.XPS.EB)
    Mg.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Mg.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Mg.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Mg.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end