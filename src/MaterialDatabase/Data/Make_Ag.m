function Ag = Make_Ag()
%% Ag
load('Ag.mat');
Ag.Mat = 'Ag';
Ag.M = 107.8682;
Ag.Z = 47;
Ag.Density = 10.5; % g/cm^3
Ag.Density = Ag.Density*10^-21/Ag.M*6.022*10^23;
Ag.NvTPP = 11;
Ag.NvSGS = 11;
Ag.Eg = 0;
Ag.Ep = 29.8;
Ag.Ef = 7.2;

Ag.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Ag.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Ag.Mat num2str(Ag.DIIMFP.E0(i)) '.diimfp']);
    Ag.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Ag.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Ag.Elastic.x = zeros(numel(E0),1);
Ag.Elastic.l_el = zeros(numel(E0),1);
Ag.Elastic.l_tr = zeros(numel(E0),1);
Ag.Elastic.x = E0;
Ag.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Ag.Z,E0);
toc
Ag.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ag.Elastic.l_el(i) = 1/data(i).sigma_el/Ag.Density;
    Ag.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ag.Density;
    Ag.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Ag.XPS.Shells = {'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'5S1/2';};
Ag.XPS.EB = [724;608;577;379;373;101;69;63;11;10;8;];
Ag.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
for i = 1:numel(Ag.XPS.EB)
    Ag.XPS.PCS(:,i) = A(:,i+3*(i-1))*10^-7;
    Ag.XPS.betta(:,i) = A(:,i+3*(i-1)+1);
    Ag.XPS.gamma(:,i) = A(:,i+3*(i-1)+2);
    Ag.XPS.delta(:,i) = A(:,i+3*(i-1)+3);
end