function Ni = Make_Ni()
%% Cu
load('Ni.mat');
Ni.Mat = 'Ni';
Ni.M = 58.6934;
Ni.Z = 28; %g/mol
Ni.Density = 8.908; %g/cm^-3
Ni.Density = Ni.Density/Ni.M*6.022*10^23*10^-21;
Ni.NvTPP = 10;
Ni.NvSGS = 10;
Ni.Eg = 0;
Ni.Ep = 35.47;
Ni.Ef = 9.1;

Ni.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Ni.DIIMFP.E0)
    WernerData = load([cd '\W_in\' Ni.Mat num2str(Ni.DIIMFP.E0(i)) '.diimfp']);
    Ni.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Ni.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Ni.Elastic.x = zeros(numel(E0),1);
Ni.Elastic.l_el = zeros(numel(E0),1);
Ni.Elastic.l_tr = zeros(numel(E0),1);
Ni.Elastic.x = E0;
Ni.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Ni.Z,E0);
toc
Ni.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ni.Elastic.l_el(i) = 1/data(i).sigma_el/Ni.Density;
    Ni.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ni.Density;
    Ni.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Ni.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Ni.XPS.EB = [8338;1015;877;860;117;75;73;10;10;8;];
Ni.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-4;
% S(:,16)=S(:,16)*10^-4;
% S(:,36)=S(:,36)*10^-4;
for i = 1:numel(Ni.XPS.EB)
    Ni.XPS.PCS(:,i) = S(1+(i-1)*4,:)'*10^-7;
    Ni.XPS.betta(:,i) = S(2+(i-1)*4,:);
    Ni.XPS.gamma(:,i) = S(3+(i-1)*4,:);
    Ni.XPS.delta(:,i) = S(4+(i-1)*4,:);
end