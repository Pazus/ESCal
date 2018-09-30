function Cu = Make_Cu()
%% Cu
load('Cu.mat');
Cu.Mat = 'Cu';
Cu.M = 63.546;
Cu.Z = 29; %g/mol
Cu.Density = 8.92; %g/cm^-3
Cu.Density = Cu.Density/Cu.M*6.022*10^23*10^-21;
Cu.NvTPP = 11;
Cu.NvSGS = 11;
Cu.Eg = 0;
Cu.Ep = 35.87;
Cu.Ef = 8.7;

Cu.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Cu.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Cu.Mat num2str(Cu.DIIMFP.E0(i)) '.diimfp']);
    Cu.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Cu.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Cu.Elastic.x = zeros(numel(E0),1);
Cu.Elastic.l_el = zeros(numel(E0),1);
Cu.Elastic.l_tr = zeros(numel(E0),1);
Cu.Elastic.x = E0;
Cu.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Cu.Z,E0);
toc
Cu.DECS.x = data(1).x;
for i = 1:numel(E0)
    Cu.Elastic.l_el(i) = 1/data(i).sigma_el/Cu.Density;
    Cu.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Cu.Density;
    Cu.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Cu.XPS.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Cu.XPS.EB = [1103;958;938;127;82;80;11;10;8;];
Cu.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
S(:,4)=S(:,4)*10^-4;
S(:,16)=S(:,16)*10^-4;
S(:,36)=S(:,36)*10^-4;
for i = 1:numel(Cu.XPS.EB)
    Cu.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Cu.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Cu.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Cu.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end