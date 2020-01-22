function Hf = Make_Hf()
%% Hf
% load('Hf.mat');
Hf.Mat = 'Hf';
Hf.M = 178.49;
Hf.Z = 72;
Hf.Density = 13.31; %g/cm^3
Hf.Density = Hf.Density*10^-21/Hf.M*6.022*10^23;
Hf.NvTPP = 4;
Hf.NvSGS = 18;
Hf.Eg = 0;
Hf.Ep = 15.73;
Hf.Ef = 7.9;

% Hf.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(Hf.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' Hf.Mat num2str(Hf.DIIMFP.E0(i)) '.diimfp']);
%     Hf.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%          Hf.DIIMFP.x = WernerData(:,1);
%     end
% end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Hf.Elastic.x = zeros(numel(E0),1);
Hf.Elastic.l_el = zeros(numel(E0),1);
Hf.Elastic.l_tr = zeros(numel(E0),1);
Hf.Elastic.x = E0;
Hf.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Hf.Z,E0);
toc
Hf.DECS.x = data(1).x;
for i = 1:numel(E0)
    Hf.Elastic.l_el(i) = 1/data(i).sigma_el/Hf.Density;
    Hf.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Hf.Density;
    Hf.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

% Hf.XPS.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
% Hf.XPS.EB = [121;77;77;10;6;];
% Hf.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-5;
% S(:,16)=S(:,16)*10^-5;
% for i = 1:numel(Hf.XPS.EB)
%     Hf.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
%     Hf.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
%     Hf.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
%     Hf.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
% end