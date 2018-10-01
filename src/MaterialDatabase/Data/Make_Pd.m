function Pd = Make_Pd()
%% Al
load('Pd.mat');
Pd.Mat = 'Pd';
Pd.M = 106.42;
Pd.Z = 46;
Pd.Density = 12.02; %g/cm^3
Pd.Density = Pd.Density*10^-21/Pd.M*6.022*10^23;
Pd.NvTPP = 10;
Pd.NvSGS = 10;
Pd.Eg = 0;
Pd.Ep = 30.61; 
Pd.Ef = 6.2; 

Pd.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Pd.DIIMFP.E0)
    WernerData = load([cd '\W_in\' Pd.Mat num2str(Pd.DIIMFP.E0(i)) '.diimfp']);
    Pd.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Pd.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Pd.Elastic.x = zeros(numel(E0),1);
Pd.Elastic.l_el = zeros(numel(E0),1);
Pd.Elastic.l_tr = zeros(numel(E0),1);
Pd.Elastic.x = E0;
Pd.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Pd.Z,E0);
toc
Pd.DECS.x = data(1).x;
for i = 1:numel(E0)
    Pd.Elastic.l_el(i) = 1/data(i).sigma_el/Pd.Density;
    Pd.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Pd.Density;
    Pd.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Pd.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2'};
Pd.XPS.EB = [24357;3611;3337;3180;677;565;537;347;342;93;63;57;9;8];
Pd.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
%S(:,4)=S(:,4)*10^-5;
%S(:,16)=S(:,16)*10^-5;
for i = 1:numel(Pd.XPS.EB)
    Pd.XPS.PCS(:,i) = S(:,i+3*(i-1)); %*10^-7
    Pd.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Pd.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Pd.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end