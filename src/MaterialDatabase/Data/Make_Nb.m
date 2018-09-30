function Nb = Make_Nb()
%% Nb
load('Nb.mat');
Nb.Mat = 'Nb';
Nb.M = 92.90638;
Nb.Z = 41; %g/mol
Nb.Density = 8.57; %g/cm^-3
Nb.Density = Nb.Density/Nb.M*6.022*10^23*10^-21;
Nb.NvTPP = 5;
Nb.NvSGS = 5;
Nb.Eg = 0;
Nb.Ep = 19.56;
Nb.Ef = 5.3;

Nb.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Nb.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Nb.Mat num2str(Nb.DIIMFP.E0(i)) '.diimfp']);
    Nb.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Nb.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Nb.Elastic.x = zeros(numel(E0),1);
Nb.Elastic.l_el = zeros(numel(E0),1);
Nb.Elastic.l_tr = zeros(numel(E0),1);
Nb.Elastic.x = E0;
Nb.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Nb.Z,E0);
toc
Nb.DECS.x = data(1).x;
for i = 1:numel(E0)
    Nb.Elastic.l_el(i) = 1/data(i).sigma_el/Nb.Density;
    Nb.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Nb.Density;
    Nb.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Nb.XPS.Shells = {'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'5S1/2';};
Nb.XPS.EB = [473;384;367;212;209;62;40;38;7;7;];
Nb.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
S(:,4)=S(:,4)*10^-4;
S(:,24)=S(:,24)*10^-4;
S(:,40)=S(:,40)*10^-4;
for i = 1:numel(Nb.XPS.EB)
    Nb.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Nb.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Nb.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Nb.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end

% Nb.DECS.x = Nb.DECS.x /180*pi;

% Nb.DECS.y = normalize(Nb.DECS.x,Nb.DECS.y);
%Nb.DIIMFP.y = normalize(Nb.DIIMFP.x,Nb.DIIMFP.y);

% MaterialData.Nb = Nb;