function V = Make_V()
%% V
load('V.mat');
V.Mat = 'V'; V.M = 50.942;
V.Z = 23; %g/mol
V.Density = 6.11; %g/cm^-3
V.Density = V.Density/V.M*6.022*10^23*10^-21;
V.NvTPP = 5;
V.NvSGS = 5;
V.Eg = 0;
V.Ep = 22.30;
V.Ef = 6.4;

V.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(V.DIIMFP.E0)
    WernerData = load([cd '\W_in\' V.Mat num2str(V.DIIMFP.E0(i)) '.diimfp']);
    V.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         V.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
V.Elastic.x = zeros(numel(E0),1);
V.Elastic.l_el = zeros(numel(E0),1);
V.Elastic.l_tr = zeros(numel(E0),1);
V.Elastic.x = E0;
V.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(V.Z,E0);
toc
V.DECS.x = data(1).x;
for i = 1:numel(E0)
    V.Elastic.l_el(i) = 1/data(i).sigma_el/V.Density;
    V.Elastic.l_tr(i) = 1/data(i).sigma_tr1/V.Density;
    V.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

V.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'4S1/2';};
V.XPS.EB = [5470;633;525;518;72;44;43;8;7;];
V.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-4;
% S(:,16)=S(:,16)*10^-4;
% S(:,36)=S(:,36)*10^-4;
for i = 1:numel(V.XPS.EB)
    V.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    V.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    V.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    V.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end

end

% Кузнецова А.В.