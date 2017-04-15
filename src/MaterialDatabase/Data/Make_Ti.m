function Ti = Make_Ti()
%% Ti
load('Ti.mat');
Ti.Mat = 'Ti';Ti.M = 47.867;%g/mol
Ti.Z = 22; 
Ti.Density = 4.51; %g/cm^-3
Ti.Density = Ti.Density/Ti.M*6.022*10^23*10^-21;
Ti.NvTPP = 4;
Ti.NvSGS = 4;
Ti.Eg = 0;
Ti.Ep = 17.68;
Ti.Ef = 6.0;

% Ti.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(Ti.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' Ti.Mat num2str(Ti.DIIMFP.E0(i)) '.diimfp']);
%     Ti.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%          Ti.DIIMFP.x = WernerData(:,1);
%     end
% end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Ti.Elastic.x = zeros(numel(E0),1);
Ti.Elastic.l_el = zeros(numel(E0),1);
Ti.Elastic.l_tr = zeros(numel(E0),1);
Ti.Elastic.x = E0;
Ti.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Ti.Z,E0);
toc
Ti.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ti.Elastic.l_el(i) = 1/data(i).sigma_el/Ti.Density;
    Ti.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ti.Density;
    Ti.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Ti.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'4S1/2';};
Ti.XPS.EB = [4970;567;465;459;64;39;38;8;7;];
Ti.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-4;
% S(:,16)=S(:,16)*10^-4;
% S(:,36)=S(:,36)*10^-4;
for i = 1:numel(Ti.XPS.EB)
    Ti.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Ti.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Ti.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Ti.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end

end

% Кузнецова А.В.