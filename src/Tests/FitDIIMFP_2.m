clear all;

%% Geometry settings
Geometry.theta0=0;
Geometry.theta=45;
Geometry.phi=0;

%% Parameters of the spectrum view
Element    =  'Al';
E0 = 2000;
Delta_max  = 150; %energy range value to calculate      
delta_CS   = 0:0.25:Delta_max; h=0.25; % h - step for DIIMFP                               

Thickness_s = 0.45; %thickness of the surface layer 
sigma_gauss = 0.38; %  eV
MKgraf      = 2.5e-4; % to fit the experimental scale to be the same with the calculation
 
%% Const for I_ion (ionization process)
Work_func    = 4.1; % a table value
% J_ion - the ionization potential (table values)
    J_ion      (1)   = 80  -Work_func;  a(1)   = 0;
    Mlambda_ion(1) = 0.0; % fit parameter
    J_ion      (2)   = 81  -Work_func;  a(2)   = 0;
    Mlambda_ion(2) = 0.0; % fit parameter
    J_ion      (3)   = 126 -Work_func;  a(3)   = 0;
    Mlambda_ion(3) = 0.15; % fit parameter
    J_ion      (4)   = 1569-Work_func;  a(4)   = 0;
    Mlambda_ion(4) = 0; % fit parameter

%% Const for Ind_bulk 

osc.A = 200;   
osc.G = 0.5;  
osc.Om = 14.8; 
osc.alpha = 1; %alpha is a constant between 1 (metals) and 0 (insulators)
osc.u = 12;
FermiEnergy = 11.2;
osc.model = 'Lindhard';

%NDIIMFP = nrm_diimfp_penn(delta_CS,osc,FermiEnergy,E0); 
NDIIMFP = F_mod(delta_CS,J_ion,Mlambda_ion,osc,E0,FermiEnergy);

%% Const for Ind_surf                  
 
osc.A = 5;   
osc.G = 0.8; 
osc.Om = 10;

%NDSEP = nrm_diimfp_penn(delta_CS,osc,FermiEnergy,E0);
NDSEP = F_mod(delta_CS,J_ion,Mlambda_ion,osc,E0,FermiEnergy);

%% DIIMFP and DSEP calculation
I_in=zeros(1, length(delta_CS));

I_in(1,:)=NDIIMFP;
I_in(2,:)=NDSEP;

Mat{1} = Material(Element,E0);
Mat{1}.SetManualDIIMFP(E0-delta_CS',NDIIMFP');

Mat{2} = Material(Element,E0);
Mat{2}.SetManualDIIMFP(E0-delta_CS',NDSEP');

x_in=struct('delta', delta_CS , 'I', I_in);
file_name = ['Data_x_in\x_in_Al_', num2str(round(E0/10)*10) ];
file_name=[file_name, '.mat'];
save(file_name,'x_in');

%% plot for DIIMFP
figure(1);
semilogy(x_in.delta,x_in.I,'-','LineWidth',2);
xlabel('Energy loss, eV'); ylabel('NDIIMFP, eV^{-1}');
legend('Bulk','Surface');

%% Specrtum calculation
Layers = [Layer(Mat{2},Thickness_s) Layer(Mat{1})];
%Layers = Layer(Mat{1});
Rez=Calculation_REELS(Layers, Geometry, Delta_max, sigma_gauss);

%% Plot for spectra
figure(2);
load Data_Al__175_40000eV_2geom; %experimental data loading
plot(E0-XY_all(:,9), XY_all(:,10)*MKgraf,'ro') %experimental data
hold on
plot(E0-Rez.energy_mesh_full-0.15,Rez.EnergyDistribution,'-','LineWidth',2) %calculation result
title(['E_0=', num2str(E0), ' eV'])
xlabel('{\itE}, eV')
ylabel('Intensity  {\itI}, rel.un.')
legend('Experimental data','Calculation')
grid on
