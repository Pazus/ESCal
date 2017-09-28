clear all;

%% Geometry settings
Geometry.theta0=0;
Geometry.theta=45;
Geometry.phi=0;

%% Parameters of the spectrum view
Element    =  'Al';
E0 = 2000;
Delta_max  = 50; %energy range value to calculate      
delta_CS   = 0:0.25:Delta_max; h=0.25; % h - step for DIIMFP                               

Thickness_s = 0.5; %thickness of the surface layer 
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
N_ion=length(J_ion);  % the number of ionization edges

%% Const for Ind_bulk 
% 1 means that these parameters are being used for the bulk
alpha_pl   = 26*0.1;
beta_pl    = 12*0.1; 
E_pl{1}    = 14.8; % the energy position of the bulk plasmon 
b{1}       = 0.18; % the width of the bulk plasmon peak  
alpha{1}   = alpha_pl; 
beta{1}    = beta_pl;  

%% Const for Ind_surf                  
% 2 means that these parameters are being used for the surface
E_pl{2}    = 10; % the energy position of the surface plasmon     
b{2}       = 0.8; % the width of the surface plasmon peak
alpha{2}   = alpha{1};
beta{2}    = beta{1};

%% Const for add Ind_q
E_pl_q    = [80 81   126  1567 10] - Work_func; % energy positions of all ionization edges
b_q       = [40 40   40   40   40]; % peak widths
alpha_q   = [0  0    0    0    0 ];           
beta_q    = [1  1    1    1    1 ];   
Mlambda_q = [0  0.02 0.04 0.03 0.07]; % the contribution of each edges

%% DIIMFP and DSEP calculation
I_in=zeros(1, length(delta_CS));
for i_layer = 1:length(E_pl)
            l_x_pl_int = true; %integral
    NDIIMFP = F1_mod(delta_CS, [E_pl{i_layer}; E_pl_q(:)],...
                               [b{i_layer}   ; b_q(:)   ],...
                               [1            ; Mlambda_q(:)],...
                               J_ion, Mlambda_ion(:),...
                               [alpha{i_layer}; alpha_q(:)],...
                               [beta{i_layer}; beta_q(:)],...
                     h, E0, l_x_pl_int);

    I_in(i_layer,:)=NDIIMFP;
    Mat{i_layer} = Material(Element,E0);
    Mat{i_layer}.SetManualDIIMFP(E0-delta_CS',NDIIMFP');
end
x_in=struct('delta', delta_CS , 'I', I_in);
file_name = ['Data_x_in\x_in_Al_', num2str(round(E0/10)*10) ];
file_name=[file_name, '.mat'];
save(file_name,'x_in');

%% plot for DIIMFP
figure(1);
plot(x_in.delta,x_in.I,'-','LineWidth',2);
xlabel('Energy loss, eV'); ylabel('NDIIMFP, eV^{-1}');
legend('Bulk','Surface');

%% Specrtum calculation
Layers = [Layer(Mat{2},Thickness_s) Layer(Mat{1})];
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
grid on
