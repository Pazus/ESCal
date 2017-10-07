function Rez = Calculation_REELS (Layers, Geometry, Delta_max, sigma_gauss)
%% Solution by ReflectionMultiLayer

% Methods = {'NS'};
Methods = {'NS','NS'};
Rml = ReflectionMultiLayer(Layers,Methods);
Rml.theta0 = Geometry.theta0;
Rml.N_in = 10; %40
% Rml.N = 160; %Al:260 %Nb:25->160; 40->200
h=round(abs(Rml.Layers(1).Material.DIIMFP_E(2)-Rml.Layers(1).Material.DIIMFP_E(1)),3,'significant');
Rml.energy_mesh_full = [-Delta_max:h:5]+Rml.Layers(1).Material.E0;
Rml.sigma_gauss = sigma_gauss;
Rml.Calculate;
Rml.CalculateEnergyDistribution(Geometry.theta,Geometry.phi);
% Rez.x = Rml.energy_mesh_full;
% Rez.y = Rml.EnergyDistribution;
Rez=Rml;
end