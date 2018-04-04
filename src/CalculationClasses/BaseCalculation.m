classdef BaseCalculation < handle
    
    properties (SetObservable = true, AbortSet = true)
        N_in                            = 0;            % Number of inelastic scatterings to calculate
        M                               = 0;            % Number of azimuthal harmonics
        N                               = 81;           % Number of -1:0 discritization points
        theta0                          = 0;            % Cosine of initial angle (needed to calculate appropriate M)
        norm                            = 1/2/pi;       % DECS norm (int [-1:1]x_el = norm)
        
        CalculationResultPropertyName   = '';                           % Name of main result property. This property will be plotted using protting functions
    end
    
    properties (SetAccess = protected)
        mu_mesh;                                                        % Angle cosine mesh
        mu_mesh_weights;                                                % Integration weights for angle cosine mesh
        
        AngularDistribution;                                            % Angular distribution for the given geometry (N_in,[phi])
        EnergyDistribution;                                             % Energy distribution for the given geometry (theta,[phi])
        EnergyDistributionByInelasticScatteringNumber;                  
        InelasticScatteringDistribution;                                % Distribution by Inelastic scattering number (not windowed by energy range, correct only for the elastically scattered electrons)
        InelasticScatteringDistributionForEnergyRange;                  % Distribution of inelastic scattering
        
        IsCalculated                    = false;                        % Object is calculated flag. If you change any parameter flag will be switched to false
    end
    
    methods
        function obj = BaseCalculation
            
            addlistener(obj,'N','PostSet',@obj.CalculateAngleMesh);
            addlistener(obj,'M','PostSet',@obj.ClearResults);
            addlistener(obj,'N_in','PostSet',@obj.ClearResults);
            addlistener(obj,'norm','PostSet',@obj.ClearResults);
            addlistener(obj,'theta0','PostSet',@obj.CalculateOptimalM);
            
            obj.CalculateMeshes;
            
        end;
        
        function plotAngularDistribution(obj,N_in,phi,theta0)
            % Plot Angular distribution plot for a chosen number of
            % inelastic scatterings, phi and theta0
            if nargin < 2 || isempty(N_in); N_in = 0; end;
            if nargin < 3; phi = []; end;
            if nargin < 4 || isempty(theta0); theta0 = obj.theta0; end;
                
            [x_plot, y_plot] = obj.prepareAngularDistributionPlot(N_in,phi,theta0);
                
            
            plot(x_plot,y_plot);
            if ~isempty(phi)
                xlim([-90 90])
            else
                xlim([0 90])
            end
            xlabel('\theta')
        end
        
        function plotAngularDistributionPolar(obj,N_in,phi,theta0)
            % Plot Angular distribution plot for a chosen number of
            % inelastic scatterings, phi and theta0
            if nargin < 2 || isempty(N_in); N_in = 0; end;
            if nargin < 3; phi = 0; end;
            if nargin < 4 || isempty(theta0); theta0 = obj.theta0; end;
            [x_plot, y_plot] = obj.prepareAngularDistributionPlot(N_in,phi,theta0);
            
            
            polar((x_plot'-min(x_plot))/180*pi,y_plot);
            
        end
        
        function plotInelasticScatteringDistribution(obj, theta, phi)
            % Plot distribution by number of inelastic scatterings for the
            % chosen theta and phi
            if nargin < 3; phi = []; end;
            CalculateInelasticScatteringDistribution(obj,theta,phi);
            plot(0:obj.N_in, obj.InelasticScatteringDistribution);
            xlabel('N_{in}');
            ylabel('Cn');
        end
        
        function plotInelasticScatteringDistributionForEnergyRange(obj,theta,phi)
            % Plot distribution by number of inelastic scatterings for the
            % chosen theta and phi. Each point is multiplyed by integral of
            % energy distribution. This takes into account that after
            % several inelastic scatterings electron can lose all its
            % initial energy
            if nargin < 3; phi = []; end;
            obj.CalculateInelasticScatteringDistributionForEnergyRange(theta,phi);
            
            plot(0:obj.N_in, obj.InelasticScatteringDistributionForEnergyRange);
            xlabel('N_{in}');
            ylabel('Cn_E');
        end
        
        function plotEnergyDistribution(obj, theta, phi, N_in)
            % Plot energy distribution for the chosen theta and phi.
            % With parameter N_in you can define which numbers of inelastic
            % scatterings to include. By default all calculated inelastic
            % scattering quantities are included
            if nargin < 3; phi = []; end;
            if nargin < 4 || isempty(N_in); N_in = 0:obj.N_in; end;
            if iscolumn(N_in); N_in = N_in'; end;
            
            obj.CalculateEnergyDistribution(theta,phi);
            
            plot(obj.energy_mesh,sum(obj.EnergyDistribution(:,N_in+1),2));
            xlabel('E, eV');
        end;
        
        function plotIntegratedEnergyDistribution(obj,intTheta)
            y_plot = obj.CalculateIntegratedEnergyDistribution(intTheta);
            x_plot = obj.energy_mesh;
            
            plot(x_plot,y_plot);
        end;
        
        function CalculateAngularDistribution(obj, N_in, phi, theta0)
            if nargin<2 || isempty(N_in); N_in = 0; end;
            if nargin<3; phi = []; end
            if nargin<4 || isempty(theta0); theta0 = obj.theta0; end;
            
            if theta0 > obj.theta0; error('theta0 has to be not more then one that was used for calculations'); end;
            if N_in > obj.N_in; error('Too big number of inelastic collision'); end;
            if ~obj.IsCalculated; error('Not calculated yet.'); end;
            
            calc_field = obj.(obj.CalculationResultPropertyName);
            cos_theta0 = cosd(theta0);
            
            if ~isempty(phi)
                
                if max(phi)>=360 || min(phi)<0
                    error('Elements of phi_mush must be in [0;360) range');
                end
                
                TempAngularDistribution = zeros(obj.N, numel(N_in), numel(phi));
                
                % ���������� ��� ����, ����� �������� ��������� ���������:
                %          \   /
                % theta0    \ /    theta
                %        ____O____
                
                phi = phi+180;
                
                for i=1:numel(N_in)
                    
                    for m=0:obj.M
                        k = 2*cosd(m*phi);
                        if m == 0; k=k/2; end
                        ind = find(obj.mu_mesh == cos_theta0);
                        
                        if ~isempty(ind)
                            r0 = calc_field(ind,:,N_in(i)+1,m+1)';
                        else
                            r0 = interp1(obj.mu_mesh,calc_field(:,:,N_in(i)+1,m+1),cos_theta0)';
                        end
                        
                        TempAngularDistribution(:,i,:) = TempAngularDistribution(:,i,:) + reshape(r0*k,[obj.N, 1, numel(phi)]);
                    end
                end
                
            else
                TempAngularDistribution = zeros(obj.N, numel(N_in));
                
                for i=1:numel(N_in)
                    ind = find(obj.mu_mesh == cos_theta0);
                    if ~isempty(ind)
                        TempAngularDistribution(:,i) = calc_field(ind,:,N_in(i)+1,1)*2*pi;
                    else
                        TempAngularDistribution(:,i) = interp1(obj.mu_mesh,calc_field(:,:,N_in(i)+1,1),cos_theta0)*2*pi;
                    end
                end
            end
            obj.AngularDistribution = TempAngularDistribution;
        end;
        
        function [x_plot, y_plot] = prepareAngularDistributionPlot(obj,N_in,phi,theta0)
            % Plot Angular distribution plot for a chosen number of
            % inelastic scatterings, phi and theta0
            if nargin < 2 || isempty(N_in); N_in = 0; end;
            if nargin < 3; phi = []; end;
            if nargin < 4 || isempty(theta0); theta0 = obj.theta0; end;
            assert(theta0<=obj.theta0 | isempty(phi),'theta0 parameter is more then initial theta0, results will be inaccurate.');
            
            if ~isempty(phi)
                if ~isscalar(phi); error('phi must be a scalar'); end;
                if max(phi)>=180 || min(phi)<0; error('Elements of phi_mush must be in [0;180) range'); end;
                
                phi_mesh = [0 180]+phi;
                obj.CalculateAngularDistribution(N_in, phi_mesh, theta0);
                
                x_plot = [-acosd(fliplr(obj.mu_mesh)), acosd(obj.mu_mesh)];
                y_plot = [flipud(obj.AngularDistribution(:,:,2)); obj.AngularDistribution(:,:,1)];
            else
                obj.CalculateAngularDistribution(N_in, [], theta0);
                
                x_plot = acosd(obj.mu_mesh);
                y_plot = obj.AngularDistribution;
            end
        end
        
        function CalculateInelasticScatteringDistribution(obj,theta,phi)
            if nargin < 3; phi = []; end;
            obj.CalculateAngularDistribution(0:obj.N_in,phi);
            for i=1:obj.N_in+1
                obj.InelasticScatteringDistribution(i) = interp1(obj.mu_mesh,obj.AngularDistribution(:,i,1),cosd(theta));
            end
        end
        
        function CalculateInelasticScatteringDistributionForEnergyRange(obj,theta,phi)
            if nargin < 3; phi = []; end;
            obj.CalculateInelasticScatteringDistribution(theta,phi);
            EnergyConvolutions = obj.CalculateEnergyConvolutions;
            
            LoweringCoeffitients = trapz(obj.Layer.Material.DIIMFP_E, EnergyConvolutions, 1);
            LoweringCoeffitients(1) = 1;
            
            obj.InelasticScatteringDistributionForEnergyRange = obj.InelasticScatteringDistribution .* LoweringCoeffitients;
        end
        
        function EnergyConvolutions = CalculateEnergyConvolutions(obj)
            if isempty(obj.Layer.Material.DIIMFP); obj.Layer.Material.CalculateDIIMFP; end;
            
            x = obj.energy_mesh;
            
            DIIMFP = obj.Layer.Material.DIIMFP;
            dE = x(2)-x(1);
            
            EnergyConvolutions = zeros(numel(x),obj.N_in+1);
            
            E0_ind = abs(x-obj.Material.E0) == min(abs(x-obj.Material.E0));
            EnergyConvolutions(E0_ind,1) = 1./dE;
            
            i_E0 = find(E0_ind);
            
            for i=2:obj.N_in+1
                EnergyConvolutions(:,i) = conv_my(EnergyConvolutions(:,i-1), DIIMFP, dE,'right',i_E0);
            end
        end
        
        function CalculateEnergyDistribution(obj,theta,phi)
            if nargin < 3; phi = []; end;
            obj.CalculateInelasticScatteringDistribution(theta,phi);
            
            Cn = obj.InelasticScatteringDistribution;
            TempEnergyDistribution = obj.CalculateEnergyConvolutions;
            
            obj.EnergyDistribution = TempEnergyDistribution * diag(Cn);
        end;
        
        function R_e = CalculateIntegratedEnergyDistribution(obj,intTheta)
            mu_ = cosd(0:intTheta/100:intTheta);
            R_k = interp1(obj.mu_mesh,squeeze(obj.Rm(1,:,:,1)),mu_);
            
            R_e = trapz(mu_,-R_k*obj.CalculateEnergyConvolutions'*2*pi);
        end
        
        function CalculateOptimalM(obj,varargin)
            obj.M = optimal_M(obj.theta0);
        end
        
        function FullEnergyDistribution = CalculateFullEnergyDistribution(obj)
            EnDistr = obj.CalculateEnergyConvolutions;
            
            tau_tot = obj.Layer.tau_tot;
            mu = obj.mu_mesh;
            lambda = obj.Layer.Material.lambda;

            FullEnergyDistribution.E = zeros(obj.N, obj.N, size(EnDistr,1),obj.M+1);
            FullEnergyDistribution.R = zeros(obj.N, obj.N, size(EnDistr,1),obj.M+1);
            FullEnergyDistribution.T = zeros(obj.N, obj.N, size(EnDistr,1),obj.M+1);
            FullEnergyDistribution.Qup = zeros(obj.N, obj.N, size(EnDistr,1),obj.M+1);
            FullEnergyDistribution.Qdown = zeros(obj.N, obj.N, size(EnDistr,1),obj.M+1);
            for m=0:obj.M
                for j=0:obj.N_in
                    FullEnergyDistribution.E(:,:,:,m+1) = FullEnergyDistribution.E(:,:,:,m+1) + bsxfun(@times,...
                        diag(exp(-tau_tot./mu).*((1-lambda)*tau_tot./mu).^(j)/factorial(j)),...
                        reshape(EnDistr(:,j+1),1,1,[]));
                    
                    if isprop(obj, 'Rm')
                        FullEnergyDistribution.R(:,:,:,m+1) = FullEnergyDistribution.R(:,:,:,m+1) + ...
                            bsxfun(@times, obj.Rm(:,:,j+1, m+1), reshape(EnDistr(:,j+1),1,1,[]));
                    end
                    if isprop(obj, 'Tm')
                        FullEnergyDistribution.T(:,:,:,m+1) = FullEnergyDistribution.T(:,:,:,m+1) + ...
                            bsxfun(@times, obj.Tm(:,:,j+1, m+1), reshape(EnDistr(:,j+1),1,1,[]));
                    end
                    if isprop(obj, 'Qm')
                        FullEnergyDistribution.Qup(:,:,:,m+1) = FullEnergyDistribution.Qup(:,:,:,m+1) + ...
                            bsxfun(@times, obj.Qm(:,:,j+1, m+1), reshape(EnDistr(:,j+1),1,1,[]));
                    end
                end
            end
            FullEnergyDistribution.Qdown = FullEnergyDistribution.Qup;
        end
        
    end
    
    methods (Access = protected)
        function ClearResults(obj, varargin)
            obj.IsCalculated    = false;
            
            obj.EnergyDistribution    = [];
            obj.AngularDistribution    = [];
            obj.InelasticScatteringDistribution = [];
        end;
        
        function CalculateAngleMesh(obj,varargin)
            [x,s] = legzo_n1(obj.N-1);
            
            x = [1,x];
            s = [0,s];
            
            obj.mu_mesh = (x+1)/2;
            obj.mu_mesh_weights = s/2/obj.norm;
            
            obj.ClearResults;
        end;
        
        function CalculateMeshes(obj)
            obj.CalculateAngleMesh();
        end;
        
    end
    
    methods (Abstract = true)
        Calculate(obj);
    end
    
end

