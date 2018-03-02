classdef BaseMultiLayer < handle
    properties (SetObservable = true)
        
        Layers@Layer; % array of @Layer
        N_in=0;
        theta0                          = 0;
        N_Layer;% ?????????? ?????, ?????????? ???? ????? x_in ??? z
        %         mesh_E;
        CalculationMethod;
        ObjectsOfLayers;
        FullEnergyDistribution;
        EnergyDistribution;
        AngularDistribution
        
%         deltaE_elastic_peak;
        energy_mesh_full;
        mu_mesh;
        sigma_gauss;
        sigma_LDS;
        alpha_LDS;
        dE;
        N;
        M;
%        Calculate_only_theta0_theta=false;
    end
    
    properties (Hidden = true)
        Result_EnergyDistribution_temp ;
    end
    
    properties (SetAccess = protected)
        
        IsCalculated                    = false;                        % Object is calculated flag. If you change any parameter flag will be switched to false
    end
    
    methods
        function obj = BaseMultiLayer(Layers, CalculationMethod)
            obj.Layers = Layers;
            obj.N_Layer = numel(Layers);
            obj.CalculationMethod = CalculationMethod;

            if length(obj.CalculationMethod)==1 && obj.N_Layer>1
                for i_layer = 2:obj.N_Layer
                    obj.CalculationMethod{i_layer}=obj.CalculationMethod{1};
                end
            end
            if length(obj.CalculationMethod)~=obj.N_Layer
                error('Number of methods and number of layers must be the same.')
            end
            % ???????? ?? ??????? ?????????? ???????
            for i_layer=1:obj.N_Layer
               if isempty(obj.Layers(i_layer).Material.DIIMFP); error(['DIIMFP for ', num2str(i_layer), ' layer is not set.']); end 
            end
            %             addlistener(obj,'N_in','PostSet',@obj.Change_N_in);
            %             addlistener(obj,'theta0','PostSet',@obj.Change_theta0);
            %             addlistener(obj,'energy_mesh_full','PostSet',@obj.Setup_energy_mesh_full);
            
            
            %             addlistener(obj,'sigma_gauss','PostSet',@obj.ClearResults);
            %             addlistener(obj,'sigma_LDS','PostSet',@obj.ClearResults);
            %             addlistener(obj,'alpha_LDS','PostSet',@obj.ClearResults);
        end
        
        
        function plotEnergyDistribution(obj)
            plot(obj.energy_mesh_full,obj.EnergyDistribution);
            xlabel('{\itE}, eV');
        end
        
        function plotAngularDistribution(obj)

            if size(obj.AngularDistribution,2)==2        
                x_plot = [-acosd(fliplr(obj.mu_mesh)), acosd(obj.mu_mesh)];
                y_plot = [flipud(obj.AngularDistribution(:,2)); obj.AngularDistribution(:,1)];
            else
                x_plot = acosd(obj.mu_mesh);
                y_plot = obj.AngularDistribution;
            end
            plot(x_plot,y_plot);
            xlabel('{\theta}, eV');
        end
        
        %         function conv_angle(obj)
        %
        %         end
        
        function convGauss (obj,sigma_gauss)
            if  ~isempty(sigma_gauss)
                if sigma_gauss~=0
                    delta_temp=-5*sigma_gauss:obj.dE:5*sigma_gauss;
                    GaussF = normpdf(delta_temp, 0, sigma_gauss);
                    if trapz(delta_temp,GaussF)<0.995 || trapz(delta_temp,GaussF)>1.005 % ? ???????????? 0.5%
                        dE_new=round(abs(sigma_gauss/10),5,'significant');
                        delta_new = [fliplr(0:-dE_new:-5*sigma_gauss) dE_new:dE_new:5*sigma_gauss]; 
                        x_new=min(obj.energy_mesh_full):dE_new:max(obj.energy_mesh_full);
                        y_new=interp1(obj.energy_mesh_full,obj.EnergyDistribution,x_new);
                        clear GaussF
                        GaussF = normpdf(delta_new, 0, sigma_gauss);
                        F_temp = conv_my(y_new, GaussF,dE_new,'same');
                        F=interp1(x_new,F_temp,obj.energy_mesh_full);
                        warning(['Gauss convolution with dE=',num2str(dE_new),' eV!'])
                    else
                        F = conv_my(obj.EnergyDistribution, GaussF,obj.dE,'same');
                    end
                    obj.EnergyDistribution = F;
                end
            end
        end
        
        function convLDS (obj,sigma_LDS, alpha_LDS)
            if ~(isempty(alpha_LDS) || isempty(sigma_LDS)) 
               if ~(alpha_LDS==0 && sigma_LDS==0)
                    mesh_LDS = -obj.Layers(1).Material.E0:obj.dE:obj.Layers(1).Material.E0; % !!!!! ?
                    LDS = gamma(1-alpha_LDS)*cos(pi*alpha_LDS/2 + (1-alpha_LDS)*atan(mesh_LDS...
                        ./(sigma_LDS/2)))./(mesh_LDS.^2 + (sigma_LDS/2).^2).^((1-alpha_LDS)./2);
                    LDS = LDS/trapz(mesh_LDS,LDS);
                    F = conv_my(obj.EnergyDistribution, LDS, obj.dE,'same');
                    obj.EnergyDistribution = F;
               end
            end
        end
        
        
        
        
        
        
        function Solution_ABCD = DoubleLayer(obj,A,B,C,D,K)
            % A + B*C*D
            if ~isempty(D)
            if nargin < 6 && isempty(K); K = 0; end
            tic
            %  Solution_ABCD = A + conv3d(conv3d(B,C,1,'same')*obj.dE,D,1,'same')*obj.dE ;
            Solution_ABCD = A + conv3d(conv3d(B,C,'right',K)*obj.dE,D,'right',K)*obj.dE ;
            else
                if isempty(A) && ~isempty(B)
                    Solution_ABCD = conv3d(B,C,'right',K)*obj.dE;
                elseif isempty(A) && isempty(B)
                    Solution_ABCD = zeros(size(C));
                elseif ~isempty(A) && isempty(B)
                    Solution_ABCD = A;
                else
                    Solution_ABCD = A + conv3d(B,C,'right',K)*obj.dE ;
                end
            end
            time=toc;
            disp(['Calculating time of ABCD: ', num2str(time),' sec.']);
            
        end
        
        
    end
    
    methods %(Access = protected)
        function CalculateEnergyDistribution_E0(obj,theta,phi)
            if nargin < 3; phi = 0; end;
            if max(phi)>=360 || min(phi)<0
                error('Elements of phi_mush must be in [0;360) range');
            end
            phi = phi+180;
            if ~obj.IsCalculated; error('Not calculated yet.'); end;
            
            cos_theta0 = cosd(obj.theta0);
            cos_theta = cosd(theta);
            TempEnergyDistribution = zeros(1,size(obj.FullEnergyDistribution{1}.R,3));
            for i=1:size(obj.FullEnergyDistribution{1}.R,3);
                for m=0:size(obj.FullEnergyDistribution{1}.R,4)-1;
                    k = 2*cosd(m*phi);
                    if m == 0; k=k/2; end
%                     r0 = interp2(obj.mu_mesh,obj.mu_mesh,obj.Fm(:,:,i,m+1),cos_theta0,cos_theta); % old
                    r0 = interp2(obj.mu_mesh,obj.mu_mesh,obj.Fm(:,:,i,m+1),cos_theta,cos_theta0); % new?????
                    TempEnergyDistribution(i) = TempEnergyDistribution(i) + r0*k;
                end
            end
            y_temp = interp1(obj.ObjectsOfLayers{1}.energy_mesh,TempEnergyDistribution,obj.energy_mesh_full);
            y_temp(isnan(y_temp)) = 0;
            obj.EnergyDistribution = y_temp;
%             obj.EnergyDistribution =TempEnergyDistribution;
        end
        
        function CalculateAngularDistribution(obj,phi,theta0)
            if nargin < 2; phi = []; end; %??????
            if nargin < 3; theta0 = []; end;
            obj.CalculateAngularDistribution_E0(phi,theta0);
        end
        
        function CalculateAngularDistribution_E0(obj, phi, theta0)
            if nargin<2 || isempty(phi); phi = []; end
            if nargin<3 || isempty(theta0); theta0 = obj.theta0; end;
            
            if theta0 > obj.theta0; error('theta0 has to be not more then one that was used for calculations'); end;
            if ~obj.IsCalculated; error('Not calculated yet.'); end;
            
            cos_theta0 = cosd(theta0);
            [~, ind_E0] = min(abs(obj.energy_mesh_full-obj.Layers(1).Material.E0));  
            ind_mu = find(obj.mu_mesh == cos_theta0);
            if ~isempty(phi)
                phi = [0 180]+phi;
                if max(phi)>=360 || min(phi)<0
                    error('Elements of phi_mush must be in [0;360) range');
                end
                
                TempAngularDistribution = zeros(obj.ObjectsOfLayers{1}.N, numel(phi));
                
                % ?????????? ??? ????, ????? ???????? ????????? ?????????:
                %          \   /
                % theta0    \ /    theta
                %        ____O____
                
                phi = phi+180; %???
                
                for m=0:obj.ObjectsOfLayers{1}.M
                    k = 2*cosd(m*phi);
                    if m == 0; k=k/2; end
                    if ~isempty(ind_mu)
                        r0 = obj.Fm(ind_mu,:,ind_E0,m+1)';
                    else
                        r0 = interp1(obj.mu_mesh,obj.Fm(:,:,ind_E0,m+1),cos_theta0)';
                    end
                    TempAngularDistribution(:,:) = TempAngularDistribution(:,:) + reshape(r0*k,[obj.ObjectsOfLayers{1}.N, numel(phi)]);
                end
                
            else
                TempAngularDistribution = zeros(obj.ObjectsOfLayers{1}.N,1);
                if ~isempty(ind_mu)
                    TempAngularDistribution(:) = 2*pi*obj.Fm(ind_mu,:,ind_E0,1)';
                else
                    TempAngularDistribution(:) = 2*pi*interp1(obj.mu_mesh,obj.Fm(:,:,ind_E0,1),cos_theta0)';
                end
            end
            obj.AngularDistribution = TempAngularDistribution*obj.dE;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         function FullEnergyDistribution = RecalculateToEnergyDistribution(obj)
        function RecalculateToEnergyDistribution(obj)
            FED = cell(1,numel(obj.N_Layer));
            for i_layer = 1:obj.N_Layer
                tic;
                N = obj.ObjectsOfLayers{i_layer}.N;
                M = obj.ObjectsOfLayers{i_layer}.M;
                FED{i_layer}.L = zeros(N, N, size(obj.ObjectsOfLayers{i_layer}.energy_mesh,1),M+1);
                FED{i_layer}.R = zeros(N, N, size(obj.ObjectsOfLayers{i_layer}.energy_mesh,1),M+1);
                FED{i_layer}.T = zeros(N, N, size(obj.ObjectsOfLayers{i_layer}.energy_mesh,1),M+1);
                if isprop(obj.ObjectsOfLayers{i_layer}, 'Qm')
                    FED{i_layer}.Q = zeros(N, N, size(obj.ObjectsOfLayers{i_layer}.energy_mesh,1),M+1);
                    %                     FullEnergyDistribution.Qdown = zeros(N, N, size(EnDistr,1),M+1);
                end
                
                
                EnDistr = obj.ObjectsOfLayers{i_layer}.CalculateEnergyConvolutions';
                
%                 EnDistr(1,end) = 0; %!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                tau_tot = obj.ObjectsOfLayers{i_layer}.Layer.tau_tot;
                mu = obj.mu_mesh;
                lambda = obj.ObjectsOfLayers{i_layer}.Layer.Material.lambda;
                
                sL = size(FED{i_layer}.L);
                L = zeros(sL);

                calcR = isprop(obj.ObjectsOfLayers{i_layer}, 'Rm');
                calcT = isprop(obj.ObjectsOfLayers{i_layer}, 'Tm');
                calcQ = isprop(obj.ObjectsOfLayers{i_layer}, 'Qm');
                
                if calcR; R = zeros(sL); end;
                if calcT; T = zeros(sL); end;
                if calcQ; Q = zeros(sL); end;
                
                curLayer = obj.ObjectsOfLayers{i_layer};
                
                curLayerLm = zeros([sL(1:2) obj.N_in+1]);
                curLayerLm(:,:,1) = diag(exp(-tau_tot./mu));
                if isfinite(tau_tot) 
                    for j=1:obj.N_in
                        curLayerLm(:,:,j+1) = curLayerLm(:,:,j).*diag((1-lambda)*tau_tot/j./mu);
                    end    
                end    
                Lm = reshape(reshape(curLayerLm,sL(1)^2,[])*EnDistr,sL(1),sL(2),[]);
                
                for m=0:M

                    L(:,:,:,m+1) = Lm;
                    if calcR; R(:,:,:,m+1)   = reshape(reshape(curLayer.Rm(:,:,:, m+1),sL(1)^2,[])*EnDistr,sL(1),sL(2),[]); end;
                    if calcT; T(:,:,:,m+1)   = reshape(reshape(curLayer.Tm(:,:,:, m+1),sL(1)^2,[])*EnDistr,sL(1),sL(2),[]); end;
                    if calcQ; Q(:,:,:,m+1) = reshape(reshape(curLayer.Qm(:,:,:, m+1),sL(1)^2,[])*EnDistr,sL(1),sL(2),[]); end;
                    
                end
                
                FED{i_layer}.L = L;
                if calcR; FED{i_layer}.R = R; end;
                if calcT; FED{i_layer}.T = T; end;
                if calcQ; FED{i_layer}.Q = Q; end;
                %                 FullEnergyDistribution.Qdown = FullEnergyDistribution.Qup;
                time=toc;
                disp(['Calculating time (FED) of the layer ', num2str(i_layer), ': ', num2str(time),' sec.']);
            end
            
            obj.FullEnergyDistribution = FED;
        end
        
%         function Change_N_in(obj, varargin)
%             for i_layer = 1:obj.N_Layer
%                 obj.ObjectsOfLayers{i_layer}.N_in = obj.N_in;
%             end
%         end
        
%         function Change_theta0(obj, varargin)
%             for i_layer = 1:obj.N_Layer
%                 obj.ObjectsOfLayers{i_layer}.theta0 = obj.theta0;
%             end
%         end
        
        function Setup_energy_mesh_full(obj, varargin)
            p = 0;
            if size(obj.energy_mesh_full)~=size(obj.Layers(1).Material.DIIMFP_E)
                p = 1;
            end
            if ~p
                p = (max(abs(obj.energy_mesh_full-obj.Layers(1).Material.DIIMFP_E))~=0);
            end
            if p 
                for i_layer = 1:obj.N_Layer
                    
                    y_in = obj.Layers(i_layer).Material.DIIMFP;
                    x_in = obj.Layers(i_layer).Material.DIIMFP_E;
                    y_new=interp1(x_in, y_in, obj.energy_mesh_full);
                    y_new(isnan(y_new))=0;

                    obj.Layers(i_layer).Material.SetManualDIIMFP(obj.energy_mesh_full',y_new');
                end
            end
%             obj.dE = abs(obj.Layers(1).Material.DIIMFP_E(2) - obj.Layers(1).Material.DIIMFP_E(1));
            obj.dE = obj.energy_mesh_full(2)-obj.energy_mesh_full(1);

        end
        
        function CalculateEachLayer(obj)
            for i_layer = 1:obj.N_Layer
                obj.ObjectsOfLayers{i_layer}.theta0 = obj.theta0;
                obj.ObjectsOfLayers{i_layer}.N_in = obj.N_in;
                % {
                if ~isempty(obj.N) % && strcmp(obj.CalculationMethod(i_layer), 'NS')
                    obj.ObjectsOfLayers{i_layer}.N = obj.N;
                end
                %}                
                tic;
                obj.ObjectsOfLayers{i_layer}.Calculate;

                time=toc;
                disp(['Calculating time of the layer ', num2str(i_layer), ': ', num2str(time),' sec.']);                    
            end
            
        end
        
    end;
    
    methods (Abstract = true)
        Calculate(obj);
    end
end
