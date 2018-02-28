classdef Material < CoreMaterial
    % MATERIAL Class loading, culculating and containing all material
    % data from P. Kaplya scross-section database
    
    properties (SetAccess = protected)
        Mat = '';

        NvTPP;
        NvSGS;
        Eg;
        Ef;
        Ep;
        
        isManualDECSLegCoefs = false;

    end
    
    properties (Access = protected, Transient = true)
        Data;
    end
    
    properties (Access = private, Constant = true)
        defaultShieldingCoefficient = 0;
    end
    
    methods
        
        function obj = Material(MatName, varargin)
            % MATERIAL constructor
            % obj = Material(MatName, E0) creates Material object for
            % material MatName and Energy E0. Energy could be scalar or vector
            % obj = Material(MatName, E0, options)
            % options could be useTransportApproximation - boolean
            % creates Material object. Elastic scattering cross section is
            % isotropic (transport approximation).
            

            ps = inputParser;
            ps.FunctionName = 'Material';
            MatNameValidation = @(x) (iscell(x) && isscalar(x) && ischar(x{1})) || ischar(x);
            ps.addRequired('MatName', MatNameValidation);
            ps.addOptional('E0', [], @(x) validateattributes(x,{'numeric'},{'vector','positive'}));
            ps.addOptional('useTransportApproximation',false,@(x) validateattributes(x, {'logical'},{'scalar'}));

            
            ps.parse(MatName,varargin{:});

            obj.Mat = char(ps.Results.MatName);
            obj.E0 = ps.Results.E0;
            obj.useTransportApproximation = ps.Results.useTransportApproximation;
            
            addlistener(obj,'E0','PostSet',@obj.RecalculateAll);
            addlistener(obj,'useTransportApproximation','PostSet',@obj.RecalculateAll);
            
            obj.InitValues;
            
            obj.RecalculateAll;
        end
        
        function InitValues(obj)           
            obj.M = obj.Data.M;
            obj.Z = obj.Data.Z;
            obj.Density = obj.Data.Density;
            obj.NvTPP = obj.Data.NvTPP;
            obj.NvSGS = obj.Data.NvSGS;
            obj.Eg = obj.Data.Eg;
            obj.Ef = obj.Data.Ef;
            obj.Ep = obj.Data.Ep;
        end
        
        function [DECS_temp] = CalculateDECS(obj, mu)
            if isempty(obj.E0)
                obj.DECS = [];
                obj.DECS_mu = [];
                return;
            end
            
            DECSData = obj.Data.DECS;
            
            if max(obj.E0) > max(DECSData.E0) || min(obj.E0)<min(DECSData.E0)
                error('E0 out of boundaries for DECS calculation');
            end
            
            E0_temp = obj.E0;
            if isrow(E0_temp); E0_temp = E0_temp'; end;
            
            if nargin < 2; mu = DECSData.x; end;
            if ~isrow(mu); mu=mu'; end;
            
            if max(mu)<=1 && min(mu)>=-1
                mu = acos(mu);
            elseif max(mu)>pi && min(mu)>=0
                mu = mu/180*pi;
            end
            
            nE = numel(E0_temp);
            nMu = numel(mu);
            
            if obj.useTransportApproximation
                DECS_temp = ones(nMu, nE);
            else
                DECS_temp = interp2(DECSData.E0, DECSData.x, DECSData.y, E0_temp, mu,'spline');
            end
            
            norm = 1./trapz(cos(mu), DECS_temp, 1);
            if mu(1)<mu(2) %in case of decreasing vector usually formed by cos(\theta)
                norm = norm * -1;
            end
                
            DECS_temp = DECS_temp*diag(norm);
            
            if nargout == 0
                obj.DECS = DECS_temp;
                obj.DECS_mu = mu;
                
                obj.IsManualDECS = false;
                notify(obj,'DECSChanged');
            end;
        end
        
        
        function [x_l] = Calculate_Leg_coefs(obj,NLeg)
            
            if isempty(obj.E0)
                obj.DECS = [];
                obj.DECS_mu = [];
                return;
            end
            
            if obj.IsManualDECS && ~obj.isManualDECSLegCoefs
                error('Legendre expansion coefficients calculation for manually set DECS is not yet implemented');
            end
            
            if obj.isManualDECSLegCoefs 
                if size(obj.DECS_l,2) < NLeg+1
                    error('Not enought DECS legandre coefs were manually set');
                end
                x_l = obj.DECS_l;
                return;
            end
            
            if nargin<2
                NLeg=1500;
            end
            
            DECSData = obj.Data.DECS;
            
            if max(obj.E0) > max(DECSData.E0) || min(obj.E0)<min(DECSData.E0)
                error('E0 out of boundaries for DECS calculation');
            end
            
            E0_temp = obj.E0;
            if isrow(E0_temp); E0_temp = E0_temp'; end;
            nE = numel(E0_temp);
            
            x_l = zeros(NLeg+1,nE);
            
            if obj.useTransportApproximation 
                
                x_l(1,:) = 1;
            else
                mu = DECSData.x;
                
                
%                 DECS_orig = interp1(DECSData.E0(:), DECSData.y', E0_temp(:), 'spline')';
%                 DECS_orig = DECS_orig / trapz(mu,diag(sin(mu))*DECS_orig);
                
                
                N = max(2000,floor(NLeg*5));
                
                % Create temp mesh with smallsize mash for small angles 
%                 [x1,w1]=legzo(N,0,0.1/180*pi);
%                 [x2,w2]=legzo(N,0.1/180*pi,pi);
                [theta_temp, w_temp]=legzo(N,0,pi);
                
                theta_temp = theta_temp(:);
                w_temp = w_temp(:);
                

                
                % interp on new mesh
                % decs_temp=interp1(DECSData.x,DECS_orig,x_temp,'spline');

                decs_temp = interp2(DECSData.E0, DECSData.x, DECSData.y, E0_temp, theta_temp,'spline');
                
                % Calculate Legendres polynomial on that mashes
                P = Legendre_mu(cos(theta_temp'),0,NLeg);
%                 P_orig = Legendre_mu(cos(DECSData.x(:)'),0,NLeg);

                
                % Search for the regular approximate solution using
                % Henyey-Greenstein phase function
                % to increase accuracy
%                 gValue = zeros(2,nE); % Henyey-Greenstein phase function coeffitient
%                 l_ = 0:NLeg;
%                 options = optimoptions(@lsqcurvefit,'FunctionTolerance',1e-9,'Display','off');
%                 xl0 = P*spdiags(sin(x_temp).*w_temp,0,2*N,2*N)*decs_temp;
%                 normValue = xl0(1,:);
%                 xl0 = xl0*diag(1./normValue);

%                 border = 20 /180*pi; % border limit for better optimization;
%                 P_orig_short = P_orig(:,mu<=border);
%                 mu_short = mu(mu<=border)';
                
%                 for i =1:nE
%                     DECS_short = DECS_orig(mu<=border,i)';
% %                    gValue(i) = lsqcurvefit(@(g,xdata)g.^xdata,0.95,l_,xl0(:,i)',0,1,options);
%                     totalCrossSection = trapz(-cos(mu), DECS_orig(:,i)) ;
%                     p = lsqcurvefit(...
%                         @(p,xdata)p(1)*(l_ + 0.5).*(p(2).^l_)*P_orig_short,...
%                         [totalCrossSection, 0.95]...
%                         ,mu_short, DECS_short...
%                         ,[1,0], [1, 1]...
%                         ,options);
%                     gValue(:,i) = p(:);
%                 end
                
%                 gl = bsxfun(@power,gValue(2,:),l_');
%                 DECS_reg = P_orig'*spdiags(l_'+0.5,0,NLeg+1,NLeg+1)*gl;
%                 
%                 DECS_irreg = diag(1./gValue(1,:))*DECS_orig-DECS_reg;
%                 
%                 decs_ir=interp1(DECSData.x,DECS_irreg,x_temp,'spline');
%                 
%                 xl_ir = P*spdiags(sin(x_temp).*w_temp,0,N,N)*decs_ir;
%                 
%                 x_l = gl + xl_ir;   
                
                
                x_l = P*spdiags(sin(theta_temp).*w_temp,0,N,N)*decs_temp;
                x_l = x_l*diag(1./x_l(1,:));
                
            end
            
            if nargout == 0
                obj.DECS_l = x_l;
            end
            
        end
        
        function setManualDECSLegCoefs(obj, coefs) 
            obj.DECS_l = coefs(:)';
            obj.isManualDECSLegCoefs = true;
            
            if ~obj.IsManualDECS
                if isempty(obj.DECS_mu)
                    obj.DECS_mu = 0:pi/360:pi;
                    obj.DECS_mu(end) = pi;
                end
                    
                P = Legendre_mu(cos(obj.DECS_mu),0,numel(coefs)-1);
                n_ = (0:(numel(coefs)-1))+0.5;
                DECS_temp = (coefs(:)'.*n_)*P;
                DECS_temp = DECS_temp / abs(trapz(cos(obj.DECS_mu),DECS_temp));
                obj.SetManualDECS(obj.DECS_mu(:), DECS_temp');
            end
        end        
        
        function CalculateDIIMFP(obj, varargin)
            % CalculateDIIMFP  Calculate Differential Inelastic Mean Free Path (DIIMFP) for the Material.
            %   CalculateDIIMFP() Calculates DIIMFP for a standard energy mesh 0:1:E0.
            %   CalculateDIIMFP(EnergyMesh) Calculates DIIMFP for a given energy mesh.
            %   CalculateDIIMFP(EnergyMesh, SheildingCoefficient) Calculates DIIMFP for a given energy mesh using given Sheilding coefficient (default 0).


            ps = inputParser;
            ps.FunctionName = 'CalculateDIIMFP';
            ps.addOptional('EnergyMesh', (0:ceil(obj.E0))', @(x) validateattributes(x,{'numeric'},{'nonnegative','nonempty'}));
            ps.addOptional('SheildingCoefficient', obj.defaultShieldingCoefficient, @(x) validateattributes(x,{'numeric'},{'scalar'}));

            ps.parse(varargin{:});
            
            if isempty(obj.E0)
                obj.DIIMFP = [];
                obj.DIIMFP_E = [];
                obj.DIIMFP_Shielding_koef = [];
                return;
            end
            
            EnergyMesh = reshape(ps.Results.EnergyMesh,[],1);
            SheildingCoefficient = ps.Results.SheildingCoefficient;
            
            E_old = EnergyMesh;
            [x_w,y_w] = obj.CalculateDIIMFP_Werner;
            
            
            % ѕродливаем сетку до нул€ с аналогичным шагом
            if min(E_old)>0
                dx = abs(EnergyMesh(end)-EnergyMesh(end-1));
                EnergyMesh = [flipud((min(EnergyMesh):-dx:0)'); EnergyMesh(2:end)];
            end
            
            prolongation_power = 2+SheildingCoefficient;
            
            % ¬ыполн€ем расчет дл€ каждого значени€ E0
            y = zeros(numel(EnergyMesh),numel(obj.E0));
            for i=1:numel(obj.E0)
                temp_E0 = obj.E0(i);
                delta_E = EnergyMesh(EnergyMesh<=temp_E0);
                ind = (delta_E>=max(x_w));
                y_l = zeros(size(delta_E));
                y_l(~ind) = interp1(x_w,y_w(:,i),delta_E(~ind),'pchip');
                y_l(ind) = y_w(end,i)*max(x_w)^prolongation_power./EnergyMesh(ind).^prolongation_power;
                
                y_l = y_l/trapz(delta_E,y_l);
                y(EnergyMesh<=temp_E0,i) = flipud(y_l);
            end
            
            obj.DIIMFP = interp1(EnergyMesh,y,E_old);
            obj.DIIMFP_E = E_old;
            obj.DIIMFP_Shielding_koef = SheildingCoefficient;
            
            obj.IsManualDIIMFP = false;
            notify(obj,'DIIMFPChanged');
        end
        
        function val = get.Data(obj)
            val = MaterialPropertiesDatabase.getData(obj.Mat);            
        end
        
        function Res =  isequal(mat1, mat2)
            Res =  strcmpi(mat1.MatName,mat2.MatName) && isequal(mat1.E0,mat2.E0);
        end
    end
    
    methods (Access = protected)
        function l_in = CalculateIMFP(obj)
            if isempty(obj.E0)
                obj.l_in = [];
                return;
            end
            
            recalculatedDensity = obj.Density*obj.M/6.022/10^23*10^21;
            
            Ep = 28.8*sqrt(obj.NvTPP*recalculatedDensity/obj.M);
            U = Ep^2/829.4;
            D = 53.4-20.8*U;
            C = 1.97 - 0.91*U;
            gamma = 0.191*recalculatedDensity^-0.5;
            betta = -0.10+0.944/sqrt(Ep^2+obj.Eg^2)+0.069*recalculatedDensity^0.1;
            
            l_in = obj.E0./( Ep^2*(betta*log(gamma*obj.E0)-C./obj.E0+(D./obj.E0.^2)) )*10^-1;
            obj.l_in = l_in;
            
        end
        
        function l_el = CalculateEMFP(obj)
            
            if isempty(obj.E0)
                obj.l_el = [];
                return;
            end
            ElasticData = obj.Data.Elastic;

            if obj.useTransportApproximation
                l_el = obj.CalculateTRMFP;
                obj.l_el = l_el;
            else
                if max(ElasticData.x) < max(obj.E0) || min(ElasticData.x) > min(obj.E0) || any(isnan(obj.E0))
                    error('E0 out of boundaries for l_el calculcation');
                end
                l_el = interp1(ElasticData.x,ElasticData.l_el,obj.E0);
                obj.l_el = l_el;
            end
        end
        
        function l_tr = CalculateTRMFP(obj)
            if isempty(obj.E0)
                obj.l_tr = [];
                return;
            end
            
            TRMFPData = obj.Data.Elastic;
            
            if max(TRMFPData.x) < max(obj.E0) || min(TRMFPData.x) > min(obj.E0) || any(isnan(obj.E0))
                error('E0 out of boundaries for l_tr calculcation');
            end
            
            l_tr = interp1(TRMFPData.x,TRMFPData.l_tr,obj.E0);
            obj.l_tr = l_tr;
        end
        
    end
    
    methods (Access = private)
        
        function [ x,y ] = CalculateDIIMFP_Werner(obj)
            DIIMFPData = obj.Data.DIIMFP;
            
            if max(obj.E0) > max(DIIMFPData.E0) || min(obj.E0)<min(DIIMFPData.E0)
                % error('E0 out of borders for DIIMFP calculation');
                disp('E0 out of borders for Werner DIIMFP calculation. Data extrapolation might lead to errors.');
            end
            E0_ = max(min(max(DIIMFPData.E0),obj.E0),min(DIIMFPData.E0));
            
            x = DIIMFPData.x;
            y = interp1(DIIMFPData.E0,log(DIIMFPData.y'),E0_)';
            y = exp(y);
            y(isnan(y)) = 0;
            y(x==0,:) = 0;
            
            x = x(1:end-1);
            y = y(1:end-1,:);
        end
        
    end
    
    methods (Static = true)
        function MaterialList = GetMaterialList()
            Data = load(Material.databaseFileName);
            
            MaterialList = fields(Data.MaterialData);
        end
    end
    
end