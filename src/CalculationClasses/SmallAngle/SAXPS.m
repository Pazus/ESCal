classdef SAXPS < SATransmition
    %SAXPS2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Qm
        flipLayer = false;
    end
    
    properties (SetAccess = protected)
        F_l
    end
    
    properties (SetAccess = private)
        Ql
    end
    
    properties (Access = private, Constant)
        FlRound = 6;
    end
    
    methods
        function obj = SAXPS(Layer)
            if ~isa(Layer.Material,'XPSMaterial'); error('Material must be of XPSMaterial class'); end;
            obj = obj@SATransmition(Layer);
            obj.LegandePolinomialsCount = 3;
            obj.CalculationResultPropertyName   = 'Qm';
        end        
    end
    
    methods (Access = protected)
        function BaseCalculation(obj)                        
            mu = obj.mu_mesh;
            
            lambda = obj.Layer.Material.lambda;
            lambda_photon = obj.Layer.Material.lambda_photon;
            x_l = obj.x_l;
            tau = obj.Layer.tau_tot;

            K2 = (1-lambda)./(1-lambda*x_l);
            Base = lambda*lambda_photon*obj.F_l./(1-lambda*x_l);
            Base = Base * obj.norm;
            Q_l = zeros(obj.N, obj.N_in + 1, obj.LegandePolinomialsCount+1);
            
            for i=1:obj.LegandePolinomialsCount+1

                if isfinite(tau)
                    pow = (1-lambda*x_l(i))./mu';
                    ex = exp(-tau*pow);

                    K = zeros(obj.N,obj.N_in+1);
                    K(:,1) = 1;
                    for n=1:obj.N_in
                        K(:,n+1) = K(:,n).*pow*tau/n;
                    end
                    S = cumsum(K,2);
                else
                    S = 0;
                    ex = 0;
                end                

                for n=1:obj.N_in+1
                    Q_l( :, n, i) = Base(i).*K2(i)^(n-1).*mu';
                    if isfinite(tau)
                        Q_l( :, n, i) = Q_l( :, n, i).*(1-ex.*S(:,n));
                    end
                end
            end
            
            obj.Ql = Q_l;
            
            obj.BaseCalculation@SATransmition;
        end
        
        function ClearResults(obj,varargin)
            obj.ClearResults@SATransmition;
            obj.Qm         = [];
            obj.Ql         = [];
        end
        
        function ExpandCrossSections(obj)
            DPCS = obj.Layer.Material.CalculateDPCS(obj.mu_mesh_temp);
            
            if obj.flipLayer
                DPCS = flip(DPCS);
            end 
            
            obj.F_l = ExpandFunctionByLegandrePolinomials(obj.mu_mesh_temp,DPCS,obj.LegandePolinomialsCount,obj.mu_weight_temp);
            obj.F_l = round(obj.F_l * 10^obj.FlRound) / 10^obj.FlRound;
            
            obj.ExpandCrossSections@SATransmition;
        end
        
        function AllocateResults(obj)
            obj.Qm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            obj.AllocateResults@SATransmition;
        end
        
        function AddExpansionCoef(obj,n,m,i,PlP)
            obj.Qm(:,:,n+1,m+1) = obj.Qm(:,:,n+1,m+1) + (-1)^(i-1+m)*PlP*diag(obj.Ql(:, n+1, i+m));
            obj.AddExpansionCoef@SATransmition(n,m,i,PlP);
        end
    end
    
end

