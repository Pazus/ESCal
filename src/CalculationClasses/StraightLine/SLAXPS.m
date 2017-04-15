classdef SLAXPS < SLATransmition
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Qm
        flipLayer = false;
    end
    
    properties (SetAccess = protected)
        F_m_m;
    end
    
    methods
        function obj = SLAXPS(Layer)
            if ~isa(Layer.Material, 'XPSMaterial'); error('Material of the Layer must be of XPSMaterial class'); end
            obj = obj@SLATransmition(Layer);
            obj.CalculationResultPropertyName   = 'Qm';
        end
    end
    
    methods (Access = protected)
        function CalculateCrossSection(obj, varargin)
            obj.CalculateCrossSection@SLATransmition;
            x = -1:0.001:1;
            
            F = obj.Layer.Material.CalculateDPCS(x);
            
            if obj.flipLayer
                F = flip(F);
            end

            [~, obj.F_m_m] = X_mp(F,x,obj.mu_mesh,obj.mu_mesh_weights,obj.M,1); %obj.F_p_m
        end
        function ClearResults(obj,varargin)
            obj.ClearResults@SLATransmition;
            obj.Qm         = [];
        end
        
        function BaseCalculation(obj)
            mu = obj.mu_mesh;
            
            lambda = obj.Layer.Material.lambda;
            lambda_photon = obj.Layer.Material.lambda_photon;
            
            tau = obj.Layer.tau_tot;
            pow = (1-lambda)./mu;
            
            K = zeros(obj.N_in+1, obj.N);
            K(1,:) = 1;
            if isfinite(tau)
                for n=1:obj.N_in
                    K(n+1,:) = K(n,:).*pow*tau/n;
                end
            end
            S = cumsum(K,1);
            
            obj.Qm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            for m=0:obj.M
                for n=0:obj.N_in
                    obj.Qm(:, :, n+1, m+1) = lambda*lambda_photon*bsxfun(@times, obj.F_m_m(:,:,m+1)/(1-lambda), (1-exp(-tau*pow).*S(n+1,:)).*mu);
                end
            end
            
            obj.BaseCalculation@SLATransmition;
        end
    end
    
    
end

