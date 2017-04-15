classdef SLAReflection < SLABase
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Rm
    end
    
    methods
        function obj = SLAReflection(Layer)
            obj = obj@SLABase(Layer);
            obj.CalculationResultPropertyName   = 'Rm';
        end
    end
    
    methods (Access = protected)
        function ClearResults(obj,varargin)
            obj.ClearResults@SLABase;
            obj.Rm         = [];
        end
        
        function BaseCalculation(obj)
            mm = bsxfun(@plus,1./obj.mu_mesh,1./obj.mu_mesh');
            
            lambda = obj.Layer.Material.lambda;
            tau = obj.Layer.tau_tot;
            pow = mm*(1-lambda);
            
            K = zeros(obj.N, obj.N, obj.N_in+1);
            K(:,:,1) = 1;
            if isfinite(tau)
                for n=1:obj.N_in
                    K(:,:,n+1) = K(:,:,n).*pow*tau/n;
                end
            end
            S = cumsum(K,3);
            
            obj.Rm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            for m=0:obj.M
                for n=0:obj.N_in
                    obj.Rm(:, :, n+1, m+1) = lambda/(1-lambda)*obj.x_m_m(:,:,m+1).*(1-exp(-tau*pow).*S(:,:,n+1))./mm;
                end
            end
        end
    end
    
    
end

