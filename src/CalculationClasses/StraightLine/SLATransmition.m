classdef SLATransmition < SLAReflection
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Tm
    end
    
    methods
        function obj = SLATransmition(Layer)
            obj = obj@SLAReflection(Layer);
            obj.CalculationResultPropertyName   = 'Tm';
        end
    end
    
    methods (Access = protected)
        function ClearResults(obj,varargin)
            obj.ClearResults@SLABase;
            obj.Tm         = [];
        end
        
        function BaseCalculation(obj)
            mu = obj.mu_mesh;
            
            lambda = obj.Layer.Material.lambda;
            
            tau = obj.Layer.tau_tot;
            
            obj.Tm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            for m=0:0%obj.M
                n_factorial = 1;
                for n=0:obj.N_in
                    n_factorial = n_factorial*max(1,n);
                    obj.Tm(:, :, n+1, m+1) = diag(exp(-(1-lambda)*tau./mu).*((1-lambda)*tau./mu).^n/n_factorial);
                end
            end
            
            obj.BaseCalculation@SLAReflection;
        end
    end
    
    
end

