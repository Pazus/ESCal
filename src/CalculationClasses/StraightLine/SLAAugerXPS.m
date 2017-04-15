classdef SLAAugerXPS < SLATransmition
    %SLAAugerXPS Класс решения задач фото-Оже-спектроскопии в
    %   квазиоднократном приближении
    
    properties
        Am
    end
    
    methods
        function obj = SLAAugerXPS(Layer)
            obj = obj@SLATransmition(Layer);
            obj.CalculationResultPropertyName   = 'Am';
        end
    end
    
    methods (Access = protected)
        
        function ClearResults(obj,varargin)
            obj.ClearResults@SLATransmition;
            obj.Am         = [];
        end
        
        function BaseCalculation(obj)
            mu = obj.mu_mesh;
            
            lambda = obj.Layer.Material.lambda;
            % Нужно брать из материала
            lambda_auge = 1; % obj.Layer.Material.lambda_photon;
            
            tau = obj.Layer.tau_tot;
            pow = (1-lambda)./mu;
            
            % Расчет коэффициентов для кратностей для слоя
            K = zeros(obj.N_in+1, obj.N);
            K(1,:) = 1;
            if isfinite(tau)
                for n=1:obj.N_in
                    K(n+1,:) = K(n,:).*pow*tau/n;
                end
            end
            S = cumsum(K,1);
            
            obj.Am = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            
            F = lambda*lambda_auge/2/pi;
            
            % Основная расчетная формула
            A = F/(1-lambda)*(1-bsxfun(@times, exp(-tau*pow), S))*diag(mu);
            
            A = repmat(A,[1,1,obj.N]);
            obj.Am(:, :, :, 1) = permute(A,[3,2,1]);
            
            obj.BaseCalculation@SLATransmition;
        end
    end
    
end

