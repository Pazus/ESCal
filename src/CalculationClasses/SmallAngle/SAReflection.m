classdef SAReflection < SABase
    %Reflection function in Small Angle Approximation
    
    properties
        Rm
    end
    
    properties (SetAccess = private)
        Rl
    end
    
    methods
        function obj = SAReflection(Layer)
            obj = obj@SABase(Layer);
            obj.CalculationResultPropertyName   = 'Rm';
        end
        
    end
    
    methods (Access = protected)
        function BaseCalculation(obj)
            
            %R_l = zeros(obj.N, obj.N, obj.N_in + 1);
            Rl_full = zeros(obj.N,obj.N,obj.N_in+1,obj.LegandePolinomialsCount+1);
            % obj.Rl = zeros(obj.N,obj.N,obj.N_in+1,obj.LegandePolinomialsCount+1);
           
            
            mm = bsxfun(@plus,1./obj.mu_mesh,1./obj.mu_mesh');
            lambda = obj.Layer.Material.lambda;
            tau = obj.Layer.tau_tot;
            x_l = obj.x_l;

            % Рассчитанные значения получаются в нормировке int( [-1:1]R=1,
            % а должны получаться int( [-1:1]R=norm
            B = obj.norm./mm;
%             Base = Base*obj.norm;
            
            if ~isinf(obj.Layer) 
                e0 = expint(tau*mm);
            end  
            
            if obj.N_in>0 && isfinite(obj.Layer)
                gi = zeros(obj.N, obj.N, obj.N_in);
                for n=1:obj.N_in
                    gi(:,:,n) = gammainc(tau*mm,n);
                end
            end
            
            for i=1:obj.LegandePolinomialsCount+1
                
                R_l = zeros(obj.N, obj.N, obj.N_in + 1);
                R_l(:,:,1) = -log(1-lambda*x_l(i))*B;
                
                if ~isinf(obj.Layer)
                    R_l(:,:,1) = R_l(:,:,1) - ( expint( tau*(1-lambda*x_l(i))*mm )-e0 ) .* B;
                end
                
                if obj.N_in>0
                    if isinf(obj.Layer)
                        % n>0
                        for n=1:obj.N_in
                            R_l(:,:,n+1) = (1-lambda)^n/n*(1/(1-lambda*x_l(i))^n-1)*B;
                        end
                    else
                        for n=1:obj.N_in
                            % R_l(:,:,n+1) = (1-lambda)^n/n*(gammainc(tau*(1-lambda*x_l(i))*mm,n)./(1-lambda*x_l(i))^n-gammainc(tau*mm,n)).*B;
                            R_l(:,:,n+1) = (1-lambda)^n/n*(gammainc(tau*(1-lambda*x_l(i))*mm,n)./(1-lambda*x_l(i))^n-gi(:,:,n)).*B;

                        end
                    end
                end
                
                Rl_full(:,:,:,i) = R_l;
                
            end
            
            obj.Rl = Rl_full;
            
        end
        
        function ClearResults(obj,varargin)
            obj.ClearResults@SABase;
            obj.Rm         = [];
            obj.Rl         = [];
        end
        
        function ExpandCrossSections(obj)
            obj.x_l = obj.Layer.Material.Calculate_Leg_coefs(obj.LegandePolinomialsCount);
%             if obj.Layer.Material.IsManualDECS
%                 obj.x_l = ExpandFunctionByLegandrePolinomials(cos(obj.Layer.Material.DECS_mu'),obj.Layer.Material.DECS,obj.LegandePolinomialsCount);
%             else
%                 obj.x_l = obj.Layer.Material.Calculate_Leg_coefs(obj.LegandePolinomialsCount);
%                 obj.x_l = obj.x_l / obj.x_l(1);
%             end
        end
        
        
        function AllocateResults(obj)
            obj.Rm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
        end
        
        function AddExpansionCoef(obj,n,m,i,PlP)
            obj.Rm(:,:,n+1,m+1) = obj.Rm(:,:,n+1,m+1) + (-1)^(i-1+m)*obj.Rl(:, : ,n+1, i+m).*PlP;
        end
    end
    
end

