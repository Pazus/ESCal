classdef SATransmition < SAReflection
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Tm
    end
    
    properties (SetAccess = private);
        Tl
    end
    
    methods
        function obj = SATransmition(Layer)
            obj = obj@SAReflection(Layer);
            obj.CalculationResultPropertyName   = 'Tm';
        end       
    end
    
    methods (Access = protected)
        function BaseCalculation(obj)
            
%             assert(~isa(obj,'SATransmition') || obj.N_in == 0, 'SmallAngle solution for Transition function for inelastic scattering is not yeat developed.');
            
            lambda = obj.Layer.Material.lambda;
            tau = obj.Layer.tau_tot;
            x_l = obj.x_l;
            
            Base = lambda*x_l;
            
            % Рассчитанные значения получаются в нормировке int( [-1:1]T=1,
            % а должны получаться int( [-1:1]T=norm
%             Base = Base*obj.norm;
            
%             norm_leg = ((0:obj.LegandePolinomialsCount)*2+1)/2;
            
            T_l = zeros(obj.N, obj.N, obj.N_in + 1, obj.LegandePolinomialsCount+1);
            obj.Tm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            
            mu = obj.mu_mesh;
            N = obj.N;
            
            if isfinite(tau)
                                
                for i=1:obj.LegandePolinomialsCount+1
                    
                    mm_bottom = bsxfun(@minus,1./mu',(1-Base(i))./mu);
                    if all(all(mm_bottom ~= 0))
    %                     mm_bottom = bsxfun(@minus,1./mu,(1-Base(i))./mu');

    %                     T_l(:,:,1,i)   = Base(i) .* bsxfun(@minus, exp(-tau*(1-Base(i))./mu), exp(-tau./mu'))   ./ mm_bottom;
                        T_l_temp   = obj.norm*Base(i) .* bsxfun(@minus, exp(-tau*(1-Base(i))./mu), exp(-tau./mu'))   ./ mm_bottom;
                        T_l(:,:,1,i)   = symmetrize(T_l_temp);
                    end

                end

            end
            
            obj.Tl = T_l;

            obj.BaseCalculation@SAReflection;
        end
        function ClearResults(obj,varargin)
            obj.ClearResults@SAReflection;
            obj.Tm         = [];
            obj.Tl         = [];
        end
        
        function AllocateResults(obj)
            obj.Tm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
            obj.AllocateResults@SAReflection;
        end
        
        function AddExpansionCoef(obj,n,m,i,PlP)
            obj.Tm(:,:,n+1,m+1) = obj.Tm(:,:,n+1,m+1) + obj.Tl(:, : ,n+1, i+m).*PlP;
            obj.AddExpansionCoef@SAReflection(n,m,i,PlP);
        end
    end
    
end

