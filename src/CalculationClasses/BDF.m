classdef BDF < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        alpha
        betta
    end
    
    properties (SetObservable = true)
        order = 5;
    end
    
    methods
        function this = BDF(order)
            this.order = order;
            this.CalculateCoefs;
            addlistener(this,'order','PostSet',@this.CalculateCoefs);
        end
        
        function CalculateCoefs(this)
            betta_full = [1;  2/3; 6/11; 12/25; 60/137; 60/147];
            alpha_full = zeros(6);
            
            alpha_full(1,1)   = 1;
            alpha_full(2,1:2) = [4, -1]/3;
            alpha_full(3,1:3) = [18, -9, 2]/11;
            alpha_full(4,1:4) = [48, -36, 16, -3]/25;
            alpha_full(5,1:5) = [300, -300, 200, -75, 12]/137;
            alpha_full(6,1:6) = [360, -450, 400, -225, 72, -10]/147;
            
            this.alpha = alpha_full(1:this.order,1:this.order);
            this.betta = betta_full(1:this.order);
        end
    end
    
end

