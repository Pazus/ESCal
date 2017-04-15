classdef Layer < handle
    %Layer  class contaiting information its thickness
    
    properties (Dependent = true, SetAccess = private)
        MaterialName = '';
    end
    properties (SetObservable = true)
        Material@CoreMaterial scalar;
        thickness@double  scalar = Inf;
    end
    
    properties (SetAccess = protected)
        tau_in@double;
        tau_tot@double;
    end
    
    properties (Access = private, Constant = true);
        defaultLayerThickness = Inf;
    end
    
    events
        LayerRecalculated;
    end
    
    methods
        function obj = Layer(Material, varargin)
            % LAYER Constructor
            % Layer(Material) constructs semi-infinite layer of the given material
            % Layer(Material, thickness)  constructs layers of THICKNESS nm of the given material 
            persistent ps;
            
            if isempty(ps)
                ps = inputParser;
                ps.FunctionName = 'Layer';
                ps.addRequired('Material', @(x) validateattributes(x, {'CoreMaterial'}, {'scalar','nonempty'}));
                ps.addOptional('thickness', obj.defaultLayerThickness, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'nonempty'}));
            end
            
            ps.parse(Material, varargin{:});

            obj.Material = ps.Results.Material;
            obj.thickness = ps.Results.thickness;
            
            addlistener(obj, 'thickness', 'PostSet', @obj.RecalculateTauThickness);
            addlistener(obj.Material, 'EnergyChanged', @obj.RecalculateTauThickness);
            
            obj.RecalculateTauThickness;
        end

        function RecalculateTauThickness(obj, varargin)
            if isempty(obj.Material.E0)
                obj.tau_in  = [];
                obj.tau_tot = [];
            else
                obj.tau_in  = obj.thickness ./ obj.Material.l_in;
                obj.tau_tot = obj.thickness ./ obj.Material.l_tot;
            end
            notify(obj, 'LayerRecalculated');
        end
        
        function val = isinf(obj)
            val = isinf(obj.thickness);
        end
        
        function val = isfinite(obj)
            val = isfinite(obj.thickness);
        end
        
        function val = ge(l1, l2)
            val = l1.thickness >= l2.thickness;
        end
        
        function val = gt(l1, l2)
            val = l1.thickness > l2.thickness;
        end
        
        function val = le(l1, l2)
            val = l1.thickness <= l2.thickness;
        end
        
        function val = lt(l1, l2)
            val = l1.thickness < l2.thickness;
        end
        
        function val = eq(l1, l2)
            val = l1.thickness == l2.thickness;
        end
        
        function val = ne(l1, l2)
            val = l1.thickness ~= l2.thickness;
        end
        
        function val = get.MaterialName(obj)
            val = obj.Material.Mat;
        end
        
        function val = isequal(l1, l2)
            val = l1.thickness == l2.thickness && isequal(l1.Material, l2.Material);
        end
        
        function set.thickness(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', 'nonempty'});
            obj.thickness = val;
        end
    end
    
end

