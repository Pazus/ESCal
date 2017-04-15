classdef XPSData
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ShellIndex
        Cathode_Mat = '';
        ShellName = '';
    end
    properties (SetAccess = protected)
        
        BindingEnergy
        E0
        PCS
        betta
        Total_PCS
    end
    
    properties (Constant = true, GetAccess = private)
        DatabaseFileName = 'MaterialData';
    end
    
    methods
        function obj = XPSData(MatName,Shell,Cathode_Mat)
            
            if iscell(MatName); MatName = char(MatName); end
            if nargin<3 || isempty(Cathode_Mat); Cathode_Mat = 'Al'; end
            
            Data = obj.LoadXPSInfoFromDatabase(MatName);
            
            if ischar(Shell)
                obj.ShellIndex = GetShellIndexByName(Data,Shell);
            elseif isfloat(Shell)
                if Shell > numel(Data.Shells); error('Shell index out of boundaries'); end
                obj.ShellIndex = Shell;
            end
            
            obj.ShellName = Data.Shells(obj.ShellIndex);
            obj.Cathode_Mat = Cathode_Mat;
            obj.BindingEnergy = Data.EB(obj.ShellIndex);
            switch obj.Cathode_Mat
                case 'Mg'
                    obj.PCS = Data.PCS_Mg(obj.ShellIndex);
                    obj.betta = Data.betta_Mg(obj.ShellIndex);
                    obj.Total_PCS = sum(Data.PCS_Mg);
                    obj.E0 = 1200 - obj.BindingEnergy;
                case 'Al'
                    obj.PCS = Data.PCS_Al(obj.ShellIndex);
                    obj.betta = Data.betta_Al(obj.ShellIndex);
                    obj.Total_PCS = sum(Data.PCS_Al);
                    obj.E0 = 1486.6-obj.BindingEnergy;
                otherwise
                    error(['Unknown cathode material: ' obj.Cathode_Mat]);
            end
        end
        function F = CalculateDPCS(this,mu,norm,M)
            if nargin < 3
                norm = 1;
            end
            if nargin < 4
                M = 0;
            end
            
            N = numel(mu);
            F = zeros(N,N,M+1);
            
            for m=0:M
                P = Legendre_mu(mu,m,3);
                F(:,:,m+1) = (P(:,1)*P(:,1)' - this.betta/2*P(:,3)*P(:,3)');
            end
            
        end
    end
    
    methods (Access = private)
        function ShellIndex = GetShellIndexByName(Data,ShellName)
            Ind = strcmpi(ShellName,Data.Shells);
            if isempty(Ind)
                error(['Shell with name ' ShellName ' not found']);
            end
            
            Enum = 1:numel(Ind);
            ShellIndex = Enum(Ind);
        end
        function Data = LoadXPSInfoFromDatabase(this, MatName)
            LoadData = load(this.DatabaseFileName);
            Data = LoadData.MaterialData.(MatName).XPS;
        end
    end
    
end

