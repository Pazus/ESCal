classdef XPSMaterial < Material
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetObservable = true)
        BindingEnergy                   % Binding energy of the shell
        Anode@XPSAnode;
    end
    
    properties
        PCS                             % Photoelectron cross-section
        Total_PCS                       % Total photoelectron cross-section
        betta                           % Assymetric parametr
        gamma                           % two nondipole parameters
        delta
        %         DPCS                            % Differential photoemission cross-section
        %         DPCS_mu                         % Differential photoemission sross-section angle cosine
    end
    
    properties (SetAccess = protected)
        ShellIndex@uint8 scalar;        % Shell Index
        ShellName        = '';          % Shell name
        
        l_photon                        % Photoemission mean free path
        l_photon_tot                    % Photoemission mean free path
        lambda_photon                   % l_photon/l_el
        % TODO: multiple energy realization
        photon_E0                       % Energy of photon
        ShellsCount                     % Number of shells
    end
    
    properties (Access = private, Constant = true)
        defaultAnodeMaterialName = 'Al';
    end
    
    methods
        function obj = XPSMaterial(MatName,Shell,varargin)
            
            persistent ps;
            
            obj = obj@Material(MatName);
            
%             if isempty(ps)
                ps = inputParser;
                ps.FunctionName = 'XPSMaterial';
                
                ps.addRequired('Shell');
                
                ps.addOptional('Anode',XPSAnode(obj.defaultAnodeMaterialName),@validateAnodeParameter);
%             end;
            
            ps.parse(Shell,varargin{:});
            
            addlistener(obj, 'BindingEnergy', 'PostSet', @obj.SetE0);
            obj.SetAnode(ps.Results.Anode);
            
            addlistener(obj, 'Anode', 'PostSet', @obj.RecalculateXPS);
            obj.SetShell(ps.Results.Shell);
        end
        
        function SetShell(obj,Shell)
%             persistent ps;
            
            function validateShell(x)
                if iscell(x) && isscalar(x) && ischar(x{1}) || ischar(x)
                    validatestring(x,obj.Data.XPS.Shells,'SetShell','Shell');
                elseif isnumeric(x)
                    assert(x<=numel(obj.Data.XPS.Shells),'Shell index is out of boundary');
                else
                    error('Shell parameter has to be of char on numeric type');
                end
            end
            
            ps = inputParser;
            ps.FunctionName = 'XPSMaterial';

            ps.addRequired('Shell',@validateShell);

            ps.parse(Shell);
            
            if ischar(Shell)
                obj.ShellIndex = obj.GetShellIndexByName(Shell);
            elseif isnumeric(Shell)
                obj.ShellIndex = uint8(Shell);
            end
            
            obj.ShellName = obj.Data.XPS.Shells(obj.ShellIndex);
            
            obj.RecalculateXPS;
            
            obj.BindingEnergy = obj.Data.XPS.EB(obj.ShellIndex);
            
        end
        
        function SetAnode(obj,Anode)
            if isa(Anode, 'char');
                obj.Anode = XPSAnode(Anode);
            elseif isa(Anode,'XPSAnode');
                obj.Anode = Anode;
            end
        end
        
        function RecalculateXPS(obj,varargin)
            XPSData = obj.Data.XPS;
            
            PhotonEnergy = obj.Anode.PhotonEnergy;
            photoelectronEnergy = PhotonEnergy - XPSData.EB(obj.ShellIndex);
            
            obj.ShellsCount = sum(XPSData.EB<=PhotonEnergy);
            
            obj.PCS = interp1(XPSData.E, XPSData.PCS(:,obj.ShellIndex), photoelectronEnergy,'pchip');
            obj.betta = interp1(XPSData.E, XPSData.betta(:,obj.ShellIndex), photoelectronEnergy,'pchip');
            obj.gamma = interp1(XPSData.E, XPSData.gamma(:,obj.ShellIndex), photoelectronEnergy,'pchip');
            obj.delta = interp1(XPSData.E, XPSData.delta(:,obj.ShellIndex), photoelectronEnergy,'pchip');
            
            tocatCrossSection = 0;
            for i =1:obj.ShellIndex
                tocatCrossSection = tocatCrossSection + interp1(XPSData.E, XPSData.PCS(:,XPSData.EB<=PhotonEnergy), PhotonEnergy-XPSData.EB(i),'pchip');
            end
            obj.Total_PCS = tocatCrossSection;
            
            obj.l_photon = 1/(obj.Density*obj.PCS);
            
            obj.SetE0;
            
        end
        function CalculateAllMFP(obj)
            obj.CalculateAllMFP@Material();
            obj.l_photon = 1./(obj.Density * obj.PCS);
            obj.lambda_photon = obj.l_el ./ obj.l_photon;
            obj.l_photon_tot = 1./(obj.Density * obj.Total_PCS);
        end
        function F = CalculateDPCS(obj,mu) %,phi,M
            persistent ps;
            if isempty(ps)
                ps = inputParser;
                ps.FunctionName = 'CalculateDPCS';
                ps.addRequired('mu',@(x) validateattributes(x,{'numeric'},{'nonempty','>=',-1,'<=',180}));
            end
            ps.parse(mu);
            %             if nargin < 3
            %                 M = 0;
            %             end
            %             F = zeros(numel(mu),numel(mu),M+1);
            P = P_Leg(mu,3);
            if any(strcmp(obj.Anode.polarization,{'unpolarized','circular'}))
                F = (1 - obj.betta/2*P(:,3)'+mu.*(obj.gamma/2*(1-mu.^2)+obj.delta))/2;
            elseif strcmp(obj.Anode.polarization,'linear')
                error('Linear polarization not yet implemented');
                F = (1 - obj.betta*P(:,3)'*ones(1, numel(phi))+sqrt(1-mu.^2).*(obj.gamma*mu.^2+obj.delta)*cos(phi))/2;
            elseif strcmp(obj.Anode.polarization,'none')
                F = (1 - obj.betta/2*P(:,3)')/2;
            else
                error('Unknown polarization!');
            end
            
        end
        
        function plotDPCS(obj, mu)
            if nargin == 1
                mu = cosd(0:0.5:180);
            end
            
            mu = mu(:)';
            
            DPCS = obj.CalculateDPCS(mu);
            
            plot(mu,DPCS);
            xlabel('\theta');
            ylabel('x_{\gamma}');
        end
        
    end
    
    methods (Access = protected)
        function SetE0(obj,varargin)
            obj.E0 = obj.Anode.PhotonEnergy - obj.BindingEnergy;
        end
    end
    
    methods (Static = true)
        function ShellsCount = GetShellCount(MatName, PhotonEnergy)
            if nargin<2
                PhotonEnergy = [];
            end
            ShellsCount = numel(XPSMaterial.GetShellList(MatName, PhotonEnergy));
        end
        
        function Shells = GetShellList(MatName,PhotonEnergy)

            Data = load(Material.databaseFileName);
            validatestring(MatName, fieldnames(Data.MaterialData));
            
            XPS = Data.MaterialData.(MatName).XPS;
            
            if nargin<2 || isempty(PhotonEnergy)
                Shells = XPS.Shells;
            else
                Shells = XPS.Shells(XPS.EB<=PhotonEnergy);
            end
        end
        
        function Shells = GetShellIndexList(MatName,PhotonEnergy)

            Data = load(Material.databaseFileName);
            validatestring(MatName, fieldnames(Data.MaterialData));
            
            XPS = Data.MaterialData.(MatName).XPS;
            
            ind = (1:numel(XPS.Shells))';
            
            if nargin<2 || isempty(PhotonEnergy)
                Shells = ind;
            else
                Shells = ind(XPS.EB<=PhotonEnergy);
            end
        end
    end
    
    methods (Access = private)
        function ShellIndex = GetShellIndexByName(obj, ShellName)
            Ind = strcmpi(ShellName,obj.Data.XPS.Shells);
            if ~any(Ind)
                error(['Subshell with name "' ShellName '" not found']);
            end
            
            Enum = 1:numel(Ind);
            ShellIndex = uint8(Enum(Ind));
        end
    end
    
end

function validateAnodeParameter(x)
    assert(isa(x,'XPSAnode') || (iscell(x) && isscalar(x) && ischar(x{1})) || ischar(x), ...
        'Anode have to be of char or XPSAnode');
    
end


