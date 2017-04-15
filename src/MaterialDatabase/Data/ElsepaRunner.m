classdef ElsepaRunner
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function [Res] = RunElsepa(Z, E0)
            old_path = cd;
            %             addpath([old_path filesep 'Elsepa']);
            
            
            %             dir_Elsepa = fullfile(old_path,'Elsepa');
            
            current_full_path = dbstack('-completenames');
            current_file_name = dbstack;
            ind = strfind(current_full_path(1).file,current_file_name(1).file);
            dir_Elsepa = [current_full_path(1).file(1:ind-2) filesep 'Elsepa'];
            f_name=fullfile(dir_Elsepa,'lub.in');
            
            if isempty(E0); return; end;
            
            fid = fopen(f_name,'w');
            fprintf(fid, 'IZ      %1.0f         atomic number                               [none]\n',Z);
            fprintf(fid, 'IELEC  -1          -1=electron, +1=positron                     [ -1]\n');
            
            for i=1:numel(E0);
                fprintf(fid, 'EV      %1.0f       kinetic energy (eV)                         [none]\n',E0(i));
            end
            
            fclose(fid);
            
            try
                
                cd(dir_Elsepa);
                !gelscata.exe > out.txt
                cd(old_path);
                
                DELIMITER = ' ';
                HEADERLINES = 31;
                
                %             Res = struct;
                for i=1:numel(E0);
                    
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    
                    
                    El = importdata(f_name_el, DELIMITER, HEADERLINES);
                    data = struct;
                    data.sigma_el = str2num(substr(El.textdata{22,1},33,11))*1e14;
                    data.sigma_tr1 = str2num(substr(El.textdata{23,1},33,11))*1e14;
                    data.sigma_tr2 = str2num(substr(El.textdata{24,1},33,11))*1e14;
                    
                    data.x = El.data(:,1)/180*pi;
                    data.y(:,i) = El.data(:,3)*10^18; %nanometers
                    Res(i) = data;
                    
                end
                
                for i=1:numel(E0);
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    delete(f_name_el);
                end
                
                
            catch ME
                for i=1:numel(E0);
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    delete(f_name_el);
                end
            end
        end
    end
    
end

