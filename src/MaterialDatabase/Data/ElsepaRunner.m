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
            
            if isempty(E0); return; end
            
            fid = fopen(f_name,'w');
            fprintf(fid, 'IZ      %1.0f         atomic number                               [none]\n',Z);
            fprintf(fid, 'IELEC  -1          -1=electron, +1=positron                     [ -1]\n');
            
            for i=1:numel(E0)
                fprintf(fid, 'EV      %1.0f       kinetic energy (eV)                         [none]\n',E0(i));
            end
            
            fclose(fid);
            
            try
                
                cd(dir_Elsepa);
                system(['source/elscata.exe < ' f_name]);
%                 !gelscata.exe > out.txt
%                 !source/elscata.exe > out.txt
                cd(old_path);
                
                DELIMITER = ' ';
                HEADERLINES = 31;
                
                %             Res = struct;
                for i=1:numel(E0)
                    
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    
                    data = struct;
                    if E0(i)<100
                        El = importdata(f_name_el, DELIMITER, 34);
                        
                        str=El.textdata{25};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_el = str2double(cell2mat(newStr))*1e14;
                        
                        str=El.textdata{26};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_tr1 = str2double(cell2mat(newStr))*1e14;
                        
                        str=El.textdata{27};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_tr2 = str2double(cell2mat(newStr))*1e14;
                    else
                        El = importdata(f_name_el, DELIMITER, HEADERLINES);
                        
                        str=El.textdata{22};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_el = str2double(cell2mat(newStr))*1e14;
                        
                        str=El.textdata{23};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_tr1 = str2double(cell2mat(newStr))*1e14;
                        
                        str=El.textdata{24};
                        newStr = extractBetween(str,"= "," cm");
                        data.sigma_tr2 = str2double(cell2mat(newStr))*1e14;
                    end
                    data.x = El.data(:,1)/180*pi;
                    data.y(:,i) = El.data(:,3)*10^18; %nanometers
                    Res(i) = data;
                    
                end
                
                for i=1:numel(E0)
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    delete(f_name_el);
                end
                
                
            catch ME
                for i=1:numel(E0)
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                    delete(f_name_el);
                end
            end
        end
    end
    
end

