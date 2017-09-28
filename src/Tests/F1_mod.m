function x_in=F1_mod(delta,epl,b,lambda_pl,J_ion,lambda_ion,alpha,beta,h,E0, l_x_pl_int)
% суммирование lambda (Грязев)
if nargin==10; l_x_pl_int = false; end
if length(lambda_pl)>1; if sum(lambda_pl(2:end))>=1; error('Wrong lambda_q!'); end; end;
if sum(lambda_ion)>=1; error('Wrong lambda_ion!'); end; 

if max(delta)>10e3; delta_loc = 0:h:10e3;
else delta_loc = delta;
end

x_pl_all = zeros(length(alpha),length(delta_loc));

for i=1:length(alpha)
%     if i>1; l_x_pl_int = false; end
    if ~l_x_pl_int
        w_pl = delta_loc .^ beta(i) ./ ((delta_loc .^ 2 - epl(i) ^ 2) .^ 2 + b(i) ^ (4 - alpha(i)) * delta_loc .^ alpha(i));
    else
        tic
        w_pl = x_in_pl_int_function(delta_loc, E0, epl(i), b(i), alpha(i), beta(i), 20000); %5000 20
        toc
    end
    x_pl = w_pl ./ trapz(delta_loc,w_pl);
    x_pl_all( i, : ) = x_pl;
end
x_ion_all = zeros(length(J_ion),length(delta_loc));

for i=1:length(J_ion)
    if ~(J_ion(i)>max(delta_loc))
        w_ion = 1 ./ (delta_loc.^2);
        w_ion(1:sum(delta_loc<J_ion(i))) = 0;
        x_ion = w_ion ./ trapz(delta_loc,w_ion);
        x_ion_all( i, : ) = x_ion;
    end
end

x_pl = x_pl_all( 1, : );
for i=2:length(lambda_pl)
    x_pl = x_pl + lambda_pl(i)*(x_pl_all(i,:)-x_pl_all(1,:));
end
    
x_in_loc = x_pl;
for i=1:length(J_ion)
   x_in_loc = x_in_loc + lambda_ion(i)*(x_ion_all( i, : )-x_pl);
end

% x_in_loc = (x_pl + x_ion)/(sum(lambda_ion)+sum(lambda_pl));
x_in = interp1(delta_loc,x_in_loc,delta);
end