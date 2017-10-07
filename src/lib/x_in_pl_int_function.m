function x_in_pl = x_in_pl_int_function(delta, E0, epsilon_pl, b, alpha, beta, N_k)
% Интегральное представление неупругого сечения
if nargin<7; N_k = 10000; end %7000
if nargin<5; alpha=0; beta= 0; end

x_in_pl=zeros(1,numel(delta)); 
% cc = zeros(numel(delta),N_k+1);
for i = 2:numel(delta)
    delta_loc = delta(i);
    k_minus = sqrt(E0)-sqrt(E0-delta_loc);
    k_plus  = sqrt(E0)+sqrt(E0-delta_loc);
    step = (k_plus-k_minus)/N_k;
    k = k_minus:step:k_plus;
%     plot(delta,x_in_pl)
    res = x_in_pl_k(k,delta_loc, epsilon_pl, b, alpha, beta);
%     cc(i,:) = res./k; 
    x_in_pl(i) = trapz(k,res./k); %sum(res./k)*step;
%     if delta(i)<50 || 1; plot3(ones(1,length(k(:)))*delta(i), k(:), res./k); hold on; end
end
% hold off;xlabel('\Delta, eV');ylabel('k_{norm}') 
end

function res = x_in_pl_k(k,delta,epsilon_pl,b, alpha, beta)
res = delta^beta./(((epsilon_pl + k.^2).^2 - delta^2 ).^2+delta^alpha*b^(4-alpha));
end