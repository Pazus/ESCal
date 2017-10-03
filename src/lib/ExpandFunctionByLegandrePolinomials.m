function y_l = ExpandFunctionByLegandrePolinomials(mu_mesh,values,LegCount,w)
PP = Legendre_mu(mu_mesh,0,LegCount);

if nargin<4
    w = 1;
end
y_l=(values.*w)*PP';
y_l = y_l/y_l(1);