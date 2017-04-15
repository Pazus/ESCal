function decs = expandLegCoefsToCrossSection(x_l, mu)
assert(isvector(x_l),'x_l has to be a vector');

x_l = x_l(:)';
nLeg = numel(x_l)-1;
legCoef = (0:nLeg)+0.5;

P = Legendre_mu(mu,0,nLeg);  
decs = x_l.*legCoef*P;