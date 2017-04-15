function [x_p_m, x_m_m] = expandLegCoefsTo2ParamCrossSection(x_l, mu, m)
if nargin < 3
    m=0;
end
calcOpposite = nargout == 2;

assert(isvector(x_l),'x_l has to be a vector');
assert(isvector(mu),'theta has to be a vector');
assert(isscalar(m),'m has to be a scalar');

x_l = x_l(:);
nLeg = numel(x_l)-1;
nMu = numel(mu);
l_ = (0:nLeg)';
legCoef = l_+0.5;

x_p_m = zeros(nMu,nMu,m+1);
x_m_m = zeros(nMu,nMu,m+1);

xlp = x_l.*legCoef;
xlm = (-1).^l_.*xlp;


for mi=0:m
    P = Legendre_mu(mu(:)',mi,nLeg);  
    x_p_m(:,:,mi+1) = symmetrize(P'*spdiags(xlp(mi+1:end),0,nLeg-mi+1,nLeg-mi+1)*P);
    
    if calcOpposite
        x_m_m(:,:,mi+1) = symmetrize(P'*spdiags(xlm(mi+1:end),0,nLeg-mi+1,nLeg-mi+1)*P);
    end
end

