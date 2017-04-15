function [x,w]=legzo_n1(n)
%       =========================================================
%       Purpose : Compute the zeros of Legendre polynomial Pn(x)
%                 in the interval [-1,1], and the corresponding
%                 weighting coefficients for Gauss-Legendre
%                 integration
%       Input :   n    --- Order of the Legendre polynomial
%       Output:   X(n) --- Zeros of the Legendre polynomial
%                 W(n) --- Corresponding weighting coefficients
%       =========================================================
if strcmp(computer('arch'),'win64') 
    [x,w] = legzo_n1_mex(n);
else
    x=zeros(1,n);
    w=zeros(1,n);
    eps=1.0e-15;
    m = (n+1)/2;
    for ii=1:m
        z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
        z1 = z+1;
        while abs(z-z1)>eps
            p1 = 1;
            p2 = 0;
            for jj = 1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
            end
            pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
            z1 = z;
            z = z1-p1/pp;
        end
        x(ii) = -z; % Build up the abscissas.
        x(n+1-ii) = z;
        w(ii) = 2/((1-z^2)*(pp^2)); % Build up the weights.
        w(n+1-ii) = w(ii);
    end
    x=-x;
end