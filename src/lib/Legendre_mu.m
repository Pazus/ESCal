function P = Legendre_mu(mu,m,L)
    smu = sqrt(1-mu.^2);
%               Legendre polynomial - Plm
    N = length(mu);
    m2 = m^2;
    P = ones(L - m, N);
    if m > 0    %               Compute Pmm
        im2 = 2: 2: 2*m;
        sim2 = realsqrt(prod(1 - 1 ./ im2));
        P(1,:) = sim2 * smu .^ m;
    end
%                               Compute Pmm+1
    sim1 = realsqrt(2*m+1);
    P(2,:) = sim1 * mu .* P(1,:);
%                               Compute Pml, l > m + 1
    for k = (m + 2): L
        i = k - m;
%         if m==0
%             k1 = (2*k-1)/k;
%             k2 = (k-1)/k;
%         else
%             iskm = 1/realsqrt(k^2-m2);
%             k1 = (2*k-1)*iskm;
%             k2 = realsqrt((k-1)^2-m2)*iskm;
%         end
        iskm = 1/realsqrt(k^2-m2);
        k1 = (2*k-1)*iskm;
        k2 = realsqrt((k-1)^2-m2)*iskm;
        P(i+1,:) = k1* mu .* P(i,:) - k2*P(i-1,:);      
    end    
end