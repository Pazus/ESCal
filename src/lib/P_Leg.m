function P = P_Leg(mu,N)
s = size(mu);
if (s(1) == 1) && (s(2) > 1)
    mu = mu';
elseif (s(1) > 1) && (s(2) > 1)
    P = 0;
    return;
end
P = ones(numel(mu),N);
P(:,2) = mu;
l = 0:N-1;
for i=2:N-1
    P(:,i+1) = ( (2*l(i)+1)*mu.*P(:,i)-l(i)*P(:,i-1) )/(l(i)+1);
end