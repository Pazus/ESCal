function M = optimal_M(theta)
M = floor(0.41*theta-6.*10^-3*theta.^2+1.75*10^-10*theta.^6+0.8);