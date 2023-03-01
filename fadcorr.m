function Omega = fadcorr(N,rho)

[x y] = meshgrid(1:N);
d = abs(x-y);
r = rho.^d;
Omega = chol(r)';