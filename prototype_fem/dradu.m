function [ output_args ] = dradu( fem, d, s, x, i, j )

g = d.ghat(s);
k = d.khat(s);
rho = d.rho(s);
beta = d.beta(s);
rhof = d.rhof;

if(i)
    phii = fem.phi1(x);
    dphii = fem.dphi1(0);
else
    phii = fem.phi0(x);
    dphii = fem.dphi0(x);
end

if(j)
    phij = fem.phi1(x);
    dphij = fem.dphi1(x);
else
    phij = fem.phi0(x);
    dphij = fem.dphi0(x);
end

output_args = g*dphii*dphij + (k+g/3)*dphii*dphij - s*(rho-beta*rhof)*phii*phij;

end
