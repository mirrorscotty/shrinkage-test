function [ output_args ] = drbdp( fem, d, s, x, i, j )

beta = d.beta(s);
rhof = d.rhof;
poro = d.poro;
r = d.r(s);

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

output_args = beta/s/rhof*dphii*phij - poro^2*s/r*phii*phij;

end

