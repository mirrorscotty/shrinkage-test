function [ output_args ] = dradp( fem, d, s, x, i, j )

alpha = d.alpha(s);
beta = d.beta(s);

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

output_args = (alpha - beta)*phij*dphii;

end

