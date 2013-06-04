function [ output_args ] = solver( fem, d )
%SOLVER Summary of this function goes here
%   Detailed explanation goes here

% R doesn't change with s, so just calculate it once
R = zeros(6,1);

s = INVLAPs(fem.tinit, fem.tend, fem.dt);

u = zeros(6,1); %Lazy way to initialize

for i = 1:length(s)
    J = AssembleJ(fem, d, s(i));
    ui = J\R; %Funky matrix division
    
    u = [u, ui]; %Add the solution at this frequency to the others.
end

u(:,1) = []; %Get rid of the hax to make the rest of the stuff easier.

ut = zeroes(1:nnt);

for i = 1:6
    ut = [ut;INVLAP(u(i,:), tini, tend, nnt)];
end

output_args = ut;

end

