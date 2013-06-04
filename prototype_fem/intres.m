function [ output_args ] = intres( fem, d, Res, s, i, j )
%INTRES Summary of this function goes here
%   Detailed explanation goes here

x3 = [-0.774596669241483, ...
       0., ...
       0.774596669241483];
w3 = [ 0.555555555555555, ...
       0.888888888888888, ...
       0.555555555555555];

for a = 1:4
    output_args = output_args + w3(a)/2 * Res( fem, d, s, (x3(a)+1)/2, i, j );
end
   
end

