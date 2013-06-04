function [ J ] = AssembleJ( fem, d, s )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

%First Element
J(1,1) = intres(fem, d, @dradu, s, 0, 0);
J(1,2) = intres(fem, d, @dradp, s, 0, 0);
J(2,1) = intres(fem, d, @drbdu, s, 0, 0);
J(2,2) = intres(fem, d, @drbdp, s, 0, 0);

J(1,3) = intres(fem, d, @dradu, s, 1, 0);
J(1,4) = intres(fem, d, @dradp, s, 1, 0);
J(2,3) = intres(fem, d, @drbdu, s, 1, 0);
J(2,4) = intres(fem, d, @drbdp, s, 1, 0);

J(3,1) = intres(fem, d, @dradu, s, 0, 1);
J(3,2) = intres(fem, d, @dradp, s, 0, 1);
J(4,1) = intres(fem, d, @drbdu, s, 0, 1);
J(4,2) = intres(fem, d, @drbdp, s, 0, 1);

J(3,3) = intres(fem, d, @dradu, s, 1, 1);
J(3,4) = intres(fem, d, @dradp, s, 1, 1);
J(4,3) = intres(fem, d, @drbdu, s, 1, 1);
J(4,4) = intres(fem, d, @drbdp, s, 1, 1);

%Second Element
J(3,3) = J(3,3) + intres(fem, d, @dradu, s, 0, 0);
J(3,4) = J(3,4) + intres(fem, d, @dradp, s, 0, 0);
J(4,3) = J(4,3) + intres(fem, d, @drbdu, s, 0, 0);
J(4,4) = J(4,4) + intres(fem, d, @drbdp, s, 0, 0);

J(3,5) = intres(fem, d, @dradu, s, 1, 0);
J(3,6) = intres(fem, d, @dradp, s, 1, 0);
J(4,5) = intres(fem, d, @drbdu, s, 1, 0);
J(4,6) = intres(fem, d, @drbdp, s, 1, 0);

J(5,3) = intres(fem, d, @dradu, s, 0, 1);
J(5,4) = intres(fem, d, @dradp, s, 0, 1);
J(6,3) = intres(fem, d, @drbdu, s, 0, 1);
J(6,4) = intres(fem, d, @drbdp, s, 0, 1);

J(5,5) = intres(fem, d, @dradu, s, 1, 1);
J(5,6) = intres(fem, d, @dradp, s, 1, 1);
J(6,5) = intres(fem, d, @drbdu, s, 1, 1);
J(6,6) = intres(fem, d, @drbdp, s, 1, 1);

end

