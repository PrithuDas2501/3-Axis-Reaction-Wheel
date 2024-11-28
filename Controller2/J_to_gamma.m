function gamma = J_to_gamma(R)
% gamma = [J11,J22,J33,J23,J13,J12]
gamma = [R(1,1);R(2,2);R(3,3);R(2,3);R(1,3);R(1,2)];
end