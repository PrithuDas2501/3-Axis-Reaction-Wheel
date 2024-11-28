function R = gamma_to_R(gamma)
% gamma = [J11,J22,J33,J23,J13,J12]
R = [gamma(1), gamma(6), gamma(5);
     gamma(6), gamma(2), gamma(4);
     gamma(5), gamma(4), gamma(3)];
end