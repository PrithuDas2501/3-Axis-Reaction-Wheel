function dS = Dynamics(t,S,U,J_sc,J_wh,tau_dist,Wlimit,Ulimit,Jbeta)

% S = [1Wx_sc; 2Wy_sc; 3Wy_sc; 4Wx_wh; 5Wy_wh; 6Wz_wh] 
% U = [1Wxdot_wh; 2Wydot_wh; 3Wzdot_wh]

w  = [S(1);S(2);S(3)];
mu = [S(4);S(5);S(6)];

% Control Input Limits
if abs(U(1))>Ulimit
    U(1) = sign(U(1))*Ulimit;
end

if abs(U(2))>Ulimit
    U(2) = sign(U(2))*Ulimit;
end

if abs(U(3))>Ulimit
    U(3) = sign(U(3))*Ulimit;
end

% Wheel Rotation Rate Limits
if abs(mu(1))>Wlimit && mu(1)*U(1)>0
    U(1) = 0;
end

if abs(mu(2))>Wlimit && mu(2)*U(2)>0
    U(2) = 0;
end

if abs(mu(3))>Wlimit && mu(3)*U(3)>0
    U(3) = 0;
end

wdot = (J_sc+Jbeta)\(cross((J_sc+Jbeta+J_wh)*w + J_wh*mu,w) - J_wh*U + tau_dist);

dS = [wdot;U];
end