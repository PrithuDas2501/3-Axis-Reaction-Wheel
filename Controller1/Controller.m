function dX = Controller(t,X,J_sc,J_wh,Rd,Wlimit,Ulimit,Jbeta)

dX = zeros(length(X),1);

[W, R] = UnFold(X);
w = W(1:3);
wstar = [0,   -w(3) w(2);
         w(3), 0,   -w(1);
        -w(2), w(1), 0];
Rdot = R*wstar;
Rbar = Rd'*R;

I = eye(3);
A = diag([1,2,3]);

e = 0.0;
S = zeros(3,1);
for i  = 1:3
    S = S + A(i,i)*cross(Rbar'*I(:,i),I(:,i)) + e;
end
% S = A*eul';
gamma = 5;
nu    = 5;
Kp    = gamma/trace(A);
Kv    = nu*diag([1/(1+abs(w(1))),1/(1+abs(w(2))),1/(1+abs(w(3)))]);
% Kv    = nu*diag([1,1,1]);
U = J_wh\(Kp*S + Kv*w);

dW = Dynamics(0,W,U,J_sc,J_wh,Wlimit,Ulimit,Jbeta);

dX(1:6) = dW;
dX(7:15) = reshape(Rdot,[9,1]);
end