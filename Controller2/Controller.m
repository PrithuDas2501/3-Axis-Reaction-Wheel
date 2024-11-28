function dX = Controller(t,X,J_sc,J_wh,Rd,wd,wddot,tau_dist,Wlimit,Ulimit,Jbeta)


dX = zeros(length(X),1);

[W, R] = UnFold(X);
w = W(1:3);
mu = W(4:6);
wstar = [0,   -w(3) w(2);
         w(3), 0,   -w(1);
        -w(2), w(1), 0];
Rdot = R*wstar;
Rbar = Rd'*R;
wbar = w - Rbar'*wd;

I = eye(3);
A = diag([1,2,3]);

e = 0.0;
S = zeros(3,1);
Sdot = zeros(3,1);
for i  = 1:3
    S = S + A(i,i)*cross(Rbar'*I(:,i),I(:,i)) + e;
    Sdot = Sdot + A(i,i)*cross(cross(Rbar'*I(:,i),wbar),I(:,i)) + e;
end

K1 = I;
gamma = 1;
nu    = 1;
Kp    = gamma/trace(A);
Kv    = nu*diag([1/(1+abs(w(1))),1/(1+abs(w(2))),1/(1+abs(w(3)))]);

J_sc_hat = gamma_to_J(X(16:21));

tau_dist_hat = X(22:24);

v1 = -cross(J_sc_hat*w + J_wh*mu,w) - J_sc_hat*(K1*Sdot + cross(wbar,w) - Rbar'*wddot);
v2 = -tau_dist_hat;
v3 = -Kv*(wbar + K1*S) - Kp*S;

U = -J_wh\(v1+v2+v3);

dW = Dynamics(0,W,U,J_sc,J_wh,tau_dist,Wlimit,Ulimit,Jbeta);

dX(1:6) = dW;
dX(7:15) = reshape(Rdot,[9,1]);

Q = eye(6);
dX(16:21) = Q\(w_to_L(w)'*wstar + w_to_L(K1*Sdot + cross(wbar,w) - Rbar'*wddot)')*(wbar + K1*S);

A_d = zeros(3);
D   = eye(3);
dX(22:24) = A_d*tau_dist_hat + D\(wbar + K1*S);

end