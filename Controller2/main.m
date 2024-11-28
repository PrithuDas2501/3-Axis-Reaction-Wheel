clear; clc

%% Simulation

T_end = 100;
tsim  = linspace(0,T_end,T_end*1000);

beta1 = 0.375;
beta2 = 0.375;
beta3 = 0.375;
Jbeta = diag([beta2+beta3,beta1+beta3,beta1+beta2]);

J_sc = diag([10,25/3,5]);
J_wh = diag([0.5,0.5,0.5]);
Rd   = diag([1,-1,-1]);
% Rd = diag([1,1,1]);
wd   = 0*[0.5;-0.5;-0.3];
wddot= [0;0;0];
R0   = diag([1,1,1]);

tau_dist = 1*[0.7;-0.3;0];

Wlimit = 100;
Ulimit = 1000;

gamma0hat = R_to_gamma(J_wh);
gamma_actual = R_to_gamma(J_sc+Jbeta);

tau_dist_hat = zeros(3,1);

X0 = [1;-1;0.5;0;0;0;reshape(R0,[9,1]);gamma0hat;tau_dist_hat];
[t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,wd,wddot,tau_dist,Wlimit,Ulimit,Jbeta),tsim,X0);

%% Data Collection

n = length(t);

R = zeros(3,3,n);
S = zeros(3,n);
Sdot = zeros(3,n);
Snorm = zeros(1,n);
EigenError = zeros(1,n);
tau_dist_error = zeros(3,n);
gamma_error = zeros(6,n);
U = zeros(3,n);
I = eye(3);
A = diag([1,2,3]);
for i = 1:n
    r = [X(i,7:9)',X(i,10:12)',X(i,13:15)'];
    R(:,:,i) = r;
    Rbar = Rd'*r;
    w = [X(i,1);X(i,2);X(i,3)];
    mu = [X(i,4);X(i,5);X(i,6)];
    wbar = w - Rbar'*wd;
    
    e = 0.0;
    s = zeros(3,1);
    sdot = zeros(3,1);

    for j  = 1:3
        s = s + A(j,j)*cross(Rbar'*I(:,j),I(:,j)) + e;
        sdot = sdot + A(j,j)*cross(cross(Rbar'*I(:,j),wbar),I(:,j)) + e;
    end

    S(:,i) = s;
    Snorm(i) = norm(s);
    Sdot(:,i) = sdot;
    
    K1 = I;
    gamma = 1;
    nu    = 1;
    Kp    = gamma/trace(A);
    Kv    = nu*diag([1/(1+abs(w(1))),1/(1+abs(w(2))),1/(1+abs(w(3)))]);

    J_sc_hat = gamma_to_J(X(i,16:21));
    gamma_error(:,i) = gamma_actual - X(i,16:21)';

    tau_dist_hat = X(i,22:24)';
    tau_dist_error(:,i) = tau_dist - tau_dist_hat;
    
    v1 = -cross(J_sc_hat*w + J_wh*mu,w) - J_sc_hat*(K1*sdot + cross(wbar,w) - Rbar'*wddot);
    v2 = -tau_dist_hat;
    v3 = -Kv*(wbar + K1*s) - Kp*s;
    
    U(:,i) = -J_wh\(v1+v2+v3);

    EigenError(i) = acos(0.5*(trace(Rbar)-1));

end


Wx_sc = X(:,1);
Wy_sc = X(:,2);
Wz_sc = X(:,3);
Wx_wh = X(:,4);
Wy_wh = X(:,5);
Wz_wh = X(:,6);

%% Plotting

figure(1); clf
grid on 
hold on
legend('on')
title('Body Rotation Rate Vs Time')
xlabel('Time (s)')
ylabel('Rotation Rates (Rad/s)')
plot(t,Wx_sc,'r',LineWidth=2,DisplayName='Wx')
plot(t,Wy_sc,'g',LineWidth=2,DisplayName='Wy')
plot(t,Wz_sc,'b',LineWidth=2,DisplayName='Wz')

figure(2); clf
grid on 
hold on
legend('on')
title('Wheel Rotation Rate Vs Time')
xlabel('Time (s)')
ylabel('Rotation Rates (Rad/s)')
plot(t,Wx_wh,'r',LineWidth=2,DisplayName='Wx^w')
plot(t,Wy_wh,'g',LineWidth=2,DisplayName='Wy^w')
plot(t,Wz_wh,'b',LineWidth=2,DisplayName='Wz^w')


figure(3); clf
grid on 
hold on
legend('on')
title('Control Input Vs Time')
xlabel('Time (s)')
ylabel('Control Input (Rad/s^2)')
plot(t,U(1,:),'r',LineWidth=2,DisplayName='\alpha_x')
plot(t,U(2,:),'g',LineWidth=2,DisplayName='\alpha_y')
plot(t,U(3,:),'b',LineWidth=2,DisplayName='\alpha_z')


figure(4); clf
grid on 
hold on
title('|S|^2 Vs Time')
xlabel('Time (s)')
ylabel('|S|^2')
plot(t,Snorm.^2,'r',LineWidth=2)

figure(5); clf
grid on 
hold on
title('Theta Error Vs Time')
xlabel('Time (s)')
ylabel('Theta')
plot(t,EigenError,'r',LineWidth=2)

figure(6); clf
grid on 
hold on
legend('on')
title('Disturbance Estimate Error Vs Time')
xlabel('Time (s)')
ylabel('Disturbance Estimate Error')
plot(t,tau_dist_error(1,:),'r',LineWidth=2,DisplayName='\tau_x')
plot(t,tau_dist_error(2,:),'g',LineWidth=2,DisplayName='\tau_y')
plot(t,tau_dist_error(3,:),'b',LineWidth=2,DisplayName='\tau_z')

figure(7); clf
grid on 
hold on
legend('on')
title('Inertia Estimate Error Vs Time')
xlabel('Time (s)')
ylabel('Inertia Estimate Error')
plot(t,gamma_error(1,:),LineWidth=1,DisplayName='J_{11}')
plot(t,gamma_error(2,:),LineWidth=1,DisplayName='J_{22}')
plot(t,gamma_error(3,:),LineWidth=1,DisplayName='J_{33}')
plot(t,gamma_error(4,:),LineWidth=1,DisplayName='J_{23}')
plot(t,gamma_error(5,:),LineWidth=1,DisplayName='J_{13}')
plot(t,gamma_error(6,:),LineWidth=1,DisplayName='J_{12}')