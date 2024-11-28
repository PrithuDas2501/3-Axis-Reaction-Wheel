clear; clc
%% Simulation Set-Up

T_end = 100;
tsim  = linspace(0,T_end,T_end*100);

beta1 = 0.375;
beta2 = 0.375;
beta3 = 0.375;
Jbeta = diag([beta2+beta3,beta1+beta3,beta1+beta2]);

J1 = diag([10,10,10]);  % Sphere
J2 = diag([10,10,5]);   % Cylinder
J3 = diag([10,25/3,5]); % Centroid
J4 = diag([10,5,5]);    % Thin Disk
J5 = diag([10,10,0.1]); % Thin Cylinder 
J_wh = diag([0.5,0.5,0.5]);

tau_dist = 1*[0.7;-0.3;0];

gamma0hat = R_to_gamma(J_wh);
gamma_actual = R_to_gamma(J3+Jbeta);

tau_dist_hat = zeros(3,1);

nl = 20;
lambda = linspace(0,1,nl);

Rd   = diag([1,-1,-1]);
wd   = [0;0;0];
wddot= [0;0;0];
R0   = diag([1,1,1]);

Wlimit = 1000;
Ulimit = 1000;

X0 = [1;-1;0.5;0;0;0;reshape(R0,[9,1]);gamma0hat;tau_dist_hat];

%% To Thin Disk
Ts1 = zeros(nl,1);
for L = 1:nl
    J_sc = (1-lambda(L))*J3 + lambda(L)*J4;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,wd,wddot,tau_dist,Wlimit,Ulimit,Jbeta),tsim,X0);
    
    n = length(t);
    EigenError = zeros(1,n);
    for i = 1:n
        r = [X(i,7:9)',X(i,10:12)',X(i,13:15)'];
        Rbar = Rd'*r;
    
        EigenError(i) = acos(0.5*(trace(Rbar)-1));
        if i>100
            SettleCheck = 1;
            for j = i-100:i
                if EigenError(j)>=0.05
                    SettleCheck = 0;
                    break
                end
            end
            if SettleCheck==1
                Ts1(L) = t(i);
                break
            end
        end
    end
end

%% To Sphere
Ts2 = zeros(nl,1);
for L = 1:nl
    J_sc = (1-lambda(L))*J3 + lambda(L)*J1;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,wd,wddot,tau_dist,Wlimit,Ulimit,Jbeta),tsim,X0);
    
    n = length(t);
    EigenError = zeros(1,n);
    for i = 1:n
        r = [X(i,7:9)',X(i,10:12)',X(i,13:15)'];
        Rbar = Rd'*r;
    
        EigenError(i) = acos(0.5*(trace(Rbar)-1));
        if i>100
            SettleCheck = 1;
            for j = i-100:i
                if EigenError(j)>=0.05
                    SettleCheck = 0;
                    break
                end
            end
            if SettleCheck==1
                Ts2(L) = t(i);
                break
            end
        end
    end
end

%% To Thin Cylinder
Ts3 = zeros(nl,1);
for L = 1:nl
    J_sc = (1-lambda(L))*J3 + lambda(L)*J5;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,wd,wddot,tau_dist,Wlimit,Ulimit,Jbeta),tsim,X0);
    
    n = length(t);
    EigenError = zeros(1,n);
    for i = 1:n
        r = [X(i,7:9)',X(i,10:12)',X(i,13:15)'];
        Rbar = Rd'*r;
    
        EigenError(i) = acos(0.5*(trace(Rbar)-1));
        if i>100
            SettleCheck = 1;
            for j = i-100:i
                if EigenError(j)>=0.05
                    SettleCheck = 0;
                    break
                end
            end
            if SettleCheck==1
                Ts3(L) = t(i);
                break
            end
        end
    end
end

%% Plot
figure(1); clf
grid on 
hold on
legend('on')
title('Settling Time Vs \lambda')
xlabel('\lambda')
ylabel('Settling Time (s)')
plot(lambda,Ts1,'r',LineWidth=2,DisplayName='Centroid to Thin Disk')
plot(lambda,Ts2,'g',LineWidth=2,DisplayName='Centroid to Sphere')
plot(lambda,Ts3,'b',LineWidth=2,DisplayName='Centroid to Thin Cylinder')