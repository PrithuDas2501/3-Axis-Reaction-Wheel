clear; clc
%% Simulation Set-Up

T_end = 100;
tsim  = linspace(0,T_end,T_end*100);

J1 = diag([10,10,10]); % Sphere
J2 = diag([10,10,5]);  % Cylinder
J3 = diag([10,25/3,5]);% Centroid
J4 = diag([10,5,5]);   % Thin Disk
J5 = diag([10,10,0.1]);% Thin Cylinder 
J_wh = diag([0.5,0.5,0.5]);

beta1 = 0.375;
beta2 = 0.375;
beta3 = 0.375;
Jbeta = diag([beta2+beta3,beta1+beta3,beta1+beta2]);

nl = 20;
theta = linspace(-pi,pi,nl);

Rd   = diag([1,-1,-1]);
R0   = diag([1,1,1]);

Wlimit = 1000;
Ulimit = 1000;

X0 = [1;-1;0.5;0;0;0;reshape(R0,[9,1])];

%% X
Ts1 = zeros(nl,1);
for L = 1:nl
    Rotmat = [1,0,0;0,cos(theta(L)),sin(theta(L));0,-sin(theta(L)),cos(theta(L))];
    J_sc = Rotmat'*J3*Rotmat;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,Wlimit,Ulimit,Jbeta),tsim,X0);
    
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

%% Y
Ts2 = zeros(nl,1);
for L = 1:nl
    Rotmat = [cos(theta(L)),0,-sin(theta(L));0,1,0;sin(theta(L)),0,cos(theta(L))];
    J_sc = Rotmat'*J3*Rotmat;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,Wlimit,Ulimit,Jbeta),tsim,X0);
    
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

%% Z
Ts3 = zeros(nl,1);
for L = 1:nl
    Rotmat = [cos(theta(L)),sin(theta(L)),0;-sin(theta(L)),cos(theta(L)),0;0,0,1];
    J_sc = Rotmat'*J3*Rotmat;
    [t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,Wlimit,Ulimit,Jbeta),tsim,X0);
    
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
title('Settling Time Vs Misalignment')
xlabel('\theta')
ylabel('Settling Time (s)')
plot(theta*180/pi,Ts1,'r',LineWidth=2,DisplayName='X')
plot(theta*180/pi,Ts2,'g',LineWidth=2,DisplayName='Y')
plot(theta*180/pi,Ts3,'b',LineWidth=2,DisplayName='Z')
