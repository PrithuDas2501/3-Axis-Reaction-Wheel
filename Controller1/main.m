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
R0   = diag([1,1,1]);

Wlimit = 20;
Ulimit = 1000;

X0 = [1;-1;0.5;0;0;0;reshape(R0,[9,1])];
% U  = [0;0;0];
% [t,S] = ode45(@(t,S)Dynamics(t,S,U,J_sc,J_wh),tsim,S0);
[t,X] = ode45(@(t,X)Controller(t,X,J_sc,J_wh,Rd,Wlimit,Ulimit,Jbeta),tsim,X0);

%% Data Collection

n = length(t);

R = zeros(3,3,n);
S = zeros(3,n);
Snorm = zeros(1,n);
EigenError = zeros(1,n);
U = zeros(3,n);

I = eye(3);
A = diag([1,2,3]);
for i = 1:n
    r = [X(i,7:9)',X(i,10:12)',X(i,13:15)'];
    R(:,:,i) = r;
    Rbar = Rd'*r;
    w = X(i,1:3)';
    mu = X(i,4:6)';
    
    % eul = rotm2eul(r);

    e = 0.0;
    s = zeros(3,1);
    for j  = 1:3
        s = s + A(j,j)*cross(Rbar'*I(:,j),I(:,j)) + e;
    end
    % s = A*eul';
    S(:,i) = s;
    Snorm(i) = norm(s);

    EigenError(i) = acos(0.5*(trace(Rbar)-1));
    % if i>100
    %     SettleCheck = 1;
    %     for j = i-100:i
    %         if EigenError(j)>0.05
    %             SettleCheck = 0;
    %             break
    %         end
    %     end
    %     if SettleCheck==1
    %         SettleIndex = i;
    %         break
    %     end
    % end

    gamma = 5;
    nu    = 5;
    Kp    = gamma/trace(A);
    Kv    = diag([1/(1+abs(X(i,1))),1/(1+abs(X(i,2))),1/(1+abs(X(i,3)))]);
    % Kv    = nu*diag([1,1,1]);
    U(:,i) = J_wh\(Kp*s + Kv*X(i,1:3)');

    if abs(U(1,i))>Ulimit
    U(1,i) = sign(U(1))*Ulimit;
    end
    
    if abs(U(2,i))>Ulimit
        U(2,i) = sign(U(2))*Ulimit;
    end
    
    if abs(U(3,i))>Ulimit
        U(3,i) = sign(U(3))*Ulimit;
    end

    if abs(mu(1))>Wlimit && mu(1)*U(1,i)>0
    U(1,i) = 0;
    end
    
    if abs(mu(2))>Wlimit && mu(2)*U(2,i)>0
        U(2,i) = 0;
    end
    
    if abs(mu(3))>Wlimit && mu(3)*U(3,i)>0
        U(3,i) = 0;
    end
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
plot(t,Snorm,'r',LineWidth=2)

figure(5); clf
grid on 
hold on
title('Theta Error Vs Time')
xlabel('Time (s)')
ylabel('Theta')
plot(t,EigenError,'r',LineWidth=2)