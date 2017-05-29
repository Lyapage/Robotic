function quadcopter
%-------------------------------------------------------------------
close all;
clear all;
global a b c K
%-------------------------------------------------------------------
%system matrix A,B, since the complicated expression for A
%there will be row vectors of A first, and then whole matrix
a1 = [0 0 0 1 0 0 0 0 0 0 0 0];
a2 = [0 0 0 0 1 0 0 0 0 0 0 0];
a3 = [0 0 0 0 0 1 0 0 0 0 0 0];
a4 = [0 0 0 -0.5342 0 0 0 0 0 0 0 0];
a5 = [0 0 0 0 -0.5342 0 0 -9.81 0 0 0 0];
a6 = [0 0 0 0 0 -0.5342 -2.45 0 0 0 0 0];
a7 = [0 0 0 0 0 0 0 0 0 1 0 0];
a8 = [0 0 0 0 0 0 0 0 0 0 1 0];
a9 = [0 0 0 0 0 0 0 0 0 0 0 1];
a10 = [0 0 0 0 0 0 0 0 0 0 0 0];
a11 = [0 0 0 0 0 0 0 0 0 0 0 0];
a12 = [0 0 0 0 0 0 0 0 0 0 0 0];
A = [a1;a2;a3;a4;a5;a6;a7;a8;a9;a10;a11;a12];

a = 0.0637;
b = 0.14;
c = 0.529;

B = [0 0 0 0;0 0 0 0;0 0 0 0;
     a 0 0 0;0 0 0 0;a a a a;
     0 0 0 0;0 0 0 0;0 0 0 0;
     b 0 -b 0; %line 10
     0 b 0 -b; %line 11
     c -c c -c]; %line 12
%-------------------------------------------------------------------
%check the controllability of quadcopter system:
C = ctrb(A,B);      %define the controllability matrix
r = rank(C);        %check the rank controllability matrix
if r >= 12
    disp('This system is controllable.') 
end
%-------------------------------------------------------------------
%check the stability of system through eigenvalues:
eig_A = eig(A);
syms s 
minpoly(A,s);
%-------------------------------------------------------------------
%quadratic cost function parameter:
Q = eye(12);
R = eye(4);

K = lqr(A,B,Q,R);      %obtaining optimal feedback gain through LQR
%-------------------------------------------------------------------
% check the stability after lqr applied:
eig_A_BK = eig(A-B*K);
if real(eig_A_BK) < 0
   disp('After LQR applied,the system is stabilized'); 
end

%-------------------------------------------------------------------
%initial condition:
x1_0 = [0 1 0]';
x2_0 = [0 0 0]';
x3_0 = [0 0 pi/4]';
x4_0 = [0 0 0]';
X0 = [x1_0;x2_0;x3_0;x4_0];
%-------------------------------------------------------------------
%using ode45 to obtain   
N=501;
t=linspace(0,10,N);
[t X] = ode45(@eom,t,X0);
    
%-------------------------------------------------------------------
%data processing:
for i = 1:N
    x(i) = X(i,1);              %x-position
    y(i) = X(i,2);              %y-position
    z(i) = X(i,3);              %z-position
    vx(i) = X(i,4);             %velocity in x-direction
    vy(i) = X(i,5);             %velocity in y-direction
    vz(i) = X(i,6);             %velocity in z-direction
    phi(i) = X(i,7);            %roll angle
    theta(i) = X(i,8);          %pitch angle
    psi(i) = X(i,9);            %yaw angle
    wx(i) = X(i,10);            %angular velocity in x-direction
    wy(i) = X(i,11);            %angular velocity in y-direction
    wz(i) = X(i,12);            %angular velocity in z-direction
    u(i,1:4) = -K*X(i,1:12)';   %four input: Omega_i^2 
end

%-------------------------------------------------------------------
%mapping for each state variables:
figure;
plot(t,x,t,y,'--',t,z,':','linewidth',8);
xlabel('t');
ylabel('Position');
legend('x','y','z','location','northeast');

figure;
plot(t,vx,t,vy,'--',t,vz,':','linewidth',8);
xlabel('t');
ylabel('Velocity');
legend('v_x','v_y','v_z','location','southeast');


figure;
plot(t,phi,t,theta,'--',t,psi,':','linewidth',8);
xlabel('t');
ylabel('Three euler angles');
legend('\phi','\theta','\psi','location','northeast');


figure;
plot(t,wx,t,wy,'--',t,wz,':','linewidth',8);
xlabel('t');
ylabel('Angular velocity');
legend('\omega_x','\omega_y','\omega_z','location','northeast');

figure;
plot(t,u(1:501,1),t,u(1:501,2),'--',t,u(1:501,3),':',t,u(1:501,4),'*'...
    ,'linewidth',8);
xlabel('t');
ylabel('Four input');
legend('u_1','u_2','u_3',...
    'u_4','location','northeast');


save quadcopter;
evalin('base','load quadcopter');
end

%-------------------------------------------------------------------
%constructing the equation of motion:
function X_dot = eom(t,X)
global a b c 
x1 = X(1);x2 = X(2);x3 = X(3);x4 = X(4);x5 = X(5);x6 = X(6);
x7 = X(7);x8 = X(8);x9 = X(9);x10 = X(10);x11 = X(11);x12 = X(12);

u = control(t,X);

x1_dot = x4;x2_dot = x5;x3_dot = x6;x4_dot = -0.5342*x4+a*u(1);
x5_dot = -0.5342*x5-9.81*x8;x6_dot = -0.5432*x6-2.4525*x7;
x7_dot = x10;x8_dot = x11;x9_dot = x12;
x10_dot = b*u(1)-b*u(3);
x11_dot = b*u(2)-b*u(4);
x12_dot = c*u(1)-c*u(2)+c*u(3)-c*u(4);

X_dot = [x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot;x7_dot;x8_dot;x9_dot;...
    x10_dot;x11_dot;x12_dot];

end

%-------------------------------------------------------------------
%constructing input u:
function u = control(t,X)
global K
u = -K*X;
end