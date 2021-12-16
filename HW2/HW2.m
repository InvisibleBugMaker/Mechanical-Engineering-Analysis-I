%Zhaoyi Jiang
%ME564 HW2
%----------------------------------------------------------------------
%P1
clc
clear all
close all

A=[0,1;-9,-10];
[V,D]=eig(A);
C=[0.2;0];
t=0:0.5:10;
Y=zeros(2,length(t));
for ii=1:length(t)
Y(:,ii)=V*[exp(D(1,1)*t(ii)),0;0,exp(D(2,2)*t(ii))]/V*C;
end
tspan=[0,10];
[t_nu,Y_nu]=ode45(@(t,Y)A*Y,tspan,C);

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.01 0.2]);
fplot(@(t)0.225*exp(-t)-0.025*exp(-9*t))
plot(t,Y(1,:),'x')
plot(t_nu,Y_nu(:,1),'b.')
title('P1 Displacement');
legend({'2nd order ODE','System of ODEs','ode45'});
xlabel('t');ylabel('y(t)');

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.16 0.01]);
fplot(@(t)-0.225*exp(-t)+0.225*exp(-9*t))
plot(t,Y(2,:),'x')
plot(t_nu,Y_nu(:,2),'b.')
title('P1 Velocity');
legend('2nd order ODE','system of ODEs','ode45','location','southeast');
xlabel('t');ylabel('dy/dt');

%----------------------------------------------------------------------
%P2
A=[0,1;-9,-6];
D=eig(A);
[H,J]=jordan(A);
C=[0.2;0];
t=0:0.5:10;
tspan=[0,10];
Y=zeros(2,length(t));
for ii=1:length(t)
Y(:,ii)=H*[exp(D(1,1)*t(ii)),t(ii)*exp(D(1,1)*t(ii));0,exp(D(1,1)*t(ii))]/H*C;
end
[t_nu,Y_nu]=ode45(@(t,Y)A*Y,tspan,C);

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.01 0.2]);
fplot(@(t)0.2*exp(-3*t)+0.6*t*exp(-3*t))
plot(t,Y(1,:),'x')
plot(t_nu,Y_nu(:,1),'b.')
title('P2 Displacement');
legend({'2nd order ODE','System of ODEs','ode45'});
xlabel('t');ylabel('y(t)');

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.23 0.01]);
fplot(@(t)-1.8*t*exp(-3*t))
plot(t,Y(2,:),'x')
plot(t_nu,Y_nu(:,2),'b.')
title('P2 Velocity');
legend('2nd order ODE','system of ODEs','ode45','location','southeast');
xlabel('t');ylabel('dy/dt');

%----------------------------------------------------------------------
%P3
A=[0,1;-9,-1];
[V,D]=eig(A);
C=[0.2;0];
t=0:0.5:10;
Y=zeros(2,length(t));
for ii=1:length(t)
Y(:,ii)=V*[exp(D(1,1)*t(ii)),0;0,exp(D(2,2)*t(ii))]/V*C;
end
tspan=[0,10];
[t_nu,Y_nu]=ode45(@(t,Y)A*Y,tspan,C);

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.13 0.2]);
fplot(@(t)exp(-0.5*t)*(0.034*sin(2.96*t)+0.2*cos(2.96*t)))
plot(t,Y(1,:),'x')
plot(t_nu,Y_nu(:,1),'b.')
title('P3 Displacement');
legend('2nd order ODE','System of ODEs','ode45');
xlabel('t');ylabel('y(t)');

figure; 
hold on;
grid on;
xlim(tspan);
ylim([-0.5 0.3]);
fplot(@(t)-0.609*exp(-0.5*t)*sin(2.96*t)+0.00064*exp(-0.5*t)*cos(2.96*t))
plot(t,Y(2,:),'x')
plot(t_nu,Y_nu(:,2),'b.')
title('P3 Velocity');
legend('2nd order ODE','system of ODEs','ode45');
xlabel('t');ylabel('dy/dt');
