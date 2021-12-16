%ME564 HW5
%Zhaoyi Jiang
clear all
close all
clc
%% P2
A=[-10 10 0;
    28 -1 0;
    0 0 -8/3]
P2_b1=eig(A)

A=[-10 10 0;
    28 -1 -6*sqrt(2);
    6*sqrt(2) 6*sqrt(2) -8/3]
P2_b2=eig(A)

A=[-10 10 0;
    28 -1 6*sqrt(2);
    -6*sqrt(2) -6*sqrt(2) -8/3]
P2_b3=eig(A)


%% P3 b(1)
alpha=0;
beta=0;
y0=[pi/2;1];
dt=0.1;
T=20;
n=T/dt;
tspan=[0,T];
Y(:,1)=y0;
yin=y0;
for ii=1:n-1
    time=(ii-1)*dt;  
    yout=rk4singlestep(@(t,y)pendulum(t,y,alpha,beta),dt,time,yin);
    Y(:,ii+1)=yout;
    yin=yout;
end
t=[0.1:0.1:20];
figure
hold on
plot(Y(1,:),Y(2,:),'*')

fun=@(t,y)[y(2);sin(y(1))-alpha*y(1)-beta*y(2)]
Y0=[pi/2;1];
tspan=[0:0.1:20];
[ode_t,ode_y]=ode45(fun,tspan,Y0);
plot(ode_y(:,1),ode_y(:,2))
xlabel('Time')
ylabel('Amplitude')
legend('RK4','ode45');
title('Problem 3b(1)')
hold off

%% P3 b(2)
alpha=3;
beta=1;
y0=[pi/2;1];
dt=0.1;
T=20;
n=T/dt;
tspan=[0,T];
Y(:,1)=y0;
yin=y0;
for ii=1:n-1
    time=(ii-1)*dt;  
    yout=rk4singlestep(@(t,y)pendulum(t,y,alpha,beta),dt,time,yin);
    Y(:,ii+1)=yout;
    yin=yout;
end
t=[0.1:0.1:20];
figure
hold on
plot(Y(1,:),Y(2,:),'*')

fun=@(t,y)[y(2);sin(y(1))-alpha*y(1)-beta*y(2)]
Y0=[pi/2;1];
tspan=[0:0.1:20];
[ode_t,ode_y]=ode45(fun,tspan,Y0);
plot(ode_y(:,1),ode_y(:,2))
xlabel('Time')
ylabel('Amplitude')
legend('RK4','ode45');
title('Problem 3b(2)')
hold off

%% P3 c
A=[0 1;-1 0];
P3_c1=eig(A)

A=[0 1;-3 -1];
P3_c2=eig(A)
%%
%Functions
function dy=pendulum(t,y,alpha,beta);
dy=[
    y(2);
    sin(y(1))-alpha*y(1)-beta*y(2);
    ];
end

function yout=rk4singlestep(fun,dt,time,yin);
k1=fun(time,yin);
k2=fun(time+dt/2, yin+(dt/2)*k1);
k3=fun(time+dt/2, yin+(dt/2)*k2);
k4=fun(time+dt, yin+dt*k3);
yout=yin+dt*(k1+2*k2+2*k3+k4)/6;
end
