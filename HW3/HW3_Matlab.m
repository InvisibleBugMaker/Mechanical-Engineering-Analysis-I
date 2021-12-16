%Zhaoyi Jiang
%ME564 HW3

%%
%P2 (1)
clc
clear all


m=10; k=90; F0=10; w0=sqrt(k/m);
c=0;
multiper_of_w=0.9;
w=multiper_of_w*w0;

A=[0,1;-k/m,-c/m];
[P,D]=eig(A);
C=[0;1];
tspan=[0,50];
dt=0.1;
time=0:dt:max(tspan);

fun1=@(t)P*[exp(D(1,1)*t),0;0,exp(D(2,2)*t)]/P;
fun2=@(t)cos(w*t);

ymatrix=zeros(2,length(time));
z1=zeros(2,length(time));
z2=zeros(1,length(time));
for i=1:length(time)
    ymatrix(:,i)=fun1(time(i))*C;
    f1m=fun1(time(i));
    z1(:,i)=f1m(:,2);
    z2(i)=fun2(time(i));
end

disptoge=dt*conv(z2,z1(1,:));
disp_s=disptoge(1:length(time));
veltoge=dt*conv(z2,z1(2,:));
vel_s=veltoge(1:length(time));

y=ymatrix(1,:)+disp_s;
dydt=ymatrix(2,:)+vel_s;

[t_ode45,Y_ode45]=ode45(@(t,Y)A*Y+[0;fun2(t)],tspan,C);

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,y,'k.',t_ode45,Y_ode45(:,1),'r');
title('Q1 Displacement');ylabel('y');xlabel('t')
legend('Matrix Exponential','ode45')
hold off

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,dydt,'k.',t_ode45,Y_ode45(:,2),'r')
title('Q1 Velocity');ylabel('dy/dt');xlabel('t')
legend('Matrix Exponential','ode45')
hold off



%%
%P2 (2)
clc
clear all

m=10; k=90; F0=10; w0=sqrt(k/m);
c=0;
multiper_of_w=1;
w=multiper_of_w*w0;

A=[0,1;-k/m,-c/m];
[P,D]=eig(A);
C=[0;1];
tspan=[0,50];
dt=0.1;
time=0:dt:max(tspan);

fun1=@(t)P*[exp(D(1,1)*t),0;0,exp(D(2,2)*t)]/P;
fun2=@(t)cos(w*t);

ymatrix=zeros(2,length(time));
z1=zeros(2,length(time));
z2=zeros(1,length(time));
for i=1:length(time)
    ymatrix(:,i)=fun1(time(i))*C;
    f1m=fun1(time(i));
    z1(:,i)=f1m(:,2);
    z2(i)=fun2(time(i));
end

disptoge=dt*conv(z2,z1(1,:));
disp_s=disptoge(1:length(time));
veltoge=dt*conv(z2,z1(2,:));
vel_s=veltoge(1:length(time));

y=ymatrix(1,:)+disp_s;
dydt=ymatrix(2,:)+vel_s;

[t_ode45,Y_ode45]=ode45(@(t,Y)A*Y+[0;fun2(t)],tspan,C);

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,y,'k.',t_ode45,Y_ode45(:,1),'r');
title('Q2 Displacement');ylabel('y');xlabel('t')
legend('Matrix Exponential','ode45')
hold off

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,dydt,'k.',t_ode45,Y_ode45(:,2),'r')
title('Q2 Velocity');ylabel('dy/dt');xlabel('t')
legend('Matrix Exponential','ode45')
hold off



%%
%P2 (3)
clc
clear all


m=10; k=90; F0=10; w0=sqrt(k/m);
c=10;
multiper_of_w=0.5;
w=multiper_of_w*w0;

A=[0,1;-k/m,-c/m];
[P,D]=eig(A);
C=[0;1];
tspan=[0,50];
dt=0.1;
time=0:dt:max(tspan);

fun1=@(t)P*[exp(D(1,1)*t),0;0,exp(D(2,2)*t)]/P;
fun2=@(t)cos(w*t);

ymatrix=zeros(2,length(time));
z1=zeros(2,length(time));
z2=zeros(1,length(time));
for i=1:length(time)
    ymatrix(:,i)=fun1(time(i))*C;
    f1m=fun1(time(i));
    z1(:,i)=f1m(:,2);
    z2(i)=fun2(time(i));
end

disptoge=dt*conv(z2,z1(1,:));
disp_s=disptoge(1:length(time));
veltoge=dt*conv(z2,z1(2,:));
vel_s=veltoge(1:length(time));

y=ymatrix(1,:)+disp_s;
dydt=ymatrix(2,:)+vel_s;

[t_ode45,Y_ode45]=ode45(@(t,Y)A*Y+[0;fun2(t)],tspan,C);

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,y,'k.',t_ode45,Y_ode45(:,1),'r');
title('Q3 Displacement');ylabel('y');xlabel('t')
legend('Matrix Exponential','ode45')
hold off

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,dydt,'k.',t_ode45,Y_ode45(:,2),'r')
title('Q3 Velocity');ylabel('dy/dt');xlabel('t')
legend('Matrix Exponential','ode45')
hold off



%%
%P2 (4)
clc
clear all

m=10; k=90; F0=10; w0=sqrt(k/m);
c=60;
multiper_of_w=0.5;
w=multiper_of_w*w0;

A=[0,1;-k/m,-c/m];
[P,D]=eig(A);
C=[0;1];
tspan=[0,50];
dt=0.1;
time=0:dt:max(tspan);

fun1=@(t)P*[exp(D(1,1)*t),0;0,exp(D(2,2)*t)]/P;
fun2=@(t)cos(w*t);

ymatrix=zeros(2,length(time));
z1=zeros(2,length(time));
z2=zeros(1,length(time));
for i=1:length(time)
    ymatrix(:,i)=fun1(time(i))*C;
    f1m=fun1(time(i));
    z1(:,i)=f1m(:,2);
    z2(i)=fun2(time(i));
end

disptoge=dt*conv(z2,z1(1,:));
disp_s=disptoge(1:length(time));
veltoge=dt*conv(z2,z1(2,:));
vel_s=veltoge(1:length(time));

y=ymatrix(1,:)+disp_s;
dydt=ymatrix(2,:)+vel_s;

[t_ode45,Y_ode45]=ode45(@(t,Y)A*Y+[0;fun2(t)],tspan,C);

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,y,'k.',t_ode45,Y_ode45(:,1),'r');
title('Q4 Displacement');ylabel('y');xlabel('t')
legend('Matrix Exponential','ode45')
hold off

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,dydt,'k.',t_ode45,Y_ode45(:,2),'r')
title('Q4 Velocity');ylabel('dy/dt');xlabel('t')
legend('Matrix Exponential','ode45')
hold off



%%
%P2 (5)
clc
clear all


m=10; k=90; F0=10; w0=sqrt(k/m);
c=100;
multiper_of_w=0.5;
w=multiper_of_w*w0;

A=[0,1;-k/m,-c/m];
[P,D]=eig(A);
C=[0;1];
tspan=[0,50];
dt=0.1;
time=0:dt:max(tspan);

fun1=@(t)P*[exp(D(1,1)*t),0;0,exp(D(2,2)*t)]/P;
fun2=@(t)cos(w*t);

ymatrix=zeros(2,length(time));
z1=zeros(2,length(time));
z2=zeros(1,length(time));
for i=1:length(time)
    ymatrix(:,i)=fun1(time(i))*C;
    f1m=fun1(time(i));
    z1(:,i)=f1m(:,2);
    z2(i)=fun2(time(i));
end

disptoge=dt*conv(z2,z1(1,:));
disp_s=disptoge(1:length(time));
veltoge=dt*conv(z2,z1(2,:));
vel_s=veltoge(1:length(time));

y=ymatrix(1,:)+disp_s;
dydt=ymatrix(2,:)+vel_s;

[t_ode45,Y_ode45]=ode45(@(t,Y)A*Y+[0;fun2(t)],tspan,C);


figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,y,'k.',t_ode45,Y_ode45(:,1),'r');
title('Q5 Displacement');ylabel('y');xlabel('t')
legend('Matrix Exponential','ode45')
hold off

figure;hold on; grid on;xlim([min(time),max(time)]);
plot(time,dydt,'k.',t_ode45,Y_ode45(:,2),'r')
title('Q5 Velocity');ylabel('dy/dt');xlabel('t')
legend('Matrix Exponential','ode45')
hold off