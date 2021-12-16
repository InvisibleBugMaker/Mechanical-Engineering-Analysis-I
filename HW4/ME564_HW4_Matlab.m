%Zhaoyi Jiang
%ME564 HW4
clc
clear all
close all
%Problem 2  Part C

fun=@(t,y)[y(2);sin(y(1))-2*y(1)-y(2)]
Y0=[0;1];
tspan=[0 20];
[T,Y]=ode45(fun,tspan,Y0);
plot(T,Y(:,1),T,Y(:,2))
title('Problem 2  Part C Solution')
xlabel('Time')
ylabel('Amplitude')
legend('Y1(Theta)','Y2(Omega)')

figure
plot(Y(:,1),Y(:,2))
title('Problem 2  Part C Phase Plane')
xlabel('Y1')
ylabel('Y2')
