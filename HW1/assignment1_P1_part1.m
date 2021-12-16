clc
clear all
close all
A=2;    B=2;
figure
y0_array=[0.1 0.5];
C=1./y0_array-B/A;
for c=C
sol=@(t) 1./(c.*exp(-A.*t)+B/A);
fplot(sol);
end
f=@(t,y) A*y-B*y.^2;
tspan=[0,5];
[t_nu,y_nu]=ode45(f,tspan,y0_array);
p2=plot(t_nu,y_nu,'b');


grid on
hold on
legend([p2],{'ode45'});
xlim([0,5]);
ylim([0,1.2]);
title('Analytical and Numerical Solutions of Logistic Equation A=2,B=2');
ylabel('y');xlabel('t');