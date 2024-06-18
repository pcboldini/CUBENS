clear all
close all

y_max=20;
y_i=3;
s=2;
eta=linspace(0,1,100);

a=(y_max-y_i)*s;
b=1+(a)./(y_max-y_i);

%a = y_max * y_i / (y_max - 2.0 * y_i);
%b = 1+ 2.0 * a / y_max;  

y_1=y_i*(1-cos(pi*eta))./2;
dy_1deta=y_i*pi*(sin(pi*eta))./2;
y_2=y_i+a.*(eta)./(b-eta);
dy_2deta=a*b./(b-eta).^2;

index=(1:1:199);

y=[y_1';y_2(2:end)'];
dy=[dy_1deta';dy_2deta(2:end)'];