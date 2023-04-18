% Call ASD_H.m to validate output from main.c

t=(0:0.01:100);
X=1024*cos(2*pi*t*1.0+1.234) + 256*cos(2*pi*t*2.0+0.123) + 64*sin(2*pi*t*4.0) + 16*cos(2*pi*t*8.0+5.67) + 4* cos(2*pi*t*16.0-5.67);  

n2u=importdata ("n2u_KAISER.test");
u2n=importdata("u2n_KAISER.test");
n2n=importdata("n2n_KAISER.test");
flt=importdata("filter_KAISER.test");
fltn2u=importdata("filter_KAISER_n2u.test");
figure();hold on;
plot(u2n(:,1),u2n(:,2),"b-","LineWidth",2);
plot(n2u(:,1),n2u(:,2),"r--");
plot(n2n(:,1),n2n(:,2),"m--");
plot(flt(:,1),flt(:,2),"c--o");
plot(fltn2u(:,1),fltn2u(:,2),"y--.");
ASD_H(X,0.01,"ks","k.");
title("KAISER");

figure();hold on;
n2u=importdata ("n2u_HANN.test");
u2n=importdata("u2n_HANN.test");
n2n=importdata("n2n_HANN.test");
u2u=importdata("u2u_HANN.test");
flt=importdata("filter_HANN.test");
fltn2u=csvread("n2u_HANN.test_flt.csv");
plot(u2n(:,1),u2n(:,2),"b-","LineWidth",2);
plot(u2u(:,1),u2u(:,2),"g-");
plot(n2u(:,1),n2u(:,2),"r--");
plot(n2n(:,1),n2n(:,2),"m--");
plot(flt(:,1),flt(:,2),"c--o");
plot(fltn2u(2:end,1),fltn2u(2:end,2),"y--.");
ASD_H(X,0.01,"hn","k.");hold on;
title("HANNING");