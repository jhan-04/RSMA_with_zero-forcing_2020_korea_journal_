clc


nt=2;%number of transmitter antenna
P=100;%20dB
i=0;

%for SNR=10:3:30
   % P=10^(SNR/10);

%for gam_dB=-10:1:0
    gam_dB=-5;
    i=i+1;
    j=0;
    % k=0;
 %  for  rho=[0.001:0.01:0.1]
        j=j+1;
        
        
      rho=0;


gam=10^(gam_dB/20);


theta=acos(1-2*rho);
Theta=acos(1-2*rho);


h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
v=0;
t=0;

v=v+1;
h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
h1_;
%  h1_=h1/sqrt(1+norm(h1)^2*P*t);

h2_=h2;
A=[h1_';h2_']*[h1_ h2_];%(10)
a_11=A(1,1);
a_12=A(1,2);
a_22=A(2,2);

mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12))%(8)
f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam)

-1i*angle(a_12);


norm(f_c)^2
abs(h1_'*f_c);
abs(h2_'*f_c);
b1=P*abs(h2_'*f_c)^2;
a1=norm(h1)^2*rho*P;
%a1=norm(h1)^2*P;
R=log2(1+a1*t)+log2(1+b1*(1- t));

Rc1_mul_p=log2(1+P*abs(h1_'*f_c)^2)

Rc2_mul_p=log2(1+P*abs(h2_'*f_c)^2)








Rs_x=R;
P2_x=0;
P1_x=0;
Pc_x=P;
%end





