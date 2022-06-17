clc
clear all
nt=2;%number of transmitter antenna
P=100;%20dB
i=0;
gam_dB=-5;
rho=0.002
%for P=1:1:10
%for gam_dB=-20:1:0
%for rho=[0.1:0.1:1]

%for  rho=[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1]
gam=10^(gam_dB/20);
theta=acos(1-2*rho);
Theta=acos(1-2*rho);


h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1





a=P^2*norm(h1)^2*norm(h2)^2*rho;
b=P*(norm(h1)^2*rho-norm(h2)^2+norm(h1)^2*rho*norm(h2)^2*P);
c=1+norm(h2)^2*P;
if a>0
    t_c1=max(0,min(b/(2*a),1));
else
    t_c1=0;
end

Rs_c1=log2(-a*t_c1^2+b*t_c1+c);
%%가정이 틀릴경우
if log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*t_c1))<=log2(1+norm(h2)^2*(1- t_c1)*P)
    Rs_c1=-10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms k real
%k=1+norm(h1)^2*rho*P*t;

Rk(k)=log2(k)+log2(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
der_Rk(k)=diff(Rk(k),k);
der2_Rkk(k)=diff(der_Rk(k),k);
%der22_Rkk(k)=diff(Rk(k),k,2);
Rk_nolog(k)=(k)*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
der_Rk_nolog(k)=diff(Rk_nolog(k),k);
der2_Rkk_nolog(k)=diff(Rk_nolog(k),k,2);




f(k)=norm(h2)^2*(norm(h1)^2*rho*P+1-k);
g(k)=norm(h1)^2+k*norm(h2)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho);
d_g(k)=diff(g(k),k);
d_f(k)=diff(f(k),k);
dd_g(k)=diff(g(k),k,2);
dd_f(k)=diff(f(k),k,2);
der_Rgf(k)=(1+f/g)+k*(d_f*g-f*d_g)/g^2;
der2_Rgf(k)=2*(d_f/g-(f*d_g)/g^2)+k*(dd_f/g-(d_f*d_g)/g^2)-k*{(d_f*d_g+f*dd_g)/g^2-2*f*d_g^2/g^3};
% der_Rk(k)=diff(Rk(k),k)
% solve(der_Rk(k)==0)
% 
%  T_Rk_nolog(k)=taylor(Rk_nolog,'ExpansionPoint', (2+norm(h1)^2*rho*P)/2, 'Order', 12);
%  plot(c,Rk_nolog(c),c,T_Rk_nolog(c),'*')
%  t(k) = taylor(Rk, 'ExpansionPoint', 1, 'Order',5);%+1/2*taylor(Rk, 'ExpansionPoint', 1, 'Order',2);(2+norm(h1)^2*rho*P)/2
% plot(c,Rk(c),c,t(c),'*')
%  c=[1:norm(h1)^2*rho*P/10:1+norm(h1)^2*rho*P];
%  bbb(k)=(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
% 
%  plot(c,bbb(c))
% plot(c,der2_Rkk_nolog(c),'*',c,der2_Rgf(c))
% plot(c,Rk(c))%,'*',c,der2_Rgf(c))
% plot(c,Rk_nolog(c),'*',T_Rk_nolog,f(c))%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   R_old=-100;
   R_new=-10;
   K=(2+norm(h1)^2*rho*P)/2;
   step=1;
   while abs((R_old-R_new)/R_old)>0.0001
       double(der_Rk(K));
%K=max(min(1+norm(h1)^2*rho*P,K+20*double(der_Rk(K))/sqrt(step)),1);
K=max(min(1+norm(h1)^2*rho*P,K+norm(h1)^2*rho*P*double(der_Rk(K))/sqrt(step)),1);

       step=step+1;
       R_old=R_new;
       R_new=double(Rk(K));

   end
   K;
   tk=(K-1)/(norm(h1)^2*rho*P);
   Rs_c30=R_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_mul_paper(gam_dB,rho,P);
Rs_c31=Rs_x;
Rs_c32=log2(1+norm(h1)^2*rho*P);



[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma_paper(gam_dB,rho,P);


a=[Rs_c1,Rs_c30,Rs_c31,Rs_c32 ,Rs_x]

%end
[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma22(gam_dB,rho,P)