clc
clear all
nt=2;%number of transmitter antenna
P=100;%20dB
i=0;
gam_dB=-5;
%for P=1:1:10
%for gam_dB=-20:1:0
%for rho=[0.1:0.1:1]

rho=0.3;
gam=10^(gam_dB/20);
theta=acos(1-2*rho);
Theta=acos(1-2*rho);


rho
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
log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*t_c1))-log2(1+norm(h2)^2*(1- t_c1)*P);
%%가정이 틀릴경우
if log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*t_c1))<=log2(1+norm(h2)^2*(1- t_c1)*P)
    Rs_c1=-10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%3
% syms t real

% h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
% %  h1_=h1/sqrt(1+norm(h1)^2*P*t);
% 
% h2_=h2;
% A=[h1_';h2_']*[h1_ h2_];%(10)
% a_11=A(1,1);
% a_12=A(1,2);
% a_22=A(2,2);
% 
% mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
% mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
% lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
% f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
% 
% 
% 
% b1=P*abs(h2'*f_c)^2;
% a1=norm(h1)^2*rho*P;
% %a1=norm(h1)^2*P;
% R(t)=log2((1+a1*t)*(1+b1*(1- t)));
% plot([0:0.01:1],R([0:0.01:1]))
% der_R(t)=diff(R(t),t);
% der_sum(t)=real((der_R(t)))+imag((der_R(t)));

%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    R_old=-100;
%    R_new=-10;
%    T=0.5;
%    step=1;
%    while abs((R_old-R_new))>0.001
% T=max(min(1,T+0.1*real(double(der_R(T)))/step),0)
%        step=step+1;
%        R_old=R_new
%        R_new=double(R(T))
%
%    end
%    T
%    v=0;
%    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms k real
%k=1+norm(h1)^2*rho*P*t;
%a=norm(h1)^2*norm(h2)^2*rho;
Rk(k)=log2(k)+log2(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
der_Rk(k)=diff(Rk(k),k);
% der_Rk(k)=diff(Rk(k),k)
% solve(der_Rk(k)==0)
%%%%%%%%%%


   R_old=-100;
   R_new=-10;
   K=(2+norm(h1)^2*rho*P)/2;
   step=1;
   while abs((R_old-R_new)/R_old)>0.0001
       double(der_Rk(K))
K=max(min(30,K+20*double(der_Rk(K))/sqrt(step)),1)
       step=step+1;
       R_old=R_new
       R_new=double(Rk(K))

   end
   K
   tk=(K-1)/(norm(h1)^2*rho*P)
   R_new
   v=0;
%    %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt=0:1/100:1
    v=v+1;
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*tt);
    %  h1_=h1/sqrt(1+norm(h1)^2*P*t);
    
    h2_=h2;
    A=[h1_';h2_']*[h1_ h2_];%(10)
    a_11=A(1,1);
    a_12=A(1,2);
    a_22=A(2,2);
    
    mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
    mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
    lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
    f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
    
    
    abs(h1_'*f_c);
    abs(h2_'*f_c);
    b1=P*abs(h2'*f_c)^2;
    a1=norm(h1)^2*rho*P;
    %a1=norm(h1)^2*P;
    RR(v)=log2(1+a1*tt)+log2(1+b1*(1- tt));
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    
end
tt=0:1/100:1;
plot(tt,RR)
Rs_noma_mul=max(RR);%0<=t<1
k=find(Rs_noma_mul==RR);
t1=tt(k);


[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_mul_paper(gam_dB,rho,P);
Rs_c31=Rs_x;
Rs_c32=log2(1+norm(h1)^2*rho*P);



[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma_paper(gam_dB,rho,P);


a=[Rs_c1,Rs_c31,Rs_c32 ,Rs_x]
