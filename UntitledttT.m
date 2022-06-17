%Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
clear all
rho=0.3;
gam_dB=-5;
P=100;
MA_x=2;
P2_x=0;


gam=10^(gam_dB/20);


theta=acos(1-2*rho);
Theta=acos(1-2*rho);


h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
v=0;
for t=0:1/1000:1
    v=v+1;
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
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
    
    
    
    b1=P*abs(h2'*f_c)^2;
    a1=norm(h1)^2*rho*P;
    %a1=norm(h1)^2*P;
    t
    abs(h1_'*f_c)^2
    abs(h2_'*f_c)^2
    (a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12))
   
    R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
 k=1+norm(h1)^2*rho*P*t;
   a=norm(h1)^2*norm(h2)^2*rho;
   Rk(v)=log2(k)+log2(1+(1-t)*P*a/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
% a/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho))
    
    
    
end


R;
t=0:1/1000:1;
plot(t,Rk)

Rs_noma_mul=max(R);%0<=t<1
%k1=find(Rs_noma_mul==R);
%t1=t(k1);
%
%             Rs_oma=log2(1+norm(h1)^2*rho*P);

tou_x=t1;
Rs_x=Rs_noma_mul;
P2_x=0;
P1_x=P*t1;
Pc_x=(1-t1)*P;
%end








