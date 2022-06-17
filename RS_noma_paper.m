%Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
function [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma_paper(P,h1,h2)
MA_x=2;
P2_x=0;


rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;
Nt=length(h1);

v=0;
for t=0:1/200:1
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
    
    
    abs(h1_'*f_c);
    abs(h2_'*f_c);
    b1=P*abs(h2'*f_c)^2;
    a1=norm(h1)^2*rho*P;
    %a1=norm(h1)^2*P;
    R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
    
    


    
    
    
end

R;
t=0:1/200:1;
Rs_noma_mul=max(R);%0<=t<1
k=find(Rs_noma_mul==R);
t1=t(k);
%
%             Rs_oma=log2(1+norm(h1)^2*rho*P);

tou_x=t1;
Rs_x=Rs_noma_mul;
P2_x=0;
P1_x=P*t1;
Pc_x=(1-t1)*P;
%end




end



