%Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
function [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_paper(P,h1,h2)


 rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;

h1_=h1/norm(h1);
h2_=h2/norm(h2);


h12_=[h1_ h2_];
%%%%caculation of p_c, f_c from (7),(8),(9),(10), with (15)
A=[h1_';h2_']*[h1_ h2_];%(10)
a_11=A(1,1);
a_12=A(1,2);
a_22=A(2,2);



mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);%*exp(1i*pi);

norm(f_c)^2;


%%%%%%caculation of

%%%%%%%P1 ~=0,P2 ~=0

Gam=(1/rho)*(1/norm(h2)^2-1/norm(h1)^2);
b=norm(h1)^2*rho*P/2;
a=1+(Gam/P)*b;
d=norm(h2)^2*rho*P/2-abs(h2'*f_c)^2*P;
c=1-(Gam/P)*d+abs(h2'*f_c)^2*(P-Gam);




t12=min(-a/(2*b)-c/(2*d),1);

Rs_12=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
a*d+b*c;

b*d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (t12*P >= (1/rho)*(1/norm(h2)^2-1/norm(h1)^2)) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.
    Rs_1=Rs_12;
    tou_1=t12;
    P1_1=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    P2_1=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    Pc_1=(1-t12)*P;
    %         Rs_x=Rs_12;
    %     tou_x=t12;
    %     P1_x=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    %     P2_x=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    %     Pc_x=(1-t12)*P;
else
    Rs_1=-100;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=0;
for t=0:1/100:1
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
    b1=P*abs(h2_'*f_c)^2;
    a1=norm(h1)^2*rho*P;
    %a1=norm(h1)^2*P;
    R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
     Rc1(v)=log2(1+abs(h1_'*f_c)^2*(1- t)*P);
    Rc2(v)=log2(1+abs(h2_'*f_c)^2*(1- t)*P);
    abs(h1_'*f_c)^2;
    abs(h2_'*f_c)^2;
    
end


t=0:1/100:1;
Rs_noma_mul=max(R);%0<=t<1
k=find(Rs_noma_mul==R);
t1=t(k);

RC1_p=Rc1(k);
RC2_p=Rc2(k);



%
%     tou_x=t1;
%     Rs_x=Rs_noma_mul;
%     P2_x=0;
%     P1_x=P*t1;
%     Pc_x=(1-t1)*P;
%end



maax= max(Rs_1,Rs_noma_mul);
if maax==Rs_noma_mul
    tou_x=t1;
    Rs_x=Rs_noma_mul;
    P2_x=0;
    P1_x=P*t1;
    Pc_x=(1-t1)*P;
    
else
    Rs_x=Rs_12;
    tou_x=t12;
    P1_x=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    P2_x=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    Pc_x=(1-t12)*P;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if tou_x==1
    MA_x=1;%SDMA
end
if (P1_x~=0)&&(P2_x==0)&&(Pc_x~=0)
    MA_x=2;%NOMA
end
if (P1_x~=0)&&(P2_x==0)&&(Pc_x==0)
    MA_x=3;%OMA
end
if (P2_x==0)&&(P1_x==0)%tou(i,j)==0
    MA_x=4;%multicasting
end
if (P1_x~=0)&&(P2_x~=0)&&(Pc_x~=0)
    MA_x=5;%RS
end




end



