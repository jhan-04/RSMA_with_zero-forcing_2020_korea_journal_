function [tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma22(P,h1,h2)

 rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-1,Rc1>Rc2

% a=P^2*norm(h1)^2*norm(h2)^2*rho;
% b=P*(norm(h1)^2*rho-norm(h2)^2+norm(h1)^2*rho*norm(h2)^2*P);
% c=1+norm(h2)^2*P;
% if a>0
%     t_2_1=max(0,min(b/(2*a),1));
% else
%     t_2_1=0;
% end
% 
% Rs_2_1=log2(-a*t_2_1^2+b*t_2_1+c);
% %%가정이 틀릴경우
% if log2(1+norm(h1)^2*(1-rho)*(1-t_2_1)*P/(1+P*norm(h1)^2*rho*t_2_1))<=log2(1+norm(h2)^2*(1- t_2_1)*P)
%     Rs_2_1=-10;
% end




e=-P^2*norm(h1)^2*norm(h2)^2*rho;
f=P*(norm(h1)^2*rho-norm(h2)^2+norm(h1)^2*rho*norm(h2)^2*P);
g=1+norm(h2)^2*P;
if e<0
    t_2_1=max(0,min(-f/(2*e),1));
else
    t_2_1=0;
end

Rs_2_1=log2(e*t_2_1^2+f*t_2_1+g);
%%가정이 틀릴경우
if log2(1+norm(h1)^2*(1-rho)*(1-t_2_1)*P/(1+P*norm(h1)^2*rho*t_2_1))<=log2(1+norm(h2)^2*(1- t_2_1)*P)
    Rs_2_1=-10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms k real
%k=1+norm(h1)^2*rho*P*t;

Rk(k)=log2(k*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho))));
der_Rk(k)=diff(Rk(k),k);
der2_Rkk(k)=diff(der_Rk(k),k);

%plot([1:norm(h1)^2*rho*P/10:1+norm(h1)^2*rho*P],der2_Rkk([1:norm(h1)^2*rho*P/10:1+norm(h1)^2*rho*P]))
% der_Rk(k)=diff(Rk(k),k)
% solve(der_Rk(k)==0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QQ(k)=k*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));

   R_old=-100;
   R_new=-10;
   K=(2+norm(h1)^2*rho*P)/2;
   step=1;
   while abs((R_old-R_new)/R_old)>0.0001
       double(der_Rk(K));
K=max(min(1+norm(h1)^2*rho*P,K+norm(h1)^2*rho*P*double(der_Rk(K))/sqrt(step)),1);
       step=step+1;
       R_old=R_new;
       R_new=double(Rk(K));

   end
   K;
   t_c30=(K-1)/(norm(h1)^2*rho*P);
   Rs_c30=R_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_mul_paper(gam_dB,rho,P);
% Rs_c31=Rs_x;
% t_c31=0;
% Rs_c32=log2(1+norm(h1)^2*rho*P);
% t_c32=1;
% 
% 
% m=max(Rs_2_1,Rs_c30);
% g=max(Rs_c31,Rs_c32);
% a=max(m,g);
a=max(Rs_2_1,Rs_c30);
if a==Rs_2_1
     tou_x=t_2_1;
     Rs_x=Rs_2_1;
end
if a==Rs_c30
    tou_x=t_c30;
    Rs_x=Rs_c30;
 end
% if a==Rs_c31
%     tou_x=t_c31;
%     Rs_x=Rs_c31;
% end
% if a==Rs_c32
%     tou_x=t_c32;
%     Rs_x=Rs_c32;
% end


    P2_x=0;
    P1_x=P*tou_x;
    Pc_x=(1-tou_x)*P;




end
