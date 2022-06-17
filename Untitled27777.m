
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
%%가정이 틀릴경우
if log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*t_c1))<=log2(1+norm(h2)^2*(1- t_c1)*P)
    Rs_c1=-10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms k real
%k=1+norm(h1)^2*rho*P*t;

Rk(k)=log2(k*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho))));
der_Rk(k)=diff(Rk(k),k);
der2_Rkk(k)=diff(der_Rk(k),k);

plot([1:norm(h1)^2*rho*P/10:1+norm(h1)^2*rho*P],der2_Rkk([1:norm(h1)^2*rho*P/10:1+norm(h1)^2*rho*P]))
% der_Rk(k)=diff(Rk(k),k)
% solve(der_Rk(k)==0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   R_old=-100;
   R_new=-10;
   K=(2+norm(h1)^2*rho*P)/2;
   step=1;
   while abs((R_old-R_new)/R_old)>0.001
       double(der_Rk(K));
K=max(min(1+norm(h1)^2*rho*P,K+20*double(der_Rk(K))/sqrt(step)),1);
       step=step+1;
       R_old=R_new;
       R_new=double(Rk(K));

   end
   K;
   t_c30=(K-1)/(norm(h1)^2*rho*P);
   Rs_c30=R_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_mul_paper(gam_dB,rho,P);
Rs_c31=Rs_x;
t_c31=0;
Rs_c32=log2(1+norm(h1)^2*rho*P);
t_c32=1;


a=max(Rs_c1,Rs_c30,Rs_c31,Rs_c32);
if a==Rs_c1
     tou_x=t_c1;
end
if a==Rs_c30
    tou_x=t_c30;
end
if a==Rs_c31
    tou_x=t_c31;
end
if a==Rs_c32
    tou_x=t_c32;
end




