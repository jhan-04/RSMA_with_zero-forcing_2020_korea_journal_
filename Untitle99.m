clc
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
    
    %%%%%%%%%%%%%%%%%%%3
    
    syms F [2,2] 
    syms f [2,1] 
 A=h1*h1';
 B=h2*h2';
   % equ1=gradient(abs(h1'*f)^2*(1/abs(h2'*f)^2+1/(norm(h1)^2*rho)+P-abs(h1'*f)^2/(abs(h2'*f)^2*norm(h1)^2*rho)),f)==0.1;
   %t=(abs(h1'*f)^2-abs(h2'*f)^2)/(norm(h1)^2*rho*P*abs(h2'*f)^2);
%    equ1=diff(log2(1+norm(h1)^2*rho*P*(abs(h1'*f)^2-abs(h2'*f)^2)/(norm(h1)^2*rho*P*abs(h2'*f)^2)+abs(h1'*f)^2*(1-(abs(h1'*f)^2-abs(h2'*f)^2)/(norm(h1)^2*rho*P*abs(h2'*f)^2))*P),f1)==0;
%    equ2=diff(log2(1+norm(h1)^2*rho*P*(abs(h1'*f)^2-abs(h2'*f)^2)/(norm(h1)^2*rho*P*abs(h2'*f)^2)+abs(h1'*f)^2*(1-(abs(h1'*f)^2-abs(h2'*f)^2)/(norm(h1)^2*rho*P*abs(h2'*f)^2))*P),f2)==0;
%    equ3=f'*f==1;
%    equ=[equ1,equ2,equ3];
%   S=solve(equ,f1,f2);
%    f_c=[eval(S.f1); eval(S.f2)]

  

% equ1=diff(trace(A*F)/trace(B*F)+trace(A*F)*(P+1/(norm(h1)^2*rho))-trace(A*F)^2/(norm(h1)^2*rho*trace(B*F)),F1_1)==0;
% equ2=diff(trace(A*F)/trace(B*F)+trace(A*F)*(P+1/(norm(h1)^2*rho))-trace(A*F)^2/(norm(h1)^2*rho*trace(B*F)),F1_2)==0;
% equ3=diff(trace(A*F)/trace(B*F)+trace(A*F)*(P+1/(norm(h1)^2*rho))-trace(A*F)^2/(norm(h1)^2*rho*trace(B*F)),F2_1)==0;
% equ4=diff(trace(A*F)/trace(B*F)+trace(A*F)*(P+1/(norm(h1)^2*rho))-trace(A*F)^2/(norm(h1)^2*rho*trace(B*F)),F2_2)==0;
% equ5=trace(F)==1;
% equ6=F==F';
% 
% equ=[equ1,equ2,equ3,equ4,equ5];
% S=solve(equ,F1_1,F1_2,F2_1,F2_2);
% F_c=round([eval(S.F1_1(1)),eval(S.F1_2(1)); eval(S.F2_1(1)),eval(S.F2_2(1))],4)
% trace(A*F_c)
% trace(B*F_c)
% t=(trace(A*F_c)-trace(B*F_c))/(norm(h1)^2*rho*P*trace(B*F_c))
% 
% equ1=A*trace(B*F)-B*trace(A*F)+A.*(P+1/(norm(h1)^2*rho))*trace(B*F)^2+B*1/(norm(h1)^2*rho)*trace(A*F)^2-2*A*trace(A*F)*trace(B*F)/(norm(h1)^2*rho)==0;
% equ5=trace(F)==1;
% equ6=F==F';%
% trace(A*F_c)
% trace(B*F_c)
% t=(trace(A*F_c)-trace(B*F_c))/(norm(h1)^2*rho*P*trace(B*F_c))
% 
% equ=[reshape(equ1,1,[]),equ5];
% S=solve(equ,F1_1,F1_2,F2_1,F2_2);
% F_c=round([eval(S.F1_1(1)),eval(S.F1_2(1)); eval(S.F2_1(1)),eval(S.F2_2(1))],4)

%     F_c=round([eval(S.F1_1(1)),eval(S.F1_2(1)); eval(S.F2_1(1)),eval(S.F2_2(1))],4)
%     trace(A*F_c)
%     trace(B*F_c)
%     t=(trace(A*F_c)-trace(B*F_c))/(norm(h1)^2*rho*P*trace(B*F_c))
%    if length(eval(S.F1_1))~=0
% 
%     F_c(1,1)=eval(S.F1_1);
%     F_c(2,1)=eval(S.F1_1);
%     n=abs(h1'*f_c)^2;
%     m=abs(h2'*f_c)^2;
%     Rs_c30=log2(n*(1/m+1/(norm(h1)^2*rho)+P-n/(m*norm(h1)^2*rho)))
%    end
    
   [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_mul_paper(gam_dB,rho,P);
    Rs_c31=Rs_x;
    Rs_c32=log2(1+norm(h1)^2*rho*P);
    
    

    [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma_paper(gam_dB,rho,P);

    
    a=[Rs_c1,Rs_c31,Rs_c32 ,Rs_x]
end


%end
%end