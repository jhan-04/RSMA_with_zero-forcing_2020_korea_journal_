clc
nt=2;%number of transmitter antenna
P=100;%20dB
i=0;
gam_dB=-5;


for rho=[0.01:0.05:0.5]
    
    rho=0.01;
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
    
    Rs_c1=log2(-a*t_c1^2+b*t_c1+c)
    log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*Rs_c1))-log2(1+norm(h2)^2*(1- t_c1)*P)
    %%가정이 틀릴경우
    if log2(1+norm(h1)^2*(1-rho)*(1-t_c1)*P/(1+P*norm(h1)^2*rho*Rs_c1))<=log2(1+norm(h2)^2*(1- t_c1)*P)
        Rs_c1=-1000;
    end
    
    %%%%%%%%%%%%%%%%%%%3
    
    
    
    a=rho*norm(h1)^2;
    A=(1/a+P)*h2*h2'-(3/a)*h1*h1';
    b=[0;0];
    X = linsolve(A,b);
    f_c=X/norm(X);
    %         syms x [2,1]
    %         equ1=x'*A==[0 0];
    %         equ2=x'*x==1;
    %         equ=[equ1,equ2];
    %         S=vpasolve(equ,x);
    %         if isempty(S.x1)==0
    %         f_c(1,1)=eval(S.x1);
    %        f_c(2,1)=eval(S.x2);
    %
    n=abs(h1'*f_c)^2;
    m=abs(h2'*f_c)^2;
    Rs_c30=log2(n*(1/m+1/(norm(h1)^2*rho)+P-n/(m*norm(h1)^2*rho)))
    %     end
    
    A=h2*h2'-h1*h1';
%     syms x [2,1]
%     equ1=x'*A*x==0;
%     equ2=x'*x==1;
%     equ=[equ1,equ2];
%     S=vpasolve(equ,x);
%     f_c(1,1)=eval(S.x1);
%     f_c(2,1)=eval(S.x2);

    abs(h2'*f_c)^2-abs(h1'*f_c)^2
    Rs_c31= log2(1+abs(h2'*f_c)^2*P)
    
    Rs_c32=log2(1+norm(h1)^2*rho*P)
    
    
%      [MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_noma(gam_dB,rho,P);
%     Rs_p
%     tou_p
%     gam_dB
    [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x]=RS_noma_paper(gam_dB,rho,P);
    Rs_x
    tou_x
end


%end
%end