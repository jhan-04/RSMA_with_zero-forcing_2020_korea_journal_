
clc
clear all
P=100;
gam_dB=-2;
rho=0.5;
%%%%%%%%%%%%%channel
gam=10^(gam_dB/20);



theta=acos(1-2*rho);


h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];

%%%%%%%%%%%%%%%%%%%%%%%
f_c=[1;1]/sqrt(2);


B=h2*h2';

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         cvx_begin quiet
%         variable t(1)
%         variable x(1)
%         variable F(2,2)
%         minimize(-x)
%         subject to
%         (norm(h1)^2*rho*P-trace(B*F))*t+trace(B*F)>=x;
%         (1+norm(h1)^2*rho-trace(B*F)+norm(h1)^2*rho*P*trace(B*F))*P*t-(norm(h1)^2*rho*P*P*trace(B*F))*t^2>=x;
%
%         cvx_end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
t=0;
f_c_old=[-100;-100];
t_old=0;
while abs(abs(f_c-f_c_old).^2+(t-t_old))>0.001
    
    
    
    
    a=norm(h1)^2*rho*P*t;
    b=(1-t)*P;
    
    f_c_old=f_c;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    variable F(2,2) hermitian
    variable x(1)
    minimize(-x)
    subject to
    a+b*trace(B*F)>=x;
    1+a+(b+b*a)*trace(B*F)>=x;
    trace(F)==1;
    F== semidefinite(2);
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    F
    
    [U,W,Z] = svds(F);%V=U*W*Z'
    obj_max=-100;
    n=length(W);
    for m=1:500000
        r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        v_=U*sqrt(W)*r;
        v_=v_/sqrt(v_'*v_);
        
        obj=min(a+b*v_'*B*v_,1+a+(b+b*a)*v_'*B*v_);
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
    end
    f_c=v_max;
    norm(f_c)^2;
    
    
    
    
    
    
    a=norm(h1)^2*rho*P-abs(h2'*f_c)^2*P;
    b=abs(h2'*f_c)^2;
    c=(norm(h1)^2*rho-abs(h2'*f_c)^2+norm(h1)^2*rho*abs(h2'*f_c)^2*P)*P;
    d=-norm(h1)^2*rho*abs(h2'*f_c)^2*P*P;
    
    t_old=t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    variable t(1)
    variable x(1)
    minimize(-x)
    subject to
    a*t+b>=x;
    1+c*t+d*t^2>=x;
    t<=1;
    t>=0;
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t
    
    
    
    
    
end
h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
% h1_=h1/sqrt(1+norm(h1)^2*P*t);
h2_=h2;


Rc=min(log2(1+abs(h1_'*f_c)^2*(1- t)*P),log2(1+abs(h2_'*f_c)^2*(1- t)*P));
R=  Rc+ log2(1+norm(h1)^2*P*t);

[MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_noma(gam_dB,rho,P);
