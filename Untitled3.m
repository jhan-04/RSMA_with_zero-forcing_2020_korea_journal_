%when imperfect CSI situation, we don't know channel error is occured.
%function [t_x,P1_x,P2_x,Pc_x,Rs_x]=RS_sdr_perf_2(P,h1,h2,beta)



rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;
Nt=length(h1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(1)P1>0,P2>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_=h1/norm(h1);
h2_=h2/norm(h2);
A=h1_*h1_';
B=h2_*h2_';
%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
variable F(Nt,Nt) hermitian
variable x(1)
minimize(-x)
subject to
trace(A*F)>=x;
trace(B*F)>=x;
trace(F)==1;
F== semidefinite(Nt);
cvx_end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
F;




[U,W,Z] = svds(F);%V=U*W*Z'
obj_max=-100;
n=length(W);
for m=1:100000
    r=sqrt(1/2)*(randn(Nt,1)+1i*randn(Nt,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    v_=U*sqrt(W)*r;
    v_=v_/sqrt(v_'*v_);
    obj=min(v_'*A*v_,v_'*B*v_);
    if obj>=obj_max
        obj_max=obj;
        v_max=v_;
    end
end
fc_1=v_max;

if (fc_1'*A*fc_1)>=(fc_1'*B*fc_1)
    ha=h1;
    hb=h2;
    Gam=(1/rho)*(1/norm(hb)^2-1/norm(ha)^2);
    b=norm(ha)^2*rho*P/2;
    a=1+(Gam/P)*b;
    d=norm(hb)^2*rho*P/2-abs(hb'*fc_1)^2*P;
    c=1-(Gam/P)*d+abs(hb'*fc_1)^2*(P-Gam);
    
else
    ha=h2;
    hb=h1;
    Gam=(1/rho)*(1/norm(hb)^2-1/norm(ha)^2);
    b=norm(ha)^2*rho*P/2;
    a=1+(Gam/P)*b;
    d=norm(hb)^2*rho*P/2-abs(hb'*fc_1)^2*P;
    c=1-(Gam/P)*d+abs(hb'*fc_1)^2*(P-Gam);
end

if (b*d)<0
    t1=min(-a/(2*b)-c/(2*d),1);
    Rs_1=log2(a*c+(a*d+b*c)*t1+b*d*t1^2);
else
    t1=1;
    Rs_121=log2(a*c+(a*d+b*c)*t1+b*d*t1^2);
    t1=0;
    Rs_120=log2(a*c+(a*d+b*c)*t1+b*d*t1^2);
    if  Rs_121>Rs_120
        t1=1;
        Rs_1=log2(a*c+(a*d+b*c)*t1+b*d*t1^2);
    else
        t1=0;
        Rs_1=log2(a*c+(a*d+b*c)*t1+b*d*t1^2);
    end
end


if (t1*P >= (1/rho)*(1/norm(h2)^2-1/norm(h1)^2)) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.

    
else
    Rs_1=-100;
end



    P1_1=t1*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    P2_1=t1*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    Pc_1=(1-t1)*P;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if (t1*P >= (1/rho)*(1/norm(h2)^2-1/norm(h1)^2)) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.

else
    Rs_1=-100;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% if P1>=0,P2=0
v=0;

for t=0:1/200:1
    v=v+1;
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
    % h1_=h1/sqrt(1+norm(h1)^2*P*t);
    h2_=h2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=h1_*h1_';
    B=h2_*h2_';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    variable F(2,2) hermitian
    variable x(1)
    minimize(-x)
    subject to
    trace(A*F)>=x;
    trace(B*F)>=x;
    trace(F)==1;
    F== semidefinite(2);
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    trace(A*F);
    trace(B*F);
    F;
    
    [U,W,Z] = svds(F);%V=U*W*Z'
    obj_max=-100;
    n=length(W);
    for m=1:10000
        r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        v_=U*sqrt(W)*r;
        v_=v_/sqrt(v_'*v_);
        obj=min(v_'*A*v_,v_'*B*v_);
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
    end
    fc_t(:,v)=v_max;
 
    Rc=min(log2(1+abs(h1_'*fc_t(:,v))^2*(1- t)*P),log2(1+abs(h2_'*fc_t(:,v))^2*(1- t)*P));
    R_t(v)=  Rc+ log2(1+norm(h1)^2*rho*P*t);

end

t=0:1/200:1;
Rs_2=max(R_t);%0<=t<1
k=find(Rs_2==R_t);
t2=t(k);
fc_2=fc_t(:,k);

P1_2=P*t2;
P2_2=0;
Pc_2=(1-t2)*P;


maax= max(Rs_1,Rs_2);
if maax==Rs_2
    t_x=t2;
    
    P2_x=P2_2;
    P1_x=P1_2;
    Pc_x=Pc_2;
    fc_x=fc_2;
    
else
    
    t_x=t1;
    P1_x=P1_1;
    P2_x=P2_1;
    Pc_x=Pc_1;
    fc_x=fc_1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculating real sum-rate

Rc_x=min(log2(1+abs(h1'*fc_x)^2*Pc_x/(1+beta*P+norm(h1)^2*rho*P1_x)),log2(1+abs(h2'*fc_x)^2*Pc_x/(1+beta*P+norm(h2)^2*rho*P2_x)));
R1_x=log2(1+beta*P+norm(h1)^2*rho*P1_x)-log2(1+beta*P);
R2_x=log2(1+beta*P+norm(h2)^2*rho*P2_x)-log2(1+beta*P);
Rs_x=Rc_x+R1_x+R2_x;

% Rs_x
% maax
% Rs_1
% Rs_2
end



