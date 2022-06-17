nt=2;%number of transmitter antenna
i=0;

SNR=20;
P=10^(SNR/10);

gam_dB=-5;
i=i+1;
j=0;

for rho=0.001:0.1:1

%rho=0.1010;
gam=10^(gam_dB/20);


theta=acos(1-2*rho);
Theta=acos(1-2*rho);


h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-1,Rc1>Rc2



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-3,Rc1==Rc2


U=1000;
L=0;
T=(U+L)/2;

k_region=[1:(norm(h1)^2*rho*P)/30:1+norm(h1)^2*rho*P];

syms k real
QQ(k)=k+k*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));


while abs((U-L)/U)>0.001
    
    T=(U+L)/2;
    cond1=QQ(k)>=T;
    cond2=1<=k;
    cond3=k<=1+norm(h1)^2*rho*P;
    conds=[cond1 cond2 cond3];
    sol = solve(conds, k, 'Real', true,'ReturnConditions', true, 'IgnoreAnalyticConstraints', true);
    %     sol.k;
    %     sol.parameters;
    %     sol.conditions;
    
    condWithValues = subs(sol.conditions, sol.k,k_region );
    feasible=isAlways(condWithValues);
    
    if isempty(feasible)||sum(feasible==1)==0
        U=T;
    else
        L=T;
        d=find(feasible==1);
        k_min=k_region(min(d));
        k_max=k_region(max(d));
        k_region=[k_min:(k_max-k_min)/30:k_max];
        if k_min==k_max
            K=k_max;
        end
    end
    
    if k_min==k_max
        break;
    end
    K= mean(k_region);
    
end


t_c30=(K-1)/(norm(h1)^2*rho*P);
Rs_c30=log2(double(QQ(K)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a=max(Rs_2_1,Rs_c30);
if a==Rs_2_1
    tou_2=t_2_1;
    Rs_2=Rs_2_1;
end
if a==Rs_c30
    tou_2=t_c30;
    Rs_2=Rs_c30;
end



P2_x=0;
P1_x=P*tou_2;
Pc_x=(1-tou_2)*P;



[tou_y, P1_y,P2_y, Pc_y,Rs_y]=RS_noma22(P,h1,h2);
Rs_2
Rs_y
end
