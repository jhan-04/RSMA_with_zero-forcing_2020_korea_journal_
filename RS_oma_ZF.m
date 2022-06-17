function [MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_oma_ZF(P,h1,h2)


 rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;

MA_p=3;


        h1_=h1/norm(h1);
        h2_=h2/norm(h2);
        %%%%%%%%%%%%%%%%%%%%%%%

P2_p=0;
Pc_p=0;
P1_p=P;

tou_p=1;


Rs_p=log2(1+norm(h1)^2*P*rho);

end
    