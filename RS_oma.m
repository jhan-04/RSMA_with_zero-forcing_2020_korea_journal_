function [MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_TDMA(P,h1,h2)
MA_p=3;


P2_p=0;
Pc_p=0;
P1_p=P;

tou_p=1;


Rs_p=1/2*log2(1+norm(h1)^2*P)+1/2*log2(1+norm(h2)^2*P);

end
    