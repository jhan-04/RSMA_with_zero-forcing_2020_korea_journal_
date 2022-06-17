%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%%% Fig 3 (a)
clc
clear all

nt=2;%number of transmitter antenna
P=100;%20dB
i=0;

%for SNR=1:2:20
SNR=10;
P=10^(SNR/10);

%for gam_dB=-10:1:0
gam_dB=-5;
i=i+1;
j=0;
% k=0;
for  rho=[0.001,0.020,0.05,0.1,0.2:0.2:1]%[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1]
    j=j+1;
    rho
    gam=10^(gam_dB/20);
    theta=acos(1-2*rho);
    Theta=acos(1-2*rho);
    
    h1=1/sqrt(2)*[1;1];
    h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
    [MA_p(j),tou_p(j), P1_p(j),P2_p(j), Pc_p(j),Rs_p(j)]=RS_paper(P,h1,h2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [MA_SDMA(j),tou_SDMA(j), P1_SDMA(j),P2_SDMA(j), Pc_SDMA(j),Rs_SDMA(j)]=RS_SDMA(P,h1,h2);
    
    [MA_oma_ZF(j),tou_oma_ZF(j), P1_oma_ZF(j),P2_oma_ZF(j), Pc_oma_ZF(j),Rs_oma_ZF(j)]=RS_oma_ZF(P,h1,h2);
   
    [MA_mul(j),tou_mul(j), P1_mul(j),P2_mul(j), Pc_mul(j),Rs_mul(j)]=RS_multi(P,h1,h2);
    
    [MA_x(j),t_x(j), P1_x(j),P2_x(j), Pc_x(j),Rs_x(j),Rs_noma(j)]=RS_final(P,h1,h2);
    
   
    
        [MA_oma(j),tou_oma(j), P1_oma(j),P2_oma(j), Pc_oma(j),Rs_tdma(j)]=RS_TDMA(P,h1,h2);
    %%%%%%%%%%%%%%%%%%
    
        [MA_np(j),tou_np(j), P1_np(j),P2_np(j), Pc_np(j),Rs_np(j)]=RS_noma_paper(P,h1,h2);
        
        [MA_mp(j),tou_mp(j), P1_mp(j),P2_mp(j), Pc_mp(j),Rs_mp(j)]=RS_mul_paper(P,h1,h2);
    
end
%
%end

figure(1)
rho=[0.001,0.020,0.05,0.1,0.2:0.2:1];%[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1];

hold on
plot(rho,Rs_x,'r-.o','LineWidth',2,'MarkerSize',8)
plot(rho,Rs_p,':bx','LineWidth',2.5,'MarkerSize',4)


plot(rho,Rs_SDMA,'-c*','LineWidth',2,'MarkerSize',8)
plot(rho,Rs_noma,'-g<','LineWidth',2,'MarkerSize',8)
plot(rho,Rs_oma_ZF,':>','LineWidth',2,'MarkerSize',4,'Color',[0 0.4470 0.7410])
plot(rho,Rs_mul,'--*','LineWidth',2,'MarkerSize',4,'Color',[0.4940 0.1840 0.5560])



plot(rho,Rs_tdma,'--','LineWidth',1.7,'MarkerSize',2,'Color',[0.8500 0.3250 0.0980])



 plot(rho,Rs_np,':mx','LineWidth',2,'MarkerSize',8); 
 plot(rho,Rs_mp,'-yx','LineWidth',1,'MarkerSize',4);

plot(rho,Rs_x,'r--o','LineWidth',2,'MarkerSize',7)
plot(rho,Rs_p,':bx','LineWidth',2,'MarkerSize',4)


legend({'Proposed RS','KKT-RS [10]','SDMA','NOMA','OMA','Multicasting','TDMA','kkt-NOMA','kkt-multicasting'},'NumColumns',1')
xlabel('\rho')
ylabel('sum of rate(bits/s/Hz)')
title(' SNR=20dB, channel strength disparity = -5dB')
hold off
%
% hold on
% plot(rho,Rs_noma2,':r*');
% %plot(rho,Rs_noma1,':y^');
% plot(rho,Rs_p,'-.k^');
% plot(rho,Rs_SDMA,'-.b^');
% plot(rho,Rs_x,'-.mo');
% %plot(rho,Rs_x1,'-.g>');
%
% hold off

% figure(2)
% rho=[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1];
%
% hold on
% plot(rho,Rs_noma2,'-.k^');
% plot(rho,Rs_noma1,'-.b^');
%
% hold off