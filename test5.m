%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%%% Fig 3 (a)
% clc
% clear all

nt=2;%number of transmitter antenna
P=1000;%20dB
i=0;

%for SNR=1:2:20
SNR=10;
    P=10^(SNR/10);

%for gam_dB=-10:1:0
    gam_dB=-5;
    i=i+1;
    j=0;
    % k=0;
   for  rho=[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1]

        j=j+1;
        rho

    gam=10^(gam_dB/20);
    theta=acos(1-2*rho);
    Theta=acos(1-2*rho);
    
    h1=1/sqrt(2)*[1;1];
    h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
%         
%         
%         [MA_p(i,j),tou_p(i,j), P1_p(i,j),P2_p(i,j), Pc_p(i,j),Rs_p(i,j)]=RS_paper(gam_dB,rho,P);
%         y=[MA_p(i,j),tou_p(i,j), P1_p(i,j),P2_p(i,j), Pc_p(i,j),Rs_p(i,j)];
        %         [MA_np(i,j),tou_np(i,j), P1_np(i,j),P2_np(i,j), Pc_np(i,j),Rs_np(i,j)]=RS_noma_paper(gam_dB,rho,P);
       % y_n=[MA_np(i,j),tou_np(i,j), P1_np(i,j),P2_np(i,j), Pc_np(i,j),Rs_np(i,j)];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         
%         [MA_SDMA(i,j),tou_SDMA(i,j), P1_SDMA(i,j),P2_SDMA(i,j), Pc_SDMA(i,j),Rs_SDMA(i,j)]=RS_SDMA(gam_dB,rho,P);
%         a=[MA_SDMA(i,j),tou_SDMA(i,j), P1_SDMA(i,j),P2_SDMA(i,j), Pc_SDMA(i,j),Rs_SDMA(i,j)];
%         
%         
%         [MA_oma(i,j),tou_oma(i,j), P1_oma(i,j),P2_oma(i,j), Pc_oma(i,j),Rs_oma(i,j)]=RS_oma(gam_dB,rho,P);
%         [MA_oma_ZF(i,j),tou_oma_ZF(i,j), P1_oma_ZF(i,j),P2_oma_ZF(i,j), Pc_oma_ZF(i,j),Rs_oma_ZF(i,j)]=RS_oma_ZF(gam_dB,rho,P);
%         [MA_mul(i,j),tou_mul(i,j), P1_mul(i,j),P2_mul(i,j), Pc_mul(i,j),Rs_mul(i,j)]=RS_multi(gam_dB,rho,P);
        

            [MA_x(i,j),tou_x(i,j), P1_x(i,j),P2_x(i,j), Pc_x(i,j),Rs_x(i,j),Rs_noma(i,j)]=RS_final(P,h1,h2);

        
        
    end
    
    
%end

figure(1)
rho=[0.001,0.020,0.05,0.1,0.17,0.25:0.125:1]

hold on
plot(rho,Rs_x,'r-.o','LineWidth',2,'MarkerSize',8)
plot(rho,Rs_p,':bx','LineWidth',2.5,'MarkerSize',4)
%plot(rho,Rs_np,':go','LineWidth',0.7,'MarkerSize',8)
%plot(rho,Rs_mp,':m^','LineWidth',0.8,'MarkerSize',5)



plot(rho,Rs_SDMA,'-c*','LineWidth',2,'MarkerSize',8)
plot(rho,Rs_noma,'-g<','LineWidth',2,'MarkerSize',5)
plot(rho,Rs_oma_ZF,':>','LineWidth',2,'MarkerSize',4,'Color',[0 0.4470 0.7410])
plot(rho,Rs_mul,'--*','LineWidth',2,'MarkerSize',4,'Color',[0.4940 0.1840 0.5560])



plot(rho,Rs_oma,'--','LineWidth',1.7,'MarkerSize',2,'Color',[0.8500 0.3250 0.0980])
plot(rho,Rs_x,'r--o','LineWidth',2,'MarkerSize',7)
plot(rho,Rs_p,':bx','LineWidth',2,'MarkerSize',4)


legend({'Proposed RS','KKT-RS [10]','SDMA','NOMA','OMA','Multicasting','TDMA'},'NumColumns',1')
xlabel('\rho')
ylabel('sum of rate(bits/s/Hz)')
title(' SNR of user-1 = 10dB, channel strength disparity = -5dB')
hold off

%title(' SNR=20dB, channel strength disparity=-5dB')

%ylim([0 10])



% 
% figure(2)
% 
% 
% x=0.001:0.0999:1;
% y=-10:1:0;
% [X Y]=meshgrid(x,y);
% 
% contourf(X,Y,Rs_x-Rs_p,10)
% title('Difference from sum rate of the previous method , SNR=20dB')
%  xlabel('rho')
% ylabel('channel strength disparity  [dB]')
% contourcbar