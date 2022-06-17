%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%%% Fig 3 (a)
clc
clear all

Nt=2;%number of transmitter antenna

%for SNR=1:2:20
SNR=10;
P=10^(SNR/10);
i=0;
for gam_dB=-10:1:0
    gam_dB

    j=0;
    i=i+1;
    for NUM=0:500
        j=j+1
        
        gam=10^(gam_dB/20);
        h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
        h2=gam*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%     gam=10^(gam_dB/20);
%     theta=acos(1-2*rho);
%     Theta=acos(1-2*rho);
%     
%     h1=1/sqrt(2)*[1;1];
%     h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rho(j)=1-abs(h1'/norm(h1)*h2/norm(h2))^2;
        
        
        
        [MA_x(j),t_x(j), P1_x(j),P2_x(j), Pc_x(j),Rs_x(j),Rs_noma(j)]=RS_final(P,h1,h2);
        
        
        [MA_p(j),tou_p(j), P1_p(j),P2_p(j), Pc_p(j),Rs_p(j)]=RS_paper(P,h1,h2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [MA_SDMA(j),tou_SDMA(j), P1_SDMA(j),P2_SDMA(j), Pc_SDMA(j),Rs_SDMA(j)]=RS_SDMA(P,h1,h2);
        
        [MA_oma_ZF(j),tou_oma_ZF(j), P1_oma_ZF(j),P2_oma_ZF(j), Pc_oma_ZF(j),Rs_oma_ZF(j)]=RS_oma_ZF(P,h1,h2);
        
        [MA_mul(j),tou_mul(j), P1_mul(j),P2_mul(j), Pc_mul(j),Rs_mul(j)]=RS_multi(P,h1,h2);
        
        
        
        
        [MA_oma(j),tou_oma(j), P1_oma(j),P2_oma(j), Pc_oma(j),Rs_tdma(j)]=RS_TDMA(P,h1,h2);
        
        
    end
    mean(rho)
    MINE(i)=mean(Rs_x);
    Paper(i)=mean(Rs_p);
    
    SDMA(i)=mean(Rs_SDMA);
    NOMA(i)=mean(Rs_noma);
    OMA(i)=mean(Rs_oma_ZF);
    MULTI(i)=mean(Rs_mul);
    
    TDMA(i)=mean(Rs_tdma);
    
end

figure(1)
gam_dB=-10:1:0

hold on
plot(gam_dB,MINE,'r-.o','LineWidth',2,'MarkerSize',8);
plot(gam_dB,Paper,':bx','LineWidth',2.5,'MarkerSize',4);

plot(gam_dB,SDMA,'-c*','LineWidth',2,'MarkerSize',8);
plot(gam_dB,NOMA,'-g<','LineWidth',2,'MarkerSize',5);
plot(gam_dB,OMA,':>','LineWidth',2,'MarkerSize',4,'Color',[0 0.4470 0.7410])
plot(gam_dB,MULTI,'--*','LineWidth',2,'MarkerSize',4,'Color',[0.4940 0.1840 0.5560])

plot(gam_dB,TDMA,'-y*','LineWidth',2,'MarkerSize',4)


%plot(rho,Rs_noma1,':y^');


%plot(rho,Rs_x1,'-.g>');


legend({'Proposed RS','KKT-RS [1]','SDMA','NOMA','OMA','Multicasting','TDMA'},'NumColumns',1')
xlabel('\gamma_{dB}')
ylabel('sum of rate(bits/s/Hz)')
title(' SNR=10dB')
hold off
