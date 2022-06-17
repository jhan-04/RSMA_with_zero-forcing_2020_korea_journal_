%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%% Fig 3 (a)
clc
clear all

Nt=2;%number of transmitter antenna

% %for SNR=1:2:20
% SNR=5;
% P=10^(SNR/10);
i=0;

MINE=0;
Paper=0;
%%%%%%%
SDMA=0;
NOMA=0;
OMA=0;
MULTI=0;

TDMA=0;
%%%%%%
NOMA_P=0;
MULTI_P=0;


for NUM=1:100
    
    j=0;
    NUM
    
    i=i+1;
    h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
    h2_s=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);

% 
%     h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
%     h2_s=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
%     h1=h1/norm(h1);%Rayleigh
%     h2_s= h2_s/norm(h2_s);


    
    rho(i)=1-abs(h1'/norm(h1)*h2_s/norm(h2_s))^2;
  %for SNR=0:2:10
    SNR=20;
      P=10^(SNR/10);
    
    for gam_dB=-10:2:0
        %gam_dB=-5;
        
        j=j+1
        
        gam=10^(gam_dB/20);
        h2=gam*h2_s;
        
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
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [MA_np(j),tou_np(j), P1_np(j),P2_np(j), Pc_np(j),Rs_np(j)]=RS_noma_paper(P,h1,h2);
        
        [MA_mp(j),tou_mp(j), P1_mp(j),P2_mp(j), Pc_mp(j),Rs_mp(j)]=RS_mul_paper(P,h1,h2);

        
     end
    MINE=MINE+Rs_x;
    Paper=Paper+(Rs_p);
    
    SDMA=SDMA+(Rs_SDMA);
    NOMA=NOMA+(Rs_noma);
    OMA=OMA+(Rs_oma_ZF);
    MULTI=MULTI+(Rs_mul);
    
    TDMA=TDMA+(Rs_tdma);
    
    NOMA_P=NOMA_P+(Rs_np);
    MULTI_P=MULTI_P+(Rs_mp);
    
end


    MINE=MINE/NUM;
    Paper=Paper/NUM;
    
    SDMA=SDMA/NUM;
    NOMA=NOMA/NUM;
    OMA=OMA/NUM;
    MULTI=MULTI/NUM;
    
    TDMA=TDMA/NUM;
    
% 
 NOMA_P= NOMA_P/NUM;
 MULTI_P=MULTI_P/NUM;



figure(1)
gam_dB=-10:2:0;


hold on
plot(gam_dB,MINE,'r-.o','LineWidth',2,'MarkerSize',8);
plot(gam_dB,Paper,':bx','LineWidth',2.5,'MarkerSize',4);

plot(gam_dB,SDMA,'-c*','LineWidth',2,'MarkerSize',8);
plot(gam_dB,NOMA,'-g<','LineWidth',2,'MarkerSize',5);
plot(gam_dB,OMA,':>','LineWidth',2,'MarkerSize',4,'Color',[0 0.4470 0.7410])
plot(gam_dB,MULTI,'--*','LineWidth',2,'MarkerSize',4,'Color',[0.4940 0.1840 0.5560])


plot(gam_dB,TDMA,'--','LineWidth',1.7,'MarkerSize',2,'Color',[0.8500 0.3250 0.0980])




 plot(gam_dB,NOMA_P,':mx','LineWidth',2,'MarkerSize',4); 
 plot(gam_dB,MULTI_P,'-yx','LineWidth',1,'MarkerSize',4);

plot(gam_dB,MINE,'r-.o','LineWidth',2,'MarkerSize',8);
plot(gam_dB,Paper,':bx','LineWidth',2.5,'MarkerSize',4);

%ylim([1,4])

legend({'Proposed RS','KKT-RS [10]','SDMA','NOMA','OMA','Multicasting','TDMA','kkt-NOMA','kkt-multicasting'},'NumColumns',2')
title('SNR of user-1 = 20dB')
ylabel('sum of rate(bits/s/Hz)')
xlabel('channel strength disparity \gamma_{dB}')
hold off

