%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%%% Fig 3 (a)
clear all
%clc
nt=2;%number of transmitter antenna
P=10;
i=0;
tic;
o=0;
%for gam_dB=-20:1:0
gam_dB=-5;
gam_dB
i=i+1;
j=0;
% k=0;
for  rho=0.001:0.05:1
    j=j+1;
   
    %rho=0.03;
     rho
   
    [MA_SDMA(i,j),tou_SDMA(i,j), P1_SDMA(i,j),P2_SDMA(i,j), Pc_SDMA(i,j),Rs_SDMA(i,j)]=RS_SDMA(gam_dB,rho,P);

    
     [MA_noma(i,j),tou_noma(i,j), P1_noma(i,j),P2_noma(i,j), Pc_noma(i,j),Rs_noma(i,j)]=RS_noma(gam_dB,rho,P);
    
    [MA_oma(i,j),tou_oma(i,j), P1_oma(i,j),P2_oma(i,j), Pc_oma(i,j),Rs_oma(i,j)]=RS_oma(gam_dB,rho,P);
    
     [MA_oma_ZF(i,j),tou_oma_ZF(i,j), P1_oma_ZF(i,j),P2_oma_ZF(i,j), Pc_oma_ZF(i,j),Rs_oma_ZF(i,j)]=RS_oma_ZF(gam_dB,rho,P);

end


rho=0.001:0.05:1;
hold on
%plot(rho,Rs_x,'r-o',rho,Rs_SDMA,'c:*',rho,Rs_noma,'-+b',rho,Rs_oma,'*:g')
plot(rho,Rs_x(6,:),'r-o','MarkerIndices',1:round(length(rho)/10):length(rho))
plot(rho,Rs_SDMA,'c-.+','MarkerIndices',1:round(length(rho)/10):length(rho))
plot(rho,Rs_noma,'-.b^','MarkerIndices',1:round(length(rho)/10):length(rho))
plot(rho,Rs_oma,'-.g>','MarkerIndices',1:round(length(rho)/10):length(rho))

plot(rho,Rs_oma_ZF,':m<','MarkerIndices',1:round(length(rho)/10):length(rho))
plot(rho,Rs_p(6,:),':k*','MarkerIndices',1:round(length(rho)/10):length(rho))

legend('RS','SDMA','NOMA','OMA without ZF','OMA','pre-RS')
xlabel('rho (-5dB(channel strength disparity), P=10W)')
ylabel('Sum of rate(bits/s/Hz)')
%end
hold off



