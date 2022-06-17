%%%%%%%%%%%%%  Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%%%%%%%%%%%%% Fig 3 (a)
clear all
clc
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
    k=0;
    for  rho=0.1:0.1:1
        tic;
        o=o+1;
        j=j+1;         
        rho
        [MA_p(i,j),tou_p(i,j), P1_p(i,j),P2_p(i,j), Pc_p(i,j),Rs_p(i,j)]=RS_paper(gam_dB,rho,P);
       reult_paper(o,:)=[MA_p(i,j),tou_p(i,j), P1_p(i,j),P2_p(i,j), Pc_p(i,j),Rs_p(i,j)]
%        [MA_x(i,j),tou_x(i,j), P1_x(i,j),P2_x(i,j), Pc_x(i,j),Rs_x(i,j)]=RS_summaxrate(gam_dB,rho,P);
%         reul_mine(o,:)=[MA_x(i,j),tou_x(i,j), P1_x(i,j),P2_x(i,j), Pc_x(i,j),Rs_x(i,j)]
        toc;
    end
    
%end
toc;


x=0.001:0.05:1;
y=-6:1:0;
[X Y]=meshgrid(x,y);

contourf(X,Y,Rs_x([5:11],:)-Rs_p([5:11],:),8)

contourcbar
title('Difference from sum rate of the previous method')
 xlabel('rho')
ylabel('channel strength disparity  [dB]')
% 
% contourf(X,Y,tou,'ShowText','on')
% 
% 
% figure(1)
% 
% 
% contourf(X,Y,MA_p,5)
%contourf(X,Y,tou_p)
%contourcbar
% 
% figure(2)
% contourf(X,Y,tou_x)
% contourcbar
% 
% xlabel('rho')
% ylabel('channel strength disparity  [dB]')
% 
% figure(3)
% contourf(X,Y,MA_p,5)
% %
% figure(4)
% contourf(X,Y,MA_x,5)
% 
% %contourf(X,Y,MA_x)
% 
% 
% xlabel('rho')
% ylabel('channel strength disparity  [dB]')
% 
% 
% % %
% % figure(5)
% % contourf(X,Y,Rs_x-Rs_p,10)
% % contourcbar