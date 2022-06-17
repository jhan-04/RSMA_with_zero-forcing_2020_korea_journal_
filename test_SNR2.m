% clear all
% clc


Nt=2;%nuber of transmitter antenna
M=2;%number of user
SNR=[0:6:30];
%Pt=10^(SNR/10);%transmit signal power
N0=1;%gaussina noise

sigma_e1=0.15; sigma_e2=0.3; sigma_e3=0.3;
beta=sigma_e1^2;
sigma_e=randn(M,1)*0+sigma_e1;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
sample=100;
for i=1:sample
    
    i
    %%%hK
            for k=1:M
                h(:,k)=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh h1
            end
    
    
 
%     ha=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
%     hb=(gam)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
    if norm(h(:,1))>=norm(h(:,2))
        h1(:,i)=h(:,1);
        h2(:,i)=h(:,2);
    else
        h1(:,i)=h(:,2);
        h2(:,i)=h(:,1);
    end
    
    
    
    for j=1:1:length(SNR)
        
        Pt=10^(SNR(j)/10);
        
        [A,B,C,D]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [X,result]=cal_X(Nt,M,Pt,A,B,C,D);
        result_set(1:length(result),i,j)=result;
        %reult/log(2)==GMI
        [GMI_X(i,j)]=cal_GMI_withX(M,A,B,C,D,X);
        [p_max(:,i,j)]=find_p(M,Pt,A,B,C,D,X);norm(p_max(:,i,j))^2;
        [GMI_L(i,j)]=cal_GMI(M,A,B,C,D,p_max(:,i,j));
        %         result(length(result))/log(2)
        %         GMI_L(i,j)
        
        
        %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [X_nc,result_nc]=cal_X_no_common(Nt,M,Pt,A,B,C,D);
        result_nc_set(1:length(result_nc),i,j)=result_nc;
        [GMI_X_nc(i,j)]=cal_GMI_withX(M,A,B,C,D,X_nc);
        [p_nc(:,i,j)]=find_p(M,Pt,A,B,C,D,X_nc);norm(p_nc(:,i,j))^2;
        [GMI_L_nc(i,j)]=cal_GMI(M,A,B,C,D,p_nc(:,i,j));
%                 result_nc(length(result_nc))/log(2)
%                 GMI_L_nc(i,j)
        
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [Xc,PK,result_zf]=cal_X_ZF(Nt,M,Pt,h,N0,sigma_e);
        result_zf_set(1:length(result_zf),i,j)=result_zf;
        pk_zf=h*inv(h'*h);
        for k=1:M
            mag(k)=norm(pk_zf(:,k));
        end
        pk(:,i,j)=reshape(pk_zf.*(sqrt(PK)./mag')',[],1);
        [p_zf(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));norm(p_zf(:,i,j))^2;
        [GMI_zf(i,j)]=cal_GMI(M,A,B,C,D,p_zf(:,i,j));
%                 result_zf(length(result_zf))/log(2)
%                 GMI_zf(i,j)
        
        
        
        
        %%%%%%%%%MRT%%%%%%%%%%
        [Xc,PK,result_mrt]=cal_X_MRT(Nt,M,Pt,h,N0,sigma_e);
        result_mrt_set(1:length(result_mrt),i,j)=result_mrt;
        for k=1:M
            mag(k)=norm(h(:,k));
        end
        pk(:,i,j)=reshape(h.*(sqrt(PK)./mag')',[],1);
        [p_mrt(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));norm(p_mrt(:,i,j))^2;
        [GMI_mrt(i,j)]=cal_GMI(M,A,B,C,D,p_mrt(:,i,j));
%                 result_mrt(length(result_mrt))/log(2)
%                 GMI_mrt(i,j)
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        for m=1:500
            for k=1:M
                hh(:,k)=h(:,k)+sigma_e(k)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
            end
            %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
            [Rs_set(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_max(:,i,j));
            ach_rate_set(i,j,m)=sum(Rs_set(:,i,j));
            
            %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
            [Rs_set_nc(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_nc(:,i,j));
            ach_rate_nc_set(i,j,m)=sum(Rs_set_nc(:,i,j));
            
            %%%%%%%%%Zero_Forcing%%%%%%%%%%
            [Rs_set_zf(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_zf(:,i,j));
            ach_rate_zf_set(i,j,m)=sum(Rs_set_zf(:,i,j));
            
                        %%%%%%%%%MRT%%%%%%%%%%
            [Rs_set_mrt(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_mrt(:,i,j));
            ach_rate_mrt_set(i,j,m)=sum(Rs_set_mrt(:,i,j));
            
        end
        
        
        ach_rate(i,j)=mean(ach_rate_set(i,j,:));
        ach_rate_nc(i,j)=mean(ach_rate_nc_set(i,j,:));
        ach_rate_zf(i,j)=mean(ach_rate_zf_set(i,j,:));
        ach_rate_mrt(i,j)=mean(ach_rate_mrt_set(i,j,:));
        
       fprintf('================sample=%d=======SNR=%d=============\n',i,SNR(j))
        fprintf('RS: GMI=%d, ach_rate=%d,\n NO_RS: GMI=%d, ach_rate=%d \n',GMI_L(i,j),ach_rate(i,j),GMI_X_nc(i,j),ach_rate_nc(i,j))
        fprintf(' RS_ZF: GMI=%d, ach_rate=%d\n  RS_MRT: GMI=%d, ach_rate=%d\n',GMI_zf(i,j),ach_rate_zf(i,j),GMI_mrt(i,j),ach_rate_mrt(i,j))


%         fprintf('RS: GMI=%d, \n NO_RS: GMI=%d \n',GMI_L(i,j),GMI_X_nc(i,j))
%         fprintf(' RS_ZF: GMI=%d\n  RS_MRT: GMI=%d\n',GMI_zf(i,j),GMI_mrt(i,j))
        


        
        
        
        
        
        %===========no information about channel error============================
        [AA,BB,CC,DD]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e*0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [X_e,result_e]=cal_X(Nt,M,Pt,AA,BB,CC,DD);
        result_set_e(1:length(result_e),i,j)=result_e;
        %reult/log(2)==GMI
        [GMI_X_e(i,j)]=cal_GMI_withX(M,AA,BB,CC,DD,X_e);
        [p_max_e(:,i,j)]=find_p(M,Pt,AA,BB,CC,DD,X_e);
        [GMI_L_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_max_e(:,i,j));
        %         result(length(result))/log(2)
        %         GMI_L(i,j)
        
        
        %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [X_nc_e,result_nc_e]=cal_X_no_common(Nt,M,Pt,AA,BB,CC,DD);
        result_nc_set_e(1:length(result_nc_e),i,j)=result_nc_e;
        [GMI_X_nc_e(i,j)]=cal_GMI_withX(M,AA,BB,CC,DD,X_nc_e);
        [p_nc_e(:,i,j)]=find_p(M,Pt,AA,BB,CC,DD,X_nc_e);
        [GMI_L_nc_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_nc_e(:,i,j));
        %         result_nc(length(result_nc))/log(2)
        %         GMI_L_nc(i,j)
        
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [Xc_e,PK_e,result_zf_e]=cal_X_ZF(Nt,M,Pt,h,N0,sigma_e*0);
        result_zf_set_e(1:length(result_zf_e),i,j)=result_zf_e;
        pk_zf_e=h*inv(h'*h);
        for k=1:M
            mag(k)=norm(pk_zf_e(:,k));
        end
        pk_e(:,i,j)=reshape(pk_zf_e.*(sqrt(PK_e)./mag')',[],1);
        [p_zf_e(:,i,j)]=find_p_ZF(M,Pt,AA,BB,CC,DD,Xc_e,pk_e(:,i,j));
        [GMI_zf_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_zf_e(:,i,j));
        %         result_zf(length(result_zf))/log(2)
        %         GMI_zf(i,j)
        
        
                %%%%%%%%%MRT%%%%%%%%%%
        [Xc_e,PK_e,result_mrt_e]=cal_X_MRT(Nt,M,Pt,h,N0,sigma_e*0);
        result_mrt_set_e(1:length(result_mrt_e),i,j)=result_mrt_e;
        for k=1:M
            mag(k)=norm(h(:,k));
        end
        pk(:,i,j)=reshape(h.*(sqrt(PK)./mag')',[],1);
        [p_mrt_e(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));
        [GMI_mrt_e(i,j)]=cal_GMI(M,A,B,C,D,p_mrt_e(:,i,j));
%                 result_mrt_e(length(result_mrt_e))/log(2)
%                 GMI_mrt_e(i,j)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        for m=1:500
            for k=1:M
                hh(:,k)=h(:,k)+sigma_e(k)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
            end
            %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
            [Rs_set_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_max_e(:,i,j));
            ach_rate_set_e(i,j,m)=sum(Rs_set_e(:,i,j));
            
            %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
            [Rs_set_nc_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_nc_e(:,i,j));
            ach_rate_nc_set_e(i,j,m)=sum(Rs_set_nc_e(:,i,j));
            
            %%%%%%%%%Zero_Forcing%%%%%%%%%%
            [Rs_set_zf_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_zf_e(:,i,j));
            ach_rate_zf_set_e(i,j,m)=sum(Rs_set_zf_e(:,i,j));
                                    %%%%%%%%%MRT%%%%%%%%%%
            [Rs_set_mrt_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_mrt_e(:,i,j));
            ach_rate_mrt_set_e(i,j,m)=sum(Rs_set_mrt_e(:,i,j));
            
        end
        
        
        ach_rate_e(i,j)=mean(ach_rate_set_e(i,j,:));
        ach_rate_nc_e(i,j)=mean(ach_rate_nc_set_e(i,j,:));
        ach_rate_zf_e(i,j)=mean(ach_rate_zf_set_e(i,j,:));
        ach_rate_mrt_e(i,j)=mean(ach_rate_mrt_set_e(i,j,:));        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('=======no error information===sample=%d=======SNR=%d=============\n',i,SNR(j))
        fprintf('RS: GMI=%d, ach_rate=%d,\n NO_RS: GMI=%d, ach_rate=%d \n',GMI_L_e(i,j),ach_rate_e(i,j),GMI_X_nc_e(i,j),ach_rate_nc_e(i,j))
        fprintf(' RS_ZF: GMI=%d, ach_rate=%d\n  RS_MRT: GMI=%d, ach_rate=%d\n',GMI_zf_e(i,j),ach_rate_zf_e(i,j),GMI_mrt_e(i,j),ach_rate_mrt_e(i,j))
        
        

        %%%%%%%%%%NRS_NOMA
        %%%%%%%%%NRS_NOMA
        [MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_sdr_perfect(Pt,h1(:,i),h2(:,i),rho);
        %KKT
        [MA_Rik(i,j),t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),Rs_Rik(i,j),fc_Rik(:,i,j),t_Nik(i,j),Rs_Nik(i,j),fc_Nik(:,i,j)]=RSNOMA_kkt_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Rpk(i,j),t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),Rs_Rpk(i,j),fc_Rpk(:,i,j),t_Npk(i,j),Rs_Npk(i,j),fc_Npk(:,i,j)]=RSNOMA_kkt_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        [MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),fc_Si(:,i,j)]=SDMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),fc_Sp(:,i,j)]=SDMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        %%%%%%%OMA
        [MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),fc_Oi(:,i,j)]=OMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),fc_Op(:,i,j)]=OMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        %%%%%%multicasting
        [MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),fc_mi(:,i,j)]=multi_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),fc_mp(:,i,j)]=multi_sdr_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%KKT
        [MA_mik(i,j),t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),Rs_mik(i,j),fc_mik(:,i,j)]=multi_kkt_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_mpk(i,j),t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),Rs_mpk(i,j),fc_mpk(:,i,j)]=multi_kkt_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        
        %%%%%%%%%%%%%%%%%%caculate GMI
      
        Rs__Ri(i,j)=cal_GMI_1(Pt,t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),fc_Ri(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Rp(i,j)=cal_GMI_1(Pt,t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),fc_Rp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Ni(i,j)=cal_GMI_1(Pt,t_Ni(i,j), t_Ni(i,j)*Pt,0, (1-t_Ni(i,j))*Pt,fc_Ni(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Np(i,j)=cal_GMI_1(Pt,t_Np(i,j), t_Np(i,j)*Pt,0, (1-t_Np(i,j))*Pt,fc_Np(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        %%%KKT
        Rs__Rik(i,j)=cal_GMI_1(Pt,t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),fc_Rik(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Rpk(i,j)=cal_GMI_1(Pt,t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),fc_Rpk(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Nik(i,j)=cal_GMI_1(Pt,t_Nik(i,j), t_Nik(i,j)*Pt, 0, (1-t_Nik(i,j))*Pt,fc_Nik(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Npk(i,j)=cal_GMI_1(Pt,t_Npk(i,j), t_Npk(i,j)*Pt, 0, (1-t_Npk(i,j))*Pt,fc_Npk(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        Rs__Si(i,j)=cal_GMI_1(Pt,t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),fc_Si(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Sp(i,j)=cal_GMI_1(Pt,t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),fc_Sp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        %%%%%%%%%%%%%%%%%%%OMA
        Rs__Oi(i,j)=cal_GMI_1(Pt,t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),fc_Oi(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Op(i,j)=cal_GMI_1(Pt,t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),fc_Op(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        %%%%%%%%%%%%%%%multi
        Rs__mi(i,j)=cal_GMI_1(Pt,t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),fc_mi(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__mp(i,j)=cal_GMI_1(Pt,t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),fc_mp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %%KKT
        Rs__mik(i,j)=cal_GMI_1(Pt,t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),fc_mik(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__mpk(i,j)=cal_GMI_1(Pt,t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),fc_mpk(:,i,j),h1(:,i),h2(:,i),rho,beta);
   

        
        
        
        fprintf('=======sample=%d===================SNR=%SNR===================\n\n',i,SNR)
        
        fprintf('RS_SDR_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),Rs__Ri(i,j))
        fprintf('NOMA_SDR_im: t=%d, Rs=%d Rs2=%d\n\n',t_Ni(i,j),Rs_Ni(i,j),Rs__Ni(i,j))
        fprintf('RS_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),Rs__Rp(i,j))
        fprintf('NOMA_SDR: t=%d, Rs=%d GMI=%d\n\n',t_Np(i,j),Rs_Np(i,j),Rs__Np(i,j))
        %%%
        fprintf('RS_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Rik(i,j),t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),Rs_Rik(i,j),Rs__Rik(i,j))
        fprintf('NOMA_KKT_im: t=%d, Rs=%d Rs2=%d\n\n',t_Nik(i,j),Rs_Nik(i,j),Rs__Nik(i,j))
        fprintf('RS_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Rpk(i,j),t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),Rs_Rpk(i,j),Rs__Rpk(i,j))
        fprintf('NOMA_KKT: t=%d, Rs=%d GMI=%d\n\n',t_Npk(i,j),Rs_Npk(i,j),Rs__Npk(i,j))
        
        fprintf('SDMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),Rs__Si(i,j))
        fprintf('SDMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),Rs__Sp(i,j))
        
        fprintf('OMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),Rs__Oi(i,j))
        fprintf('OMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),Rs__Op(i,j))
        
        fprintf('multi_SDR_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),Rs__mi(i,j))
        fprintf('multi_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),Rs__mp(i,j))
        
        fprintf('multi_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mik(i,j),t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),Rs_mik(i,j),Rs__mik(i,j))
        fprintf('multi_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_mpk(i,j),t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),Rs_mpk(i,j),Rs__mpk(i,j))
        











        
        
        
    end
    
    
    
    
    
end

%i=1:sample;
hold on
plot(SNR,mean(GMI_L),'-r','LineWidth',2,...
    'MarkerSize',8)
plot(SNR,mean(GMI_L_nc),'b-','LineWidth',2,...
    'MarkerSize',8)
plot(SNR,mean(GMI_zf),'g-','LineWidth',2,...
    'MarkerSize',10)
plot(SNR,mean(GMI_mrt),'c-','LineWidth',2,...
    'MarkerSize',10)

plot(SNR,mean(ach_rate),'r--o','LineWidth',2,...
    'MarkerSize',8)
plot(SNR,mean(ach_rate_nc),'b--o','LineWidth',2,...
    'MarkerSize',7)
plot(SNR,mean(ach_rate_zf),'g--o','LineWidth',2,...
    'MarkerSize',7)
plot(SNR,mean(ach_rate_mrt),'c--o','LineWidth',2,...
    'MarkerSize',6)

plot(SNR,mean(ach_rate_e),':r>','LineWidth',2,...
    'MarkerSize',5)
plot(SNR,mean(ach_rate_nc_e),'b:^','LineWidth',2,...
    'MarkerSize',4)
plot(SNR,mean(ach_rate_zf_e),'g:<','LineWidth',2,...
    'MarkerSize',4)
plot(SNR,mean(ach_rate_mrt_e),'c:<','LineWidth',2,...
    'MarkerSize',4)

legend('GMI RS','GMI no-RS','GMI RS-ZF','GMI RS-MRT','rate RS','rate no-RS','rate RS-ZF','rate RS-MRT','rate RS no info','rate no-RS  no info','rate RS-ZF no info','rate RS-MRT no info')
xlabel('SNR')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(M),', \sigma_{e}=',num2str(sigma_e(1)),', sample=',num2str(sample)])

% figure(2)
% hold on
% plot(SNR,mean(GMI_L_e),'-.r')
% plot(SNR,mean(GMI_L_nc_e),'b-.')
% plot(SNR,mean(GMI_zf_e),'g-.')
% plot(SNR,mean(GMI_mrt_e),'c-.')
% 
% legend('GMI RS','GMI no-RS','GMI RS-ZF','GMI RS-MRT')
% xlabel('SNR')
% ylabel('rate (bits/s/Hz)')
% title(['Nt=',num2str(Nt),', M=',num2str(M),', \sigma_{e}=',num2str(0),', sample=',num2str(sample)])


%pp=sqrt(Pt)*p_max/norm(p_max);

