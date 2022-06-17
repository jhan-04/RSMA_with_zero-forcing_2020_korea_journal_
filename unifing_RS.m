%Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
clear all
clc
nt=2;%number of transmitter antenna
P=10;
i=0;
for gam_dB=-20:1:-10
    gam_dB
    gam=10^(gam_dB/20);
    i=i+1;
    j=0;
    k=0;
    for  rho=0:0.2:1
        j=j+1;
        rho;
        
        theta=acos(1-2*rho);
        Theta(i,j)=acos(1-2*rho);
        
        
        h1=1/sqrt(2)*[1;1];
        h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
        h1_=h1;
        h2_=h2/gam;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A=h1_*h1_';
        B=h2_*h2_';
        %%%%%%%%%%%%%%%%%%%%
        cvx_begin quiet
        variable F(2,2) hermitian
        variable x(1)
        minimize(-x)
        subject to
        trace(A*F)>=x;
        trace(B*F)>=x;
        trace(F)==1;
        F== semidefinite(2);
        cvx_end
        %%%%%%%%%%%%%%%%
        
        [U,W,Z] = svds(P);%V=U*W*Z'
        obj_max=-100;
        n=length(W);
        for m=1:500
            r=sqrt(1/2)*(randn(2,1)+1i*randn(2,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
            v_=U*sqrt(W)*r;
            v_=v_/sqrt(v_'*v_);
            obj=min(v_'*A*v_,v_'*B*v_);
            if obj>=obj_max
                obj_max=obj;
                v_max=v_;
            end
        end
        f_c=v_max;
        
        
        %%%%%%%P1 ~=0,P2 ~=0
        
        Gam=(1/rho)*(1/norm(h2)^2-1/norm(h1)^2);
        b=norm(h1)^2*rho*P/2;
        a=1+(Gam/P)*b;
        d=norm(h2)^2*rho*P/2-abs(h2'*f_c)^2*P;
        c=1-(Gam/P)*d+abs(h2'*f_c)^2*(P-Gam);
        
        
        
        
        t12(i,j)=min(-a/(2*b)-c/(2*d),1);
        
        Rs_12(i,j)=log2(a*c+(a*d+b*c)*t12(i,j)+b*d*t12(i,j)^2);
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
        
        
        if (t12(i,j)*P >= Gam) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.
            Rs_x(i,j)=Rs_12(i,j);
            tou_x(i,j)=t12(i,j);
            P1_x=t12(i,j)*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
            P2_x=t12(i,j)*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
            Pc_x=(1-t12(i,j))*P;
            
        else
            v=0;
            for t=0:1/100:1
                v=v+1;
                h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
                h2_=h2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                A=h1_*h1_';
                B=h2_*h2_';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cvx_begin quiet
                variable F(2,2) hermitian
                variable x(1)
                minimize(-x)
                subject to
                trace(A*F)>=x;
                trace(B*F)>=x;
                trace(F)==1;
                F== semidefinite(2);
                cvx_end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                [U,W,Z] = svds(F);%V=U*W*Z'
                obj_max=-100;
                n=length(W);
                for m=1:500
                    r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
                    v_=U*sqrt(W)*r;
                    v_=v_/sqrt(v_'*v_);
                    obj=min(v_'*A*v_,v_'*B*v_);
                    if obj>=obj_max
                        obj_max=obj;
                        v_max=v_;
                    end
                end
                f_c=v_max;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                b1=P*abs(h2'*f_c)^2;
                a1=norm(h1)^2*rho*P;
                R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
                
            end
            
            
            t=0:1/100:1;
            Rs_1(i,j)=max(R);
            k=find(Rs_1(i,j)==R);
            t1(i,j)=t(k);
            
            tou_x(i,j)=t1(i,j);
            P2_x=0;
            P1_x=P*t1(i,j);
            Pc_x=(1-t1(i,j))*P;
            Rs_x(i,j)=Rs_1(i,j);
            
            
            
            
        end
        
        
        if tou_x(i,j)==1
            MA_x(i,j)=1;%SDMA
        end
        if (P1_x~=0)&&(P2_x==0)&&(Pc_x~=0)
            MA_x(i,j)=2;%NOMA
        end
        if (P1_x~=0)&&(P2_x==0)&&(Pc_x==0)
            MA_x(i,j)=3;%OMA
        end
        if (P2_x==0)&&(P1_x==0)%tou(i,j)==0
            MA_x(i,j)=4;%multicasting
        end
        if (P1_x~=0)&&(P2_x~=0)&&(Pc_x~=0)
            MA_x(i,j)=5;%RS
        end
        
        
        
        
    end
    
end



figure(1)
x=0:0.2:1;
y=-20:1:0;
[X Y]=meshgrid(x,y);
% contourf(X,Y,tou,'ShowText','on')


figure(1)
contourf(X,Y,MA_x,5)

xlabel('rho')
ylabel('channel strength disparity  [dB]')
figure(2)
contourf(X,Y,tou_x)
contourcbar
xlabel('rho')
ylabel('channel strength disparity  [dB]')

%
