syms h2 [2,1] 
syms h1 [2,1] 
syms P rho 
syms k  
f=norm(h2)^2*(norm(h1)^2*rho*P+1-k);
g=norm(h1)^2+k*norm(h2)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho);
d_g=diff(g,k);
d_f=diff(f,k);
dd_g=diff(g,k,2);
dd_f=diff(f,k,2);
der_Rgf=(1+f/g)+k*(d_f*g-f*d_g)/g^2;
der2_Rgf=2*(d_f/g-(f*d_g)/g^2)+k*(dd_f/g-(d_f*d_g)/g^2)-k*{(d_f*d_g+f*dd_g)/g^2-2*f*d_g^2/g^3};




Rk_nolog=(k)*(1+norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
der_Rk_nolog=diff(Rk_nolog,k);
der2_Rkk_nolog=diff(Rk_nolog,k,2);