% The code allows to simulate the model in order to obtain the results
% concerning the effects of PARP on EC proliferation

clear all
close all

%Model parameters
global k_E Pj nPj K_E

Parameter_ECstudy 

%Temporal parameters
Tmax=3600*24; %Dimensional final time
dt=10^(0);   %Adimensional temporal step corresponding to dt=1sec

%Initial condition (estimated from the data in Rajesh et al. (2006))
E=zeros(Tmax,1);
E_0=0.065;
E(1)=E_0;

X_CelMod=E_0;

T=0;
step=1;
while T<Tmax
   step=step+1;
   
   par_4step=[k_E,K_E,Pj,nPj,dt];
   K1=RK4_1step_ECstudy(par_4step,X_CelMod);
   X1=X_CelMod+K1./2;
   
   K2=RK4_1step_ECstudy(par_4step,X1);
   X2=X_CelMod+K2./2;
    
   K3=RK4_1step_ECstudy(par_4step,X2);
   X3=X_CelMod+K3;

   K4=RK4_1step_ECstudy(par_4step,X3);
   X_CelMod=X_CelMod+1/6.*(K1+2*K2+2*K3+K4);
   
   clear K1 K2 K3 K4 X1 X2 X3
   E(step)=X_CelMod;
      
   T=T+dt;
  
end 

%Saving option for final results
save('ECstudy_final.mat')

%Observed values
EC_density=E(end);

return

