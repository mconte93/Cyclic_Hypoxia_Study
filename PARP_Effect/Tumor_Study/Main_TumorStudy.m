% The code allows to simulate the model in order to obtain the results
% concerning the effects of PARP on tumor proliferation

clear all
close all
 
%Model parameters
global O2_value K_C k_C th_SH th_H Pj nj

Parameter_TumorStudy

%Temporal scaling parameters
tau=4;

kminus1=0.16*10^(-3);
Tmax=72/24*86400*kminus1;   %Adimensional final time
T_real=Tmax/kminus1/86400;  %Dimensional final time

dt=0.0096/tau;              %Adimensional temporal step corresponding to dt=0.25sec

%Initial condition (estimated from the data)
O2=210/O2_value;
C_0=1.805;
X_Cel(1)=C_0;
C(1)=C_0;

T=0;
step=1;

while T<Tmax
   step=step+1;
   
   %%%Fourth step of the time splitting: cellular module
   par=[k_C,K_C,th_SH,th_H,Pj,nj,dt];
   K1=RK4_1step_Cell(par,X_Cel,O2);
   X1=X_Cel+K1./2;
   
   K2=RK4_1step_Cell(par,X1,O2);
   X2=X_Cel+K2./2;
    
   K3=RK4_1step_Cell(par,X2,O2);
   X3=X_Cel+K3;

   K4=RK4_1step_Cell(par,X3,O2);
   X_Cel=X_Cel+1/6.*(K1+2*K2+2*K3+K4).*tau;
   
   clear K1 K2 K3 K4 X1 X2 X3
   C(step)=X_Cel;

   T=T+dt;
   
   %Saving option for partial results
   if mod(step,1000)==0
   save('TumorStudy_partial.mat')
   end
   
end 
toc

%Saving option for final results
save('TumorStudy_final.mat')

%Observed values
Tumor_density=C(end);

return
