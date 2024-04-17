% The code allows to simulate the model in order to obtain the results
% concerning the effects of PARP on VEGF production

clear all
close all

%Model parameters
global K_C k_C th_SH th_H d_C k_cte d_cte k_V d_V  

Parameter_VEGFstudy

%Temporal scaling parameter
kminus1=0.16*10^(-3);

Tmax=4/24*86400*kminus1;     %Adimensional final time
T_real=Tmax/kminus1/86400;   %Dimensional final time

dt=0.0096;   %\hat{dt}=kminus1*dt se dt=1min, allora \hat{dt}=0.0096

%Function accounting for PJ-34 effects on VEGF production, estimated from
%Figure 2H of Marti et al. (2021)
x_Pj=[2/24*86400*kminus1,4/24*86400*kminus1,86400*kminus1];
y_Pj=[0.3494,0.5712,0.6024];
x_v=0:dt:Tmax;
p_Pj=interp1(x_Pj,y_Pj,x_v,'linear');

for i=1:size(p_Pj,2)
     if isnan(p_Pj(i))==1
        p_Pj(i)=y_Pj(1);
    end
end

%Model parameters
InCon_VEGFstudy

X_CelMod=zeros(1,2);
X_CelMod(1,:)=[C_0;Cyto_0];

X_VegfMod=zeros(1,1);
X_VegfMod(1,:)=V_0;

T=0;
step=1;
while T<Tmax
   step=step+1;
  
   %%%Firts step of the time splitting: C, Cyt
   par_1step=[k_C,K_C,th_SH,th_H,d_C,k_cte,d_cte,dt];
   
   K1=RK4_VEGFCell_1step(par_1step,X_CelMod(1,1:2),O2);
   X1=X_CelMod(1,1:2)+K1./2;
   
   K2=RK4_VEGFCell_1step(par_1step,X1,O2);
   X2=X_CelMod(1,1:2)+K2./2;
    
   K3=RK4_VEGFCell_1step(par_1step,X2,O2);
   X3=X_CelMod(1,1:2)+K3;

   K4=RK4_VEGFCell_1step(par_1step,X3,O2);
   X_CelMod(1,1:2)=X_CelMod(1,1:2)+1/6.*(K1+2*K2+2*K3+K4);
   
   clear K1 K2 K3 K4 X1 X2 X3
   C(step)=X_CelMod(1,1);
   Cyto(step)=X_CelMod(1,2);

   %%%Secondo step of the time splitting: VEGF
   par_2step=[k_V,d_V,dt];

   K1=RK4_VEGFCell_2step(par_2step,X_VegfMod(1,1),p_V,p_Pj(step),C(step),O2);
   X1=X_VegfMod(1,1)+K1./2;
   
   K2=RK4_VEGFCell_2step(par_2step,X1,p_V,p_Pj(step),C(step),O2);
   X2=X_VegfMod(1,1)+K2./2;
    
   K3=RK4_VEGFCell_2step(par_2step,X2,p_V,p_Pj(step),C(step),O2);
   X3=X_VegfMod(1,1)+K3;

   K4=RK4_VEGFCell_2step(par_2step,X3,p_V,p_Pj(step),C(step),O2);
   X_VegfMod(1,1)=X_VegfMod(1,1)+1/6.*(K1+2*K2+2*K3+K4);
   
   clear K1 K2 K3 K4 X1 X2 X3
   
   V(step)=X_VegfMod(1,1);

   T=T+dt;
   
   %Saving option for partial results
   if mod(step,1000)==0
   save('VEGFstudy_partial.mat')
   end
   
end 

%Saving option for final results
save('VEGFstudy_final.mat')

%Observed quantity
VEGF_final=V(end);

return
