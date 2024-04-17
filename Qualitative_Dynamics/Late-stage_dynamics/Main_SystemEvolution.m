% The code allows to simulate the full model in the advanced stage tumor, based on the interaction 
% between the three different modules.

clear all
close all

%Model parameters
global delta1 alfa1 delta2 alfa2 alfam2 kAng2 gamma1 gamma2 degAng2 dos
global m b k_R K_R d_R K_1 K_2 gamma3
global k_proO2 k_conO2 d_O2 d_C K_C k_C th_SH th_H th_x th_pT k_E K_E k_cte d_cte
global k_V d_V Delta rho gammax gammav nR Kbarx Kbar

RescalePar_xstar2

%Temporal scaling parameter
tau=4;
kminus1=0.16*10^(-3);

Tmax=360/tau*86400*kminus1;  %Adimensional final time
T_real=Tmax/kminus1/86400;   %Dimensional final time

dt=0.0096/tau;   %Adimensional temporal step corresponding to dt=1min

%Initial conditions
InCon_AngVegfCell

X_AngMod=zeros(1,13);
X_AngMod(1,:)=[Ang1_0;Ang1Tie2_4_0;Ang2_0;Ang2Tie2_2_0;Tie2_0;y11_0;y12_0;y21_0;y22_0;z1_0;z2_0;pTie2_0;R_0];

X_CelMod=zeros(1,4);
X_CelMod(1,:)=[E_0;C_0;O2_0;Cyto_0];

X_VegfMod=zeros(1,4);
X_VegfMod(1,:)=[V_0;U_0;B_0;x_0];

%Equilibrium distribution for dimerized VEGRF
y_0=Kbar*gammav*V_0;
ustar_0=(-(1+y_0)+sqrt((1+y_0).^2+8.*nR.*pi.*Delta.^2.*rho.*Kbarx.*y_0))./(4*pi*Delta^2*rho.*Kbarx.*y_0);
xstar(1)=1./(2*gammax).*(nR-ustar_0.*(1+y_0));


T=0;
step=1;
while T<Tmax
   step=step+1;
   
   %%%% First step of the time splitting: Ang1, Ang1Tie2_4, Ang2, Ang2Tie2_2, Tie2, 
   par_1step=[delta1,alfa1,delta2,alfa2,alfam2,gamma1,gamma2,kAng2,degAng2,dt];
   K1=RK4_1step_AngCell_xstar(par_1step,X_AngMod(1,1:5),E(step-1),V(step-1));
   X1=X_AngMod(1,1:5)+K1./2;
   
   K2=RK4_1step_AngCell_xstar(par_1step,X1,E(step-1),V(step-1));
   X2=X_AngMod(1,1:5)+K2./2;
    
   K3=RK4_1step_AngCell_xstar(par_1step,X2,E(step-1),V(step-1));
   X3=X_AngMod(1,1:5)+K3;

   K4=RK4_1step_AngCell_xstar(par_1step,X3,E(step-1),V(step-1));
   X_AngMod(1,1:5)=X_AngMod(1,1:5)+1/6.*(K1+2*K2+2*K3+K4).*tau;

   clear K1 K2 K3 K4 X1 X2 X3
   
   Ang1(step)=X_AngMod(1,1);
   Ang1Tie2_4(step)=X_AngMod(1,2);
   Ang2(step)=X_AngMod(1,3);
   Ang2Tie2_2(step)=X_AngMod(1,4);
   Tie2(step)=X_AngMod(1,5);
   

   %%%% Second step of the time splitting: yij and zi
   par_2step=[alfa1,alfa2,alfam2,gamma1,gamma2,Ang1_0,dos,dt];
   K1=RK4_2step_AngCell_xstar(par_2step,X_AngMod(1,6:11),Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X1=X_AngMod(1,6:11)+K1./2;
   
   K2=RK4_2step_AngCell_xstar(par_2step,X1,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X2=X_AngMod(1,6:11)+K2./2;
    
   K3=RK4_2step_AngCell_xstar(par_2step,X2,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X3=X_AngMod(1,6:11)+K3;

   K4=RK4_2step_AngCell_xstar(par_2step,X3,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X_AngMod(1,6:11)=X_AngMod(1,6:11)+1/6.*(K1+2*K2+2*K3+K4).*tau;
   
   clear K1 K2 K3 K4 X1 X2 X3
   
   y11(step)=X_AngMod(1,6);
   y12(step)=X_AngMod(1,7);
   y21(step)=X_AngMod(1,8);
   y22(step)=X_AngMod(1,9);
   z1(step)=X_AngMod(1,10);
   z2(step)=X_AngMod(1,11);
   

   %%%Third step of the time splitting: pTie2 and R
   par_3step=[m,b,k_R,K_R,d_R,K_1,K_2,alfa1,gamma1,delta1,alfa2,gamma2,delta2,alfam2,Ang1_0,Ang2_0,Tie2_0,gamma3,dos,dt];
   K1=RK4_3step_AngCell_xstar(par_3step,X_AngMod(1,12:13),Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X1=X_AngMod(1,12:13)+K1./2;
   
   K2=RK4_3step_AngCell_xstar(par_3step,X1,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X2=X_AngMod(1,12:13)+K2./2;
    
   K3=RK4_3step_AngCell_xstar(par_3step,X2,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X3=X_AngMod(1,12:13)+K3;

   K4=RK4_3step_AngCell_xstar(par_3step,X3,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step)); 
   X_AngMod(1,12:13)=X_AngMod(12:13)+1/6.*(K1+2*K2+2*K3+K4).*tau;
   
   clear K1 K2 K3 K4 X1 X2 X3

   pTie2(step)=X_AngMod(1,12);
   R(step)=X_AngMod(1,13);


   %%%Fourth step of the time splitting: C, E, O_2, Cyt
   par_4step=[k_C,K_C,K_E,th_pT,k_E,th_SH,th_H,th_x,k_proO2,k_conO2,d_O2,d_C,k_cte,d_cte,xstar(step-1),dt];
   
   K1=RK4_4step_AngCell_xstar(par_4step,X_CelMod(1,1:4),pTie2(step));
   X1=X_CelMod(1,1:4)+K1./2;
   
   K2=RK4_4step_AngCell_xstar(par_4step,X1,pTie2(step));
   X2=X_CelMod(1,1:4)+K2./2;
    
   K3=RK4_4step_AngCell_xstar(par_4step,X2,pTie2(step));
   X3=X_CelMod(1,1:4)+K3;

   K4=RK4_4step_AngCell_xstar(par_4step,X3,pTie2(step));
   X_CelMod(1,1:4)=X_CelMod(1,1:4)+1/6.*(K1+2*K2+2*K3+K4).*tau;
   
   clear K1 K2 K3 K4 X1 X2 X3

   E(step)=X_CelMod(1,1);
   C(step)=X_CelMod(1,2);
   O2(step)=X_CelMod(1,3);
   Cyto(step)=X_CelMod(1,4);


   %%%Fifth step of the time splitting: VEGF and xstar
   par_5step=[k_V,d_V,dt];
   K1=RK4_5step_AngCell_xstar(par_5step,X_VegfMod(1,1),p_V,C(step),O2(step));
   X1=X_VegfMod(1,1)+K1./2;
   
   K2=RK4_5step_AngCell_xstar(par_5step,X1,p_V,C(step),O2(step));
   X2=X_VegfMod(1,1)+K2./2;
    
   K3=RK4_5step_AngCell_xstar(par_5step,X2,p_V,C(step),O2(step));
   X3=X_VegfMod(1,1)+K3;

   K4=RK4_5step_AngCell_xstar(par_5step,X3,p_V,C(step),O2(step));
   X_VegfMod(1,1)=X_VegfMod(1,1)+1/6.*(K1+2*K2+2*K3+K4).*tau;
   
   clear K1 K2 K3 K4 X1 X2 X3
   
   V(step)=X_VegfMod(1,1);


   %%%Evaluation of xstar evolution
   y=Kbar*gammav*V(step);
   ustar=(-(1+y)+sqrt((1+y).^2+8.*nR.*pi.*Delta.^2.*rho.*Kbarx.*y))./(4*pi*Delta^2*rho.*Kbarx.*y);
   xstar(step)=1./(2*gammax).*(nR-ustar.*(1+y));

   T=T+dt;
   
   %Saving option for partial results
   if mod(step,1000)==0
   save('AdvancedStage_partial.mat')
   end
   
end 

%Saving option for final results
save('AdvancedStage_final.mat')

return