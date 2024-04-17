% The code allows to simulate the angiopotin module for the results in Fig. 2

clear all
close all

%Model parameters
global k1 kminus1 k2 kminus2 dos 
global m b k_R K_R d_R K_1 K_2

%Conversion of the initial concentration of angiopoietin from ng/ml to nM: 50ng/ml=0.179 nM 
%Ang_{1,0}=100 ng/ml=2*0.179;
%Ang_{2,0}=[0-4000]ng/ml=dos*0.179 
%dos" varies in [0,10,20,40,80] to obtain the results shown in Figure 2

dos=10;

Parameters_Angio

%Temporal parameters
Tmax=3600*0.5; %Dimensional final time
dt=10^(0);   %Adimensional temporal step corresponding to dt=1sec

%Initial conditions
InCon_Ang

X_AngMod=zeros(1,13);
X_AngMod(1,:)=[Ang1_0;Ang1Tie2_4_0;Ang2_0;Ang2Tie2_2_0;Tie2_0;y11_0;y12_0;y21_0;y22_0;z1_0;z2_0;pTie2_0;R_0];


T=0;
step=1;
while T<Tmax
   step=step+1;
   
   %%%% First step of the time splitting: Ang1, Ang1Tie2_4, Ang2, Ang2Tie2_2, Tie2 
   par_1step=[k1,kminus1,k2,kminus2,dt]; 
   K1=RK4_1step_Ang(par_1step,X_AngMod(1,1:5));
   X1=X_AngMod(1,1:5)+K1./2;
   
   K2=RK4_1step_Ang(par_1step,X1);
   X2=X_AngMod(1,1:5)+K2./2;
    
   K3=RK4_1step_Ang(par_1step,X2);
   X3=X_AngMod(1,1:5)+K3;

   K4=RK4_1step_Ang(par_1step,X3);
   X_AngMod(1,1:5)=X_AngMod(1,1:5)+1/6.*(K1+2*K2+2*K3+K4);

   clear K1 K2 K3 K4 X1 X2 X3
   
   Ang1(step)=X_AngMod(1,1);
   Ang1Tie2_4(step)=X_AngMod(1,2);
   Ang2(step)=X_AngMod(1,3);
   Ang2Tie2_2(step)=X_AngMod(1,4);
   Tie2(step)=X_AngMod(1,5);

   %%%% Second step of the time splitting: yij and zi
   par_2step=[k1,kminus1,k2,kminus2,Ang1_0,Ang2_0,dt];
   K1=RK4_2step_Ang(par_2step,X_AngMod(1,6:11),Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X1=X_AngMod(1,6:11)+K1./2;
   
   K2=RK4_2step_Ang(par_2step,X1,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X2=X_AngMod(1,6:11)+K2./2;
    
   K3=RK4_2step_Ang(par_2step,X2,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   X3=X_AngMod(1,6:11)+K3;

   K4=RK4_2step_Ang(par_2step,X3,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step);
   
   X_AngMod(1,6:11)=X_AngMod(1,6:11)+1/6.*(K1+2*K2+2*K3+K4);

   clear K1 K2 K3 K4 X1 X2 X3
   
   y11(step)=X_AngMod(1,6);
   y12(step)=X_AngMod(1,7);
   y21(step)=X_AngMod(1,8);
   y22(step)=X_AngMod(1,9);
   z1(step)=X_AngMod(1,10);
   z2(step)=X_AngMod(1,11);
   
   %%%Third step of the time splitting: pTie2 and R
   par_3step=[m,b,k_R,K_R,d_R,K_1,K_2,k2,kminus2,k1,kminus1,Ang1_0,Ang2_0,dt];
   K1=RK4_3step_Ang(par_3step,X_AngMod(1,12:13),Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X1=X_AngMod(1,12:13)+K1./2;
   
   K2=RK4_3step_Ang(par_3step,X1,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X2=X_AngMod(1,12:13)+K2./2;
    
   K3=RK4_3step_Ang(par_3step,X2,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X3=X_AngMod(1,12:13)+K3;

   K4=RK4_3step_Ang(par_3step,X3,Tie2(step),Ang1(step),Ang2(step),z1(step),z2(step),Ang2Tie2_2(step),Ang1Tie2_4(step));
   X_AngMod(1,12:13)=X_AngMod(12:13)+1/6.*(K1+2*K2+2*K3+K4);
   
   clear K1 K2 K3 K4 X1 X2 X3

   pTie2(step)=X_AngMod(1,12);
   R(step)=X_AngMod(1,13);
   
   T=T+dt;
  
end 

%Saving option for final results
save('Angiopoietin_study.mat')

%Observed values
ratio=R(end)/R_contr;

