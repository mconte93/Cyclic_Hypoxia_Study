%Adimensional initial conditions
 
%Angiopoietin module   
Ang1_0=1;           
Ang1Tie2_4_0=1;     
Ang2_0=1;           
Ang2Tie2_2_0=1;     
Tie2_0=1;           
y11_0=-0.001;
y12_0=-0.001;
y21_0=-0.001;
y22_0=-0.001;
z1_0=-0.001;
z2_0=-0.001;
pTie2_0=1;          
R_0=0.17/2;       

Ang1(1)=Ang1_0;
Ang1Tie2_4(1)=Ang1Tie2_4_0;
Ang2(1)=Ang2_0;
Ang2Tie2_2(1)=Ang2Tie2_2_0;
Tie2(1)=Tie2_0;
y11(1)=y11_0;
y12(1)=y12_0;
y21(1)=y21_0;
y22(1)=y22_0;
z1(1)=z1_0;
z2(1)=z2_0;
pTie2(1)=pTie2_0;
R(1)=pTie2(1)/Tie2(1);

%Cellular module
E_0=1;%  
C_0=1;
Cyto_0=0;
O2_0=1;

E(1)=E_0;
C(1)=C_0;
Cyto(1)=Cyto_0;
O2(1)=O2_0;

%VEGF module
V_0=1;
U_0=1;
B_0=1;
x_0=1;

V(1)=V_0;
U(1)=U_0;
B(1)=B_0;
x(1)=x_0;

%Function for oxygen-dependent VEGF expression - estimated from Steinbrech
%et al. (1999)
global O2_value 
x_O2=[10,50,100,150,210]./O2_value;
y_V=[3.98,2.29,1.197,0.704,1];
p_V=polyfit(x_O2,y_V,4);

 