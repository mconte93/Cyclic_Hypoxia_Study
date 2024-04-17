%Adimensional initial conditions

global O2_value 

%% Cellular module
C_0=1;
Cyto_0=0;
O2_0=10/O2_value;

C(1)=C_0;
Cyto(1)=Cyto_0;
O2=O2_0;
%% VEGF module
%Function for oxygen-dependent VEGF expression - estimated from Steinbrech
%et al. (1999)

x_O2=[10,50,100,150,210]./O2_value;
y_V=[3.98,2.29,1.197,0.704,1];
p_V=polyfit(x_O2,y_V,4);

V_0=1;
V(1)=V_0;

