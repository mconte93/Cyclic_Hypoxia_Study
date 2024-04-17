%Dimensionless values of the parameter for obtaining the simulations shown in Figure 3,
%i.e., the early-stage tumor evolution.

global k1 kminus1 k2 kminus2 
global m b k_R K_R d_R K_1 K_2

%% Angiopoietin module parameter

%Ang1 parameter
k1=0.116*10^(-3);      %1/nM*s 
kminus1=0.16*10^(-3);  %1/s

%Ang2 parameter
k2=0.116*10^(-3);      %1/nM*s  
kminus2=0.35*10^(-3);  %1/s

%pTie2 parameter                                              
m=8;                            
k_R=600;                           
b=m/k_R-10^(-4);                    
K_R=min(max([0.5 2*m b]),1);  
d_R=0.000684;                     
K_1=1.34;              %nM
K_2=3;                 %nM               
