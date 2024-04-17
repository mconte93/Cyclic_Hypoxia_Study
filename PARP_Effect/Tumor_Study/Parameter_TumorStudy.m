%Dimensionless values of the parameter for obtaining the simulations shown in Figure 6,
%i.e., PARP effects on tumor proliferation

global O2_value K_C k_C th_SH th_H Pj nj

O2_value=210; %210 muM=normoxia

%% Cell module parameters

%Tumor cells
k_C=2.15*10^(-2);  
K_C=100;
th_H=50/O2_value; 
th_SH=8.2/O2_value; 

%Therapy
Pj=1;
nj=0.2;