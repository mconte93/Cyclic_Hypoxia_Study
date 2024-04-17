%Dimensionless values of the parameter for obtaining the simulations shown in Figure 5,
%i.e., PARP effects on VEGF production

global K_C k_C th_SH th_H O2_value d_C k_cte d_cte k_V d_V Pj 
 
O2_value=10; %10 muM=severe hypoxia

%% Cell module parameters

%Tumor cells parameter
k_C=0.00625;
K_C=3;
th_H=50/O2_value; 
th_SH=8.2/O2_value; 

%Cytokine parameter
d_C=0.0125*10^(-2);
k_cte=0.7733*10^(-3);
d_cte=0.3866*10^(-3) ;

%% VEGF module parameter

% %VEGF
k_V=1.9332;
d_V=1.135;

%Therapy
Pj=10;