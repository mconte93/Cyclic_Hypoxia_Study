%Dimensionless values of the parameter for obtaining the simulations shown in Figure 3,
%i.e., the early-stage tumor evolution.

global delta1 alfa1 delta2 alfa2 alfam2 kAng2 gamma1 gamma2 degAng2 dos
global m b k_R K_R d_R K_1 K_2 gamma3
global k_proO2 k_conO2 O2_value d_O2 
global d_C K_C k_C th_SH th_H th_x th_pT k_E K_E k_cte d_cte
global k_V d_V barkon_v barkon_u barkonx barkoff barkoffx Delta rho gammax gammav gammab nR Kbarx Kbar

O2_value=210; %210 muM=normoxia, reference value for O2

%% Angiopoietin module parameter

%Conversion of the initial concentration of angiopoietin from ng/ml to nM: 50ng/ml=0.179 nM 
%Ang_{1,0}=100 ng/ml=2*0.179;
%Ang_{2,0}=[0-4000]ng/ml=dos*0.179 for different values of "dos" in [0,10,20,40,80].
dos=10;

%Ang1 parameter
delta1=0.005;      
alfa1=0.725*2^4; 
gamma1=0.179;

%Ang2 parameter
delta2=0.01/dos;      
alfa2=0.725*2^2;  
alfam2=2.1875;
kAng2=1.53/dos;                  
gamma2=0.0895*dos;
degAng2=0.069;  

%pTie2 parameter                                              
m=8;                            
k_R=600;                           
b=m/k_R-10^(-4);                    
K_R=min(max([0.5 2*m b]),1);  
d_R=0.000684;                     
K_1=0.67;              
K_2=1.5;               
gamma3=11.7647;

%% Cell module parameters

%Tumor cells parameter
k_C=0.00625;
K_C=3;
th_H=50/O2_value; 
th_SH=8.2/O2_value; 

%ECs parameter
th_x=16;    
th_pT=5.8824; 
k_E=0.0031;  
K_E=10; 

%Oxygen parameter
k_proO2=0.0244; 
k_conO2=0.3735;  
d_O2=0.05; 

%Cytokine parameter
d_C=0.0125*10^(-2);
k_cte=0.7733*10^(-3);
d_cte=0.3866*10^(-3) ;

%% VEGF module parameter

k_V=1.9332;
d_V=1.135;
barkon_u=3.316;   
barkon_v=3.316;   
barkonx=14.924*10^(2); 
barkoff=2.625;
barkoffx=6.125;
Kbar=barkon_v/barkoff;
Kbarx=barkonx/barkoffx;
Delta=0.0025;
rho=50;
gammax=0.07;
gammav=0.05; 
gammab=0.35;
nR=7.79;

