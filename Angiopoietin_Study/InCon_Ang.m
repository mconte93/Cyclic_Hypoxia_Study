%Adimensional initial conditions
global dos

%Angiopoietin module   
Ang1_0=0.179*2;           
Ang1Tie2_4_0=0.179*0.001;     
Ang2_0=0.179*dos; 
if dos==0
Ang2Tie2_2_0=0;
else
    Ang2Tie2_2_0=0.179*0.001;    
end
Tie2_0=2;           
y11_0=-0.001;
y12_0=-0.001;
y21_0=-0.001;
y22_0=-0.001;
z1_0=-0.001;
z2_0=-0.001;
pTie2_0=0.17;          
R_0=pTie2_0/Tie2_0;       
R_contr=pTie2_0/Tie2_0;

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