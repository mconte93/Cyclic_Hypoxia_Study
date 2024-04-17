function [Fun,ftilde]=pTie2_interpolation_adim(z1,z2,Tie2,dos,Ang1_0,Ang2_0)

%Function for the evaluation of the coefficient functions involved in the
%evolution of R and pTie2, starting from the data in Bogdanovic (2016)

%Data extract and interpolated
load('Int_figBogdanovic_pchip1_pchip2_098_06_032.mat','f1','f2','df1','df2')

%Initial values for Tie2, Ang1, and Ang2
Tie2_0=2;
Ang10=Ang1_0*0.179*2;
Ang20=Ang2_0*0.179*dos;

%Scaling parameter for df2/dAng20, estimated from the Ph.D. thesis of Alawo
%(2017)
sc_2=0.25; 

%Evaluation of ftilde, which is an weighted sum of the data about Ang1
%and Ang2 evolution, and the weighted sum of the derivative of f w.r.t. Ang10 and
%Ang20 (Fun)

x1=Ang10;
x2=Ang20;
d1=1/(x2);  %Inverse of the distance from the Ang1_0 axis
d2=1/(x1);  %Inverse of the distance from the Ang2_0 axis

if x1==0
    f_2=sc_2*ppval(f2,x2);
    ftilde=f_2;
    f=ftilde/(Tie2*Tie2_0);
    Fun=1/(sc_2*ppval(df2,x2)-z2*f);

elseif x2==0
    
    f_1=ppval(f1,x1);
    ftilde=f_1;
    f=ftilde/(Tie2*Tie2_0);
    Fun=1/(ppval(df1,x1)-z1*f);
    
else
       r1=x1/(x1+x2);
       r2=x2/(x1+x2);
       f_1=ppval(f1,x1);
       f_2=sc_2*ppval(f2,x2);
       
       ftilde=r1*f_1+r2*f_2;
       f=ftilde/(Tie2*Tie2_0);
       
       Fun=1/(d2*r2*(f_1-f_2)+ppval(df1,x1)-1/r1*z1*f+d1*r1*(f_2-f_1)+sc_2*ppval(df2,x2)-1/r2*z2*f);

end