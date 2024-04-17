function [K] = RK4_1step_Ang(par,X)

%Function implementing an RK4 method to solve the ODE system for the
%angiopoietin module: Ang1, Ang2, (Ang1.Tie2)4, (Ang2.Tie2)2, Tie2

   %Paramenters
   k1=par(1); 
   kminus1=par(2);
   k2=par(3);
   kminus2=par(4);
   dt=par(5);

   Ang1=X(1);
   Ang1Tie2_4=X(2);
   Ang2=X(3);
   Ang2Tie2_2=X(4);
   Tie2=X(5);
   
   %Vectorial fields
   V1=k1*Ang1*(Tie2)^4 - kminus1*Ang1Tie2_4; 
   V2=k2*Ang2*(Tie2)^2 - kminus2*Ang2Tie2_2;

   b_Ang1=-V1;
   b_Ang1Tie2_4= V1; 
   b_Ang2=-V2;
   b_Ang2Tie2_2=V2;
   b_Tie2=-V1-V2;
   
   B=[b_Ang1 b_Ang1Tie2_4 b_Ang2 b_Ang2Tie2_2 b_Tie2];
   
   K=dt.*B;
    
end