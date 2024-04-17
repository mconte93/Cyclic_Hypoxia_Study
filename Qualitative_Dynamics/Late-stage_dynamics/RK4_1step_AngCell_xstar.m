function [K] = RK4_1step_AngCell_xstar(par,X,E,V)

%Function implementing an RK4 method to solve the ODE system for the
%angiopoietin module: Ang1, Ang2, (Ang1.Tie2)4, (Ang2.Tie2)2, Tie2

   %Paramenters
   delta1=par(1); 
   alfa1=par(2);
   delta2=par(3);
   alfa2=par(4);
   alfam2=par(5);
   gamma1=par(6);
   gamma2=par(7);
   kAng2=par(8);
   dAng2=par(9);
   dt=par(10);

   Ang1=X(1);
   Ang1Tie2_4=X(2);
   Ang2=X(3);
   Ang2Tie2_2=X(4);
   Tie2=X(5);
   
   %Function for Ang2 production by VEGF
   h_Ang2=kAng2*V*tanh(V);
   
   %Vectorial fields  
   b_Ang1=-alfa1*Ang1*(Tie2)^4+delta1*Ang1Tie2_4;
   b_Ang1Tie2_4=-Ang1Tie2_4+alfa1/delta1*Ang1*(Tie2)^4; 
   b_Ang2=delta2*alfam2*Ang2Tie2_2-alfa2*Ang2*(Tie2)^2+h_Ang2*E-dAng2*Ang2;
   b_Ang2Tie2_2=-alfam2*Ang2Tie2_2+alfa2/delta2*Ang2*(Tie2)^2;
   b_Tie2=-alfa1*gamma1*Ang1*(Tie2)^4+delta1*gamma1*Ang1Tie2_4+alfam2*delta2*gamma2*Ang2Tie2_2-alfa2*gamma2*Ang2*(Tie2)^2;
   
   B=[b_Ang1 b_Ang1Tie2_4 b_Ang2 b_Ang2Tie2_2 b_Tie2];
   
   K=dt.*B;
    
end