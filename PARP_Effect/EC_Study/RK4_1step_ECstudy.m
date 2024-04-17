function [K] = RK4_1step_ECstudy(par,X)

%Function implementing an RK4 method to solve the reduced ODE system for the
%cellular module: E

   %Paramenters
   k_E=par(1);
   K_E=par(2);
   Pj=par(3); 
   nPj=par(4);
   dt=par(5);
     
   E=X(1);
      
   %Vectorial fields 
   b_E=k_E/(1+Pj)^nPj*E*(1-E/K_E);
   
   B=b_E;
   
   K=dt.*B;
    
end