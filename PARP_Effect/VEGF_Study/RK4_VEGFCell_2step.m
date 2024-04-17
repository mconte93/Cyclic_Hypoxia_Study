function [K] = RK4_VEGFCell_2step(par,X,p_V,h2_Pj,C,O2)

%Function implementing an RK4 method to solve the reduced ODE system for the
%VEGF module to study the effect of PARP on VEGF production: V

   %Paramenters
   k_V=par(1);
   d_V=par(2);
   dt=par(3);
   
   V=X(1);

   %Function for VEGF production, influenced by oxygen levels 
   h1_O2=polyval(p_V,O2);

   %Vectorial field
   b_V=k_V*h1_O2*h2_Pj*C-d_V*V;  
   
   B_vec=b_V;
   
   K=dt.*B_vec;
    
end