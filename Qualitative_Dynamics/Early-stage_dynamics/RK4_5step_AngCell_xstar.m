function [K] = RK4_5step_AngCell_xstar(par,X,p_V,C,O2)

%Function implementing an RK4 method to solve the ODE system for the
%VEGF module: V

   %Paramenters

   k_V=par(1);
   d_V=par(2);
   dt=par(3);
   
   V=X(1);

   %Function for VEGF production, influenced by oxygen levels 
   h_1V=polyval(p_V,O2);
   
   %Vectorial field
   b_V=k_V*h_1V*C-d_V*V;  
   
   B_vec=b_V;
   
   K=dt.*B_vec;
    
end