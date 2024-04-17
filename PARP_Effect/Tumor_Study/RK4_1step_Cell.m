function [K] = RK4_1step_Cell(par,X,O2)

%Function implementing an RK4 method to solve the reduced ODE system for the
%cellular module: C 

   %Paramenters
   k_C=par(1);
   K_C=par(2);
   th_SH=par(3); 
   th_H=par(4);  
   Pj=par(5);
   nj=par(6);
   dt=par(7);
     
   C=X;
      
%Function for tumor cell production, influenced by oxygen levels

   %Coefficient to ensure continuity of h{1,C}(O2) at th_H
   alfa=0.0314;
    if O2>=th_H
        h_1C=k_C;
    else
        h_1C=tanh(alfa*(O2-th_SH));
    end

%Function accounting for the effect of PARP
    if Pj==0
        f_Pj=1;
    else
        f_Pj=1/((1+Pj)^nj);
    end

   %Vectorial fields 
   b_C=h_1C*f_Pj*C*(1-C/K_C);
   
   B=b_C;
   
   K=dt.*B;
    
end