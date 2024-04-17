function [K] = RK4_VEGFCell_1step(par,X,O2)

%Function implementing an RK4 method to solve the reduced ODE system for the
%cellular module to study the effect of PARP on VEGF production: C, Cyt

   %Paramenters

   k_C=par(1);
   K_C=par(2);
   th_SH=par(3); 
   th_H=par(4);
   d_C=par(5);
   k_cte=par(6);
   d_cte=par(7);
   dt=par(8);
     
   C=X(1);
   Cyto=X(2);
   
   %Function for tumor cell production, influenced by oxygen levels

   %Coefficient to ensure continuity of h{1,C}(O2) at th_H
   alfa=0.0314;
    if O2>=th_H
        h_1C=k_C;
    else
        h_1C=tanh(alfa*(O2-th_SH));
    end

    %Function for cytokines production, influenced by oxygen levels
    if sign(th_H-O2)<0
        k_cte=0;
    end

   %Vectorial fields 
   b_C=h_1C*C*(1-C/K_C)-d_C*Cyto*C;
   b_Cyto=k_cte*(C/(O2+th_SH/5))-d_cte*Cyto;
   
   B=[b_C b_Cyto];
   
   K=dt.*B;
    
end