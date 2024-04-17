function [K] = RK4_4step_AngCell_xstar(par,X,pTie2)

%Function implementing an RK4 method to solve the ODE system for the
%cellular module: C, E, Cyt, O2

   %Paramenters
   k_C=par(1);
   K_C=par(2);
   K_E=par(3);
   thpT=par(4);
   k_E=par(5);
   th_SH=par(6); 
   th_H=par(7);
   th_x=par(8);
   k_proO2=par(9);
   k_conO2=par(10);
   d_O2=par(11);
   d_C=par(12);
   k_cte=par(13);
   d_cte=par(14);
   xstar=par(15);
   dt=par(16);
     
   E=X(1);
   C=X(2);
   O2=X(3);
   Cyto=X(4);
      
   %Function for endothelial cell production, influenced by xstar and pTie2  
   k_E=k_E*thpT^100./(thpT.^100+pTie2.^100); 
   h_1E=k_E*tanh((xstar-th_x));
   
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
   b_E=h_1E*E*(1-E/K_E);
   b_C=h_1C*C*(1-C/K_C)-d_C*Cyto*C;
   b_O2=k_proO2*E-k_conO2*O2*C-d_O2*O2;
   b_Cyto=k_cte*C/(O2+th_SH/5)-d_cte*Cyto;
   
   B=[b_E b_C b_O2 b_Cyto];
   
   K=dt.*B;
    
end