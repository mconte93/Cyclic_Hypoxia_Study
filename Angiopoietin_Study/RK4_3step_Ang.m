function [K] = RK4_3step_Ang(par,X,Tie2,Ang1,Ang2,z1,z2,Ang2Tie2_2,Ang1Tie2_4)

%Function implementing an RK4 method to solve the ODE system for the
%angiopoietin module: pTie2, R

   %Paramenters
   m=par(1);
   b=par(2);
   k_R=par(3);
   K_R=par(4);
   d_R=par(5);
   K_1=par(6);
   K_2=par(7);
   k2=par(8);
   kminus2=par(9);
   k1=par(10);
   kminus1=par(11);
   Ang1_0=par(12);
   Ang2_0=par(13);
   dt=par(14);
   
   R=X(2);
   
   %Evaluation of the Beware function
   if Ang2==0
      sum=0; 
           for j1=0:4
                j2=0;
               sum=sum+Z_beware(j1,j2,4,Ang1,Ang2,Tie2,K_1,K_2);
           end
      g2=Z_beware(4,0,4,Ang1,Ang2,Tie2,K_1,K_2)/sum; 
  
   elseif Ang1==0
       
   M1=9/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));
   M2=16/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));
   M3=75/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));

      sum=0; 
           for j2=0:1
                j1=0;
               sum=sum+Z_beware(j1,j2,4,Ang1,Ang2,Tie2,K_1,K_2);
           end
       g2=M1/(sum+M1+M2+M3); 
  
   else
       
   M1=9/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));
   M2=16/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));
   M3=75/100*(Z_beware(0,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(1,2,4,Ang1,Ang2,Tie2,K_1,K_2)+Z_beware(2,2,4,Ang1,Ang2,Tie2,K_1,K_2));

   sum=0; 
    for j1=0:4
        for j2=0:1
            if j1+j2<=4
               sum=sum+Z_beware(j1,j2,4,Ang1,Ang2,Tie2,K_1,K_2);
            end
        end
    end

   g2=(Z_beware(4,0,4,Ang1,Ang2,Tie2,K_1,K_2)+M1)/(sum+M1+M2+M3);
   end
  
   %Evaluation of the functions involved in the evolution of R and pTie2
   [Fun,ftilde]=pTie2_interpolation_adim(Ang1_0,Ang2_0,z1,z2,Tie2);
   
   %Evaluation of F(R)
   F=g2*R*(k_R*(1-R/K_R)-m/(R+b))-d_R*R;

   %Vectorial fields  
   V1=k1*Ang1*(Tie2)^4 - kminus1*Ang1Tie2_4;
   V2=k2*Ang2*(Tie2)^2 - kminus2*Ang2Tie2_2; 
   dTie2dt=-V1-V2;
   
   b_pTie2=dTie2dt*(R-F*Fun*ftilde);
   b_R=-F*Fun*(ftilde/Tie2)*dTie2dt;
   
   B=[b_pTie2 b_R];
   
   K=dt.*B;
    
end