function [K] = RK4_2step_AngCell_xstar(par,X,Tie2,y11,y12,y21,y22,z1,z2,Ang1,Ang2,step)

%Function implementing an RK4 method to solve the ODE system for the
%angiopoietin module: yij, zii for i,j=1,2

   %Paramenters
   alfa1=par(1);
   alfa2=par(2);
   alfam2=par(3);
   gamma1=par(4);
   gamma2=par(5);
   Ang1_0=par(6);
   dos=par(7);
   dt=par(8);
  
   y11_n=X(1);   %
   y12_n=X(2);
   y21_n=X(3);
   y22_n=X(4);
   z1_n=X(5);
   z2_n=X(6);
   
   %We approximate the integral at the temporal step (m+1) using all the information until the temporal step m
   H=dt;
   m=step-1;
   int11=0;
   int12=0;
   int21=0;
   int22=0;
   
   if m==1
     k=m-1;
     indx=k+1;
     int11_step=H/2*exp((dt*k-dt*(step-1))).*(4.*gamma1.*Tie2(indx).^3.*Ang1(indx).*z1(indx)+Tie2(indx).^4.*y11(indx));
     V_y11=-alfa1*Tie2(step).^4*y11_n-4*alfa1*gamma1*Tie2(step)^3*Ang1(step)*z1_n+alfa1*(int11_step+H/2.*(4.*gamma1*Tie2(step).^3.*Ang1(step).*z1_n+Tie2(step).^4.*y11_n));
 
     int12_step=H/2*exp((dt*k-dt*(step-1))).*(4.*gamma1*Tie2(indx).^3.*Ang1(indx).*z2(indx)+Tie2(indx).^4.*y12(indx));
     V_y12=-alfa1*Tie2(step).^4*y12_n-4*alfa1*gamma1*Tie2(step)^3*Ang1(step)*z2_n+alfa1*(int12_step+H/2.*(4.*gamma1*Tie2(step).^3.*Ang1(step).*z2_n+Tie2(step).^4.*y12_n));

     int21_step=H/2*exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z1(indx)+Tie2(indx).^2.*y21(indx));
     V_y21=-alfa2*Tie2(step).^2*y21_n-2*alfa2*gamma2*Tie2(step)*Ang2(step)*z1_n+alfa2*alfam2*(int21_step+H/2.*(2.*gamma2*Tie2(step).*Ang2(step).*z1_n+Tie2(step).^2.*y21_n));

     int22_step=H/2*exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z2(indx)+Tie2(indx).^2.*y22(indx));
     V_y22=-alfa2*Tie2(step).^2*y22_n-2*alfa2*gamma2*Tie2(step)*Ang2(step)*z2_n+alfa2*alfam2*(int22_step+H/2.*(2.*gamma2*Tie2(step).*Ang2(step).*z2_n+Tie2(step).^2.*y22_n));
              
   else
       
       for k=0:m-2
       indx=k+1;
       int11_1=exp((dt*k-dt*(step-1))).*(4.*gamma1*Tie2(indx).^3.*Ang1(indx).*z1(indx)+Tie2(indx).^4.*y11(indx));
       int11_2=exp((dt*(k+1)-dt*(step-1))).*(4.*gamma1*Tie2(indx+1).^3.*Ang1(indx+1).*z1(indx+1)+Tie2(indx+1).^4.*y11(indx+1));
       int11=int11+H/2*(int11_1+int11_2);
       
       int12_1=exp((dt*k-dt*(step-1))).*(4.*gamma1*Tie2(indx).^3.*Ang1(indx).*z2(indx)+Tie2(indx).^4.*y12(indx));
       int12_2=exp((dt*(k+1)-dt*(step-1))).*(4.*gamma1*Tie2(indx+1).^3.*Ang1(indx+1).*z2(indx+1)+Tie2(indx+1).^4.*y12(indx+1));
       int12=int12+H/2*(int12_1+int12_2);
       
       int21_1=exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z1(indx)+Tie2(indx).^2.*y21(indx));
       int21_2=exp(alfam2*(dt*(k+1)-dt*(step-1))).*(2.*gamma2*Tie2(indx+1).*Ang2(indx+1).*z1(indx+1)+Tie2(indx+1).^2.*y21(indx+1));
       int21=int21+H/2*(int21_1+int21_2);
       
       int22_1=exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z2(indx)+Tie2(indx).^2.*y22(indx));
       int22_2=exp(alfam2*(dt*(k+1)-dt*(step-1))).*(2.*gamma2*Tie2(indx+1).*Ang2(indx+1).*z2(indx+1)+Tie2(indx+1).^2.*y22(indx+1));
       int22=int22+H/2*(int22_1+int22_2);   
       end
  
   k=m-1;
   indx=k+1;
   int11_step=int11+H/2*exp((dt*k-dt*(step-1))).*(4.*gamma1*Tie2(indx).^3.*Ang1(indx).*z1(indx)+Tie2(indx).^4.*y11(indx));
   V_y11=-alfa1*Tie2(step).^4*y11_n-4*alfa1*gamma1*Tie2(step)^3*Ang1(step)*z1_n+alfa1*(int11_step+H/2.*(4.*gamma1*Tie2(step).^3.*Ang1(step).*z1_n+Tie2(step).^4.*y11_n));

   int12_step=int12+H/2*exp((dt*k-dt*(step-1))).*(4.*gamma1*Tie2(indx).^3.*Ang1(indx).*z2(indx)+Tie2(indx).^4.*y12(indx));
   V_y12=-alfa1*Tie2(step).^4*y12_n-4*alfa1*gamma1*Tie2(step)^3*Ang1(step)*z2_n+alfa1*(int12_step+H/2.*(4.*gamma1*Tie2(step).^3.*Ang1(step).*z2_n+Tie2(step).^4.*y12_n));

   int21_step=int21+H/2*exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z1(indx)+Tie2(indx).^2.*y21(indx));
   V_y21=-alfa2*Tie2(step).^2*y21_n-2*alfa2*gamma2*Tie2(step)*Ang2(step)*z1_n+alfa2*alfam2*(int21_step+H/2.*(2.*gamma2*Tie2(step).*Ang2(step).*z1_n+Tie2(step).^2.*y21_n));

   int22_step=int22+H/2*exp(alfam2*(dt*k-dt*(step-1))).*(2.*gamma2*Tie2(indx).*Ang2(indx).*z2(indx)+Tie2(indx).^2.*y22(indx));
   V_y22=-alfa2*Tie2(step).^2*y22_n-2*alfa2*gamma2*Tie2(step)*Ang2(step)*z2_n+alfa2*alfam2*(int22_step+H/2.*(2.*gamma2*Tie2(step).*Ang2(step).*z2_n+Tie2(step).^2.*y22_n));
   end
   
   %Vectorial fields  
    if Ang1_0==0
        b_y11=0;
        b_y21=0;
        b_y12=V_y12;
        b_y22=V_y22;
    elseif dos==0
        b_y11=V_y11;
        b_y21=V_y21;
        b_y12=0;
        b_y22=0;
    else
    b_y11=V_y11;
    b_y12=V_y12;
    b_y21=V_y21;
    b_y22=V_y22;
    end
   
    b_z1=b_y11+b_y21;
    b_z2=b_y12+b_y22;
    
   B=[b_y11 b_y12 b_y21 b_y22 b_z1 b_z2];
   
   K=dt.*B;
    
end