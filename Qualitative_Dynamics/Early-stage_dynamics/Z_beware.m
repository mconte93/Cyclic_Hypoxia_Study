function [Zn] = Z_beware(j1,j2,n,Ang1,Ang2,Tie2,K1,K2)

%Function that calculate the value of Z^n(j1,j2) for the beware coefficient

if Ang2==0
    j0=n-j1;
    Zn = factorial(n)/(factorial(j1)*factorial(j0))*Ang1*(Tie2/K1)^j1;
elseif Ang1==0
    j0=n-j2;
    Zn = factorial(n)/(factorial(j2)*factorial(j0))*Ang2*(Tie2/K2)^j2;
else
    j0=n-j1-j2;
    Zn = factorial(n)/(factorial(j1)*factorial(j2)*factorial(j0))*Ang1*Ang2*(Tie2/K1)^j1*(Tie2/K2)^j2;
end

end
