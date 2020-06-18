function [income,earnings]=individual_income(j,a,eta,educ,xi,b)

global theta r agrid epsilon eta_grid SS

if nargin==4
    income=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
    earnings=epsilon(j,educ)*theta*exp(eta_grid(eta));
elseif nargin==6
    income=r*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ);
    earnings=( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi));
else
   error('Inccorect number of function inputs') 
end

end