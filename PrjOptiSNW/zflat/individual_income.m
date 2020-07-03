function [income,earnings]=individual_income(j,a,eta,educ,xi,b)

global theta r agrid epsilon eta_H_grid SS Bequests bequests_option

if nargin==4
    income=r*(agrid(a)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ);
    earnings=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
elseif nargin==6
    income=r*(agrid(a)+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ);
    earnings=( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi));
else
   error('Inccorect number of function inputs') 
end

end