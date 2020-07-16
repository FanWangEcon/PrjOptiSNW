function [income,earnings]=snw_hh_individual_income(j,a,eta,educ,...
    theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
    xi,b)

if nargin==12
    income=r*(agrid(a)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ);
    earnings=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
elseif nargin==14
    income=r*(agrid(a)+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ);
    earnings=epsilon(j,educ)*theta*exp(eta_H_grid(eta)); % What earnings are before we account for the drop in earnings due to unemployment
else
   error('Incorrect number of function inputs') 
end

end