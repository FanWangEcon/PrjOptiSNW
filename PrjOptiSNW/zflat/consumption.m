function consumption=consumption(j,a,eta,educ,married,kids,aprime,xi,b)

global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option

if nargin==7
	[inc,earn]=individual_income(j,a,eta,educ);
    spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                                                              
    consumption=(1+r)*(agrid(a)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))-max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))))-aprime;
elseif nargin==9
 	[inc,earn]=individual_income(j,a,eta,educ,xi,b);
    spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

    consumption=(1+r)*(agrid(a)+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi))-max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi))))-aprime;
else
    error('Incorrect number of function inputs')
end

end