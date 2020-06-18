function consumption=consumption(j,a,eta,educ,married,kids,aprime,xi,b)

global theta r agrid epsilon eta_grid SS

if nargin==7
    [inc,earn]=individual_income(j,a,eta,educ);
    spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

    consumption=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))-aprime;
elseif nargin==9
    [inc,earn]=individual_income(j,a,eta,educ,xi,b);
    spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

    consumption=(1+r)*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))-aprime;
else
    error('Inccorect number of function inputs')
end

end