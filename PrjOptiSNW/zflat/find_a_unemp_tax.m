function F=find_a_unemp_tax(a_aux,j,a,eta,educ,married,kids,TR,xi,b,welf_checks,a2_COVID)

global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option

[inc,earn]=individual_income(j,a,eta,educ,xi,b);
spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

resource_target=(1+r)*(agrid(a)+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi))-max(0,Tax_COVID(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi)),a2_COVID))+TR*welf_checks;

inc=r*(a_aux+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ);
earn=epsilon(j,educ)*theta*exp(eta_H_grid(eta)); % What earnings are before we account for the drop in earnings due to unemployment
spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

resource_aux=(1+r)*(a_aux+Bequests*(bequests_option-1))+( epsilon(j,educ)*theta*exp(eta_H_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi))-max(0,Tax_COVID(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi)),a2_COVID));

F=resource_target-resource_aux;

end