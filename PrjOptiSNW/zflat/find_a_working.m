function F=find_a_working(a_aux,j,a,eta,educ,married,kids,TR,welf_checks)

global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option

[inc,earn]=individual_income(j,a,eta,educ);
spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

resource_target=(1+r)*(agrid(a)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))-max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))))+TR*welf_checks;

inc=r*(a_aux+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ);
earn=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

resource_aux=(1+r)*(a_aux+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta))+SS(j,educ)+(married-1)*spouse_inc*exp(eta_S_grid(eta))-max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))));

F=resource_target-resource_aux;

end


