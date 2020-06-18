function F=find_a_working(a_aux,j,a,eta,educ,married,kids,TR,welf_checks)

global theta r agrid epsilon eta_grid SS

inc=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));

resource_target=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))+TR*welf_checks;

inc=r*a_aux+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));

resource_aux=(1+r)*a_aux+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc));

F=resource_target-resource_aux;

end


