function F=value_func_aux_working_tax(ap_aux,j,a,eta,educ,married,kids,V_ss,a2_COVID)

global beta theta r agrid epsilon eta_grid SS pi_eta pi_kids psi n_etagrid n_kidsgrid

% inc=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
% spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));

[inc,earn]=individual_income(j,a,eta,educ);
spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

c_aux=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax_COVID(inc,(married-1)*spouse_inc,a2_COVID))-ap_aux;

ind_aux=find(agrid<=ap_aux,1,'last');

% Linear interpolation
vals(1)=1-((ap_aux-agrid(ind_aux))/(agrid(ind_aux+1)-agrid(ind_aux)));
vals(2)=1-vals(1);

cont=0;
for etap=1:n_etagrid
    for kidsp=1:n_kidsgrid
        cont=cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married)*(vals(1)*V_ss(j+1,ind_aux,etap,educ,married,kidsp)+vals(2)*V_ss(j+1,ind_aux+1,etap,educ,married,kidsp));
    end
end

F=-(utility(c_aux,married,kids)+beta*psi(j)*cont);

end