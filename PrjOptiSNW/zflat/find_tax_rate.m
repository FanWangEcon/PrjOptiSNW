function F=find_tax_rate(a2_guess,Phi_true,Gov_cons,Covid_checks,xi,b,cutoffs)

global theta g_cons epsilon agrid eta_H_grid eta_S_grid SS r a2 pi_unemp n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid Bequests bequests_option

disp('Assuming spouse (if married) is also getting UI benefits')

% Aggregate variables
Y_inc_agg=0;
SS_spend=0;
UI_benefits=0;

for j=1:n_jgrid
   for a=1:n_agrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid
                       
                       wages=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
                       if wages<=cutoffs(1)
                           wage_ind=1;
                       elseif wages>cutoffs(1) && wages<=cutoffs(2)
                           wage_ind=2;
                       elseif wages>cutoffs(2) && wages<=cutoffs(3)
                           wage_ind=3;
                       elseif wages>cutoffs(3) && wages<=cutoffs(4)
                           wage_ind=4;
                       elseif wages>cutoffs(4)
                           wage_ind=5;
                       end
                       
                       [~,earn]=individual_income(j,a,eta,educ); % What individual earnings are before we account for the drop in earnings due to unemployment
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ)); % What average spousal earnings are before we account for the drop in earnings due to unemployment
                       
					   inc_aux=pi_unemp(j,wage_ind)*epsilon(j,educ)*theta*exp(eta_H_grid(eta))*xi+(1-pi_unemp(j,wage_ind))*epsilon(j,educ)*theta*exp(eta_H_grid(eta))+r*(agrid(a)+Bequests*(bequests_option-1)); % Income (excluding Social Security benefits) after accounting for potential earnings drop in case of unemployment
                       
		               Y_inc_agg=Y_inc_agg+Phi_true(j,a,eta,educ,married,kids)*( inc_aux+pi_unemp(j,wage_ind)*(married-1)*spouse_inc*exp(eta_S_grid(eta))*xi+(1-pi_unemp(j,wage_ind))*(married-1)*spouse_inc*exp(eta_S_grid(eta)) ); % Aggregate income (labor earnings, spousal income, interest earnings)
					                                                 
                       SS_spend=SS_spend+sum(Phi_true(j,a,eta,educ,married,kids)*SS(j,educ)); % Total spending on Social Security
                                              
                       UI_benefits=UI_benefits+sum(Phi_true(j,a,eta,educ,married,kids))*pi_unemp(j,wage_ind)*( epsilon(j,educ)*theta*exp(eta_H_grid(eta))+(married-1)*spouse_inc*exp(eta_S_grid(eta)) )*b*(1-xi); % Total spending on unemployment insurance benefits
					   
                   end                       
               end
           end
       end
   end
end
						   


% Find value of a2 that balances government budget
tol=10^-4;
err=1;

a2_guess_orig=a2_guess;

while err>tol
    
    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
       for a=1:n_agrid
           for eta=1:n_etagrid
               for educ=1:n_educgrid
                   for married=1:n_marriedgrid
                       for kids=1:n_kidsgrid
                           
                           [inc,earn]=individual_income(j,a,eta,educ,xi,b);
                           spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                           
                           Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi)))); % Tax revenues

                       end                       
                   end
               end
           end
       end
    end
    
    a2=a2*(((UI_benefits+SS_spend+Gov_cons+Covid_checks)/Tax_revenues_aux)^0.75); % Find value of a2 that balances government budget
    
    err=abs((Tax_revenues_aux/(UI_benefits+SS_spend+(g_cons+omega)*Y_inc_agg))-1);
    
    disp(err)
    
end

a2_update=[a2_guess_orig,a2];

name='Old and updated value of a2=';
name2=[name,num2str(a2_update)];
disp(name2);

F=a2;

end