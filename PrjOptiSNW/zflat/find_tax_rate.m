function F=find_tax_rate(a2_guess,Phi_true,omega,xi,b,cutoffs)

global theta r g_cons agrid epsilon eta_grid SS pi_unemp n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

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

                       spouse_inc=spousal_income(j,educ,kids,( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi)),SS(j,educ));
                       Y_inc_agg=Y_inc_agg+Phi_true(j,a,eta,educ,married,kids)*( r*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc ); % Aggregate income
                                              
                       SS_spend=SS_spend+Phi_true(j,a,eta,educ,married,kids)*SS(j,educ); % Total spending on Social Security
                       
                       wages=epsilon(j,educ)*theta*exp(eta_grid(eta));
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
                                              
                       UI_benefits=UI_benefits+Phi_true(j,a,eta,educ,married,kids)*pi_unemp(j,wage_ind)*( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*b*(1-xi); % Total spending on unemployment insurance benefits
                       
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
                           
                           inc_aux=r*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ); % Income excluding spousal income (if married)
                           spouse_inc=spousal_income(j,educ,kids,( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi)),SS(j,educ));
                           Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc_aux,(married-1)*spouse_inc)); % Tax revenues

                       end                       
                   end
               end
           end
       end
    end
    
    a2_guess=a2_guess*(((UI_benefits+SS_spend+(g_cons+omega)*Y_inc_agg)/Tax_revenues_aux)^0.5); % Find value of a2 that balances government budget
    
    err=abs((Tax_revenues_aux/(UI_benefits+SS_spend+(g_cons+omega)*Y_inc_agg))-1);
    
    disp(err)
    
end

a2_update=[a2_guess_orig,a2_guess];

name='Old and updated value of a2=';
name2=[name,num2str(a2_update)];
disp(name2);

F=a2_guess;

end