function [Phi_true,Phi_adj,A_agg,Y_inc_agg]=Aggregation(ap_ss,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids)

%% Aggregation

global a2 g_cons agrid SS pi_eta pi_kids Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

Phiss=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Assume everyone starts with 0 assets
% Use stationary distribution of productivity shock and initial
% distribution for educational attainment, marital status, and number of
% kids from PSID
for eta=1:n_etagrid % Productivity
   for educ=1:n_educgrid % Fixed effects
       for married=1:n_marriedgrid % Marital status
           for kids=1:n_kidsgrid % No. of kids
               Phiss(1,1,eta,educ,married,kids)=stat_distr_eta(eta)*stat_distr_educ(educ)*stat_distr_married(married)*stat_distr_kids(educ,married,kids);
           end
       end
   end
end

% Use policy functions and survival probabilities to get distribution of remaining idiosyncratic states
for j=1:(n_jgrid-1) % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids

                       if ap_ss(j,a,eta,educ,married,kids)==0
                           inds(1)=1;
                           inds(2)=1;
                           vals(1)=1;
                           vals(2)=0;

                       elseif ap_ss(j,a,eta,educ,married,kids)==agrid(n_agrid)
                           inds(1)=n_agrid;
                           inds(2)=n_agrid;
                           vals(1)=1;
                           vals(2)=0;

                       else

                           ind_aux=find(agrid<=ap_ss(j,a,eta,educ,married,kids),1,'last');

                           inds(1)=ind_aux;
                           inds(2)=ind_aux+1;

                           % Linear interpolation
                           vals(1)=1-((ap_ss(j,a,eta,educ,married,kids)-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                           vals(2)=1-vals(1);

                       end

                       for etap=1:n_etagrid
                           for kidsp=1:n_kidsgrid
                               Phiss(j+1,inds(1),etap,educ,married,kidsp)=Phiss(j+1,inds(1),etap,educ,married,kidsp)+Phiss(j,a,eta,educ,married,kids)*vals(1)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married);
                               Phiss(j+1,inds(2),etap,educ,married,kidsp)=Phiss(j+1,inds(2),etap,educ,married,kidsp)+Phiss(j,a,eta,educ,married,kids)*vals(2)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married);
                           end
                       end

                   end
               end
           end
       end
   end 
end

% Normalize distribution of idiosyncratic states to sum to 1 for each age
Phi_adj=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
Phi_true=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid
    
    dummy=sum(sum(sum(sum(sum(Phiss(j,:,:,:,:,:))))));

    for a=1:n_agrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid

                       if dummy>0
                           Phi_adj(j,a,eta,educ,married,kids)=Phiss(j,a,eta,educ,married,kids)/dummy;
                       else
                           Phi_adj(j,a,eta,educ,married,kids)=0;
                       end
                       
                       Phi_true(j,a,eta,educ,married,kids)=Phi_adj(j,a,eta,educ,married,kids)*Pop(j);

                   end
               end
           end
       end
    end
    
end


% Check if the upper bound on assets binds
check_asset_distr=sum(sum(sum(sum(sum(Phi_true(:,n_agrid,:,:,:,:))))));

name='Share of population with assets equal to upper bound on asset grid=';
name2=[name,num2str(check_asset_distr/sum(Pop))];
disp(name2);

% Aggregate variables
A_agg=0;
Y_inc_agg=0;
Tax_revenues=0;
SS_spend=0;

for j=1:n_jgrid
   for a=1:n_agrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid

                       A_agg=A_agg+Phi_true(j,a,eta,educ,married,kids)*agrid(a); % Aggregate wealth
                       
                       [inc_aux,earn]=individual_income(j,a,eta,educ);
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                       
                       Y_inc_agg=Y_inc_agg+Phi_true(j,a,eta,educ,married,kids)*( inc_aux+(married-1)*spouse_inc ); % Aggregate income
                                              
%                      inc_aux=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ); % Income excluding spousal income (if married)
                       Tax_revenues=Tax_revenues+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc_aux,(married-1)*spouse_inc)); % Tax revenues

                       SS_spend=SS_spend+Phi_true(j,a,eta,educ,married,kids)*SS(j,educ); % Total spending on Social Security

                   end                       
               end
           end
       end
   end
end

% Update guess for a2 (determines average level of income taxation)
% Assuming government balances its budget period-by-period

tol=10^-4;
err=abs((Tax_revenues/(SS_spend+g_cons*Y_inc_agg))-1);

a2_update=a2;

while err>tol
    
    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
       for a=1:n_agrid
           for eta=1:n_etagrid
               for educ=1:n_educgrid
                   for married=1:n_marriedgrid
                       for kids=1:n_kidsgrid
                           
                           [inc_aux,earn]=individual_income(j,a,eta,educ);
                           spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                           Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc_aux,(married-1)*spouse_inc)); % Tax revenues

                       end                       
                   end
               end
           end
       end
    end
    
    a2=a2*(((SS_spend+g_cons*Y_inc_agg)/Tax_revenues_aux)^0.5); % Find value of a2 that balances government budget
    
    err=abs((Tax_revenues_aux/(SS_spend+g_cons*Y_inc_agg))-1);
    
    disp(err)
    
end

a2_update=[a2_update,a2];

name='Old and updated value of a2=';
name2=[name,num2str(a2_update)];
disp(name2);

end