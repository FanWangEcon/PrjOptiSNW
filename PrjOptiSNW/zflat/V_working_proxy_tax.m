function [V_W,C_W,exitflag_fsolve]=V_working_proxy_tax(welf_checks,TR,V_working_tax,cons_working_tax,a2_COVID,options2)

%% Solve optimization problem

global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
C_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
exitflag_fsolve=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % Number of kids

                       % Find value of assets that approximates the value of the welfare checks                   
                       x0=agrid(a)+TR*welf_checks; % Initial guess for a

                       [a_aux,~,exitflag_fsolve(j,a,eta,educ,married,kids)]=fsolve(@(x)find_a_working_tax(x,j,a,eta,educ,married,kids,TR,welf_checks,a2_COVID),x0,options2);

                       if a_aux<0
						   disp(a_aux)
                           error('Check code! Should not allow for negative welfare checks')
                       elseif a_aux>agrid(n_agrid)
                           a_aux=agrid(n_agrid);
                       end

                       % Linear interpolation
                       ind_aux=find(agrid<=a_aux,1,'last');
                       
                       if a_aux==0
                           inds(1)=1;
                           inds(2)=1;                       
                           vals(1)=1;
                           vals(2)=0;

                       elseif a_aux==agrid(n_agrid)
                           inds(1)=n_agrid;
                           inds(2)=n_agrid;                       
                           vals(1)=1;
                           vals(2)=0;

                       else
                           inds(1)=ind_aux;
                           inds(2)=ind_aux+1;                       
                           vals(1)=1-((a_aux-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                           vals(2)=1-vals(1);

                       end
                       
                       V_W(j,a,eta,educ,married,kids)=vals(1)*V_working_tax(j,inds(1),eta,educ,married,kids)+vals(2)*V_working_tax(j,inds(2),eta,educ,married,kids);
%                       C_W(j,a,eta,educ,married,kids)=vals(1)*cons_working_tax(j,inds(1),eta,educ,married,kids)+vals(2)*cons_working_tax(j,inds(2),eta,educ,married,kids);
                       
                       C_W(j,a,eta,educ,married,kids)=vals(1)*(cons_working_tax(j,inds(1),eta,educ,married,kids)/(married+kids-1))+vals(2)*(cons_working_tax(j,inds(2),eta,educ,married,kids)/(married+kids-1));
   
                   end
               end
           end
       end
   end
   
%  disp(j)
   
end


end