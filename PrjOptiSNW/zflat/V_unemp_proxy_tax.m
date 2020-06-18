function [V_U,exitflag_fsolve]=V_unemp_proxy_tax(welf_checks,TR,xi,b,V_unemp_tax,a2_COVID,options2)
%function [V_U,exitflag_fsolve]=V_unemp_proxy_tax(welf_checks,TR,xi,b,V_ss,cons_unemp,ap_unemp,options2)

%% Solve optimization problem

%global beta agrid pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

V_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
exitflag_fsolve=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % Number of kids

                       % Find value of assets that approximates the value of the welfare checks                   
                       x0=agrid(a)+TR*welf_checks; % Initial guess for a

                       [a_aux,~,exitflag_fsolve(j,a,eta,educ,married,kids)]=fsolve(@(x)find_a_unemp_tax(x,j,a,eta,educ,married,kids,TR,xi,b,welf_checks,a2_COVID),x0,options2);

                       if a_aux<0
                           a_aux=0;
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
                       
                       V_U(j,a,eta,educ,married,kids)=vals(1)*V_unemp_tax(j,inds(1),eta,educ,married,kids)+vals(2)*V_unemp_tax(j,inds(2),eta,educ,married,kids);
   
                   end
               end
           end
       end
   end
   
%  disp(j)
   
end


end