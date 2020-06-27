function [V_VFI,ap_VFI,cons_VFI,exitflag_VFI]=VFI_grid_search

%% Solve optimization problem

global beta theta r agrid epsilon eta_grid SS pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

exitflag_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Solve for value function and policy functions by means of backwards induction
for j=n_jgrid:(-1):1 % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % Number of kids

                       if j==n_jgrid

                           inc=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
                           spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));
                           
                           ap_VFI(j,a,eta,educ,married,kids)=1;
                           cons_VFI(j,a,eta,educ,married,kids)=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc));

                           if cons_VFI(j,a,eta,educ,married,kids)<=0
                              disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                              error('Non-positive consumption')
                           end

                           V_VFI(j,a,eta,educ,married,kids)=utility(cons_VFI(j,a,eta,educ,married,kids),married,kids);

                       else
                           
                           inc=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
                           spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));
                           
                           V_VFI_aux=NaN(n_agrid,1);
                           
                           for aa=1:n_agrid
                               
                               cont=0;
                               for etap=1:n_etagrid
                                   for kidsp=1:n_kidsgrid
                                       cont=cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)*V_VFI(j+1,aa,etap,educ,married,kidsp);
                                   end
                               end

                               c_aux=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))-agrid(aa);
                               
                               V_VFI_aux(aa)=utility(c_aux,married,kids)+beta*psi(j)*cont;
                               
                           end
                           
                           [max_val,max_ind]=max(V_VFI_aux);
                           
                           ap_VFI(j,a,eta,educ,married,kids)=max_ind;
                           cons_VFI(j,a,eta,educ,married,kids)=(1+r)*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))-agrid(max_ind);

                           V_VFI(j,a,eta,educ,married,kids)=max_val;
                           
                           if cons_VFI(j,a,eta,educ,married,kids)<=0
                              disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                              error('Non-positive consumption')
                           end

                       end

                   end
               end
           end
       end
   end
   
   disp(j)
   
end


end