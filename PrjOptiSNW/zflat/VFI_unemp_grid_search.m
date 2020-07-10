function [V_VFI,ap_VFI,cons_VFI,exitflag_VFI]=VFI_unemp_grid_search(V_ss,xi,b)

%% Solve optimization problem

global beta pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid	

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

exitflag_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % Number of kids
                        
                       if j<n_jgrid
                           
                           c_aux=consumption_grid_search(j,a,eta,educ,married,kids,1:n_agrid,xi,b);
								   
		                   cont=zeros(n_agrid,n_etagrid);
		                   for kidsp=1:n_kidsgrid
							   cont(1:n_agrid,1:n_etagrid)=cont(1:n_agrid,1:n_etagrid)+pi_kids(kids,kidsp,j,educ,married)*reshape(V_ss(j+1,1:n_agrid,1:n_etagrid,educ,married,kidsp),n_agrid,n_etagrid);
		                   end

		                   cont=pi_eta(eta,:)*cont(1:n_agrid,:)';	   

						   V_VFI_aux=utility_grid_search(c_aux,married,kids)+beta*psi(j)*cont';
							  
						   [max_val,max_ind]=max(V_VFI_aux);

						   ap_VFI(j,a,eta,educ,married,kids)=max_ind;
						   cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids),xi,b);

						   V_VFI(j,a,eta,educ,married,kids)=max_val;

						   if cons_VFI(j,a,eta,educ,married,kids)<=0
							   disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
   							   error('Non-positive consumption')
						   end
						   
                       elseif j==n_jgrid
                                                      
                           ap_VFI(j,a,eta,educ,married,kids)=1;
                           cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids),xi,b);

                           if cons_VFI(j,a,eta,educ,married,kids)<=0
                              disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                              error('Non-positive consumption')
                           end

                           V_VFI(j,a,eta,educ,married,kids)=utility_grid_search(cons_VFI(j,a,eta,educ,married,kids),married,kids);
                           
                       end

                   end
               end
           end
       end
   end
   
   disp(j)
   
end


end