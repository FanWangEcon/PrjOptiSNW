function [V_W,exitflag_fsolve]=V_working_proxy_tax(welf_checks,TR,V_tax,a2_COVID,options2)
%function [V_W,exitflag_fsolve]=V_working_proxy_tax(welf_checks,TR,V_ss,cons_ss,ap_ss,a2_COVID,options2)

%% Solve optimization problem

%global beta agrid pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
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
                       
                       V_W(j,a,eta,educ,married,kids)=vals(1)*V_tax(j,inds(1),eta,educ,married,kids)+vals(2)*V_tax(j,inds(2),eta,educ,married,kids);
%                        flow=vals(1)*utility(cons_ss(j,inds(1),eta,educ,married,kids),married,kids)+vals(2)*utility(cons_ss(j,inds(2),eta,educ,married,kids),married,kids);
%                        
%                        % Linear interpolation
%                        ind_aux=find(agrid<=ap_ss(j,inds(1),eta,educ,married,kids),1,'last');
%                        
%                        if ap_ss(j,inds(1),eta,educ,married,kids)==0
%                            indsp1(1)=1;
%                            indsp1(2)=1;                       
%                            valsp1(1)=1;
%                            valsp1(2)=0;
% 
%                        elseif ap_ss(j,inds(1),eta,educ,married,kids)>=agrid(n_agrid)
%                            indsp1(1)=n_agrid;
%                            indsp1(2)=n_agrid;                       
%                            valsp1(1)=1;
%                            valsp1(2)=0;
% 
%                        else
%                            indsp1(1)=ind_aux;
%                            indsp1(2)=ind_aux+1;                       
%                            valsp1(1)=1-((ap_ss(j,inds(1),eta,educ,married,kids)-agrid(indsp1(1)))/(agrid(indsp1(2))-agrid(indsp1(1))));
%                            valsp1(2)=1-valsp1(1);
% 
%                        end
%                        
%                        % Linear interpolation
%                        ind_aux=find(agrid<=ap_ss(j,inds(2),eta,educ,married,kids),1,'last');
%                                               
%                        if ap_ss(j,inds(2),eta,educ,married,kids)==0
%                            indsp2(1)=1;
%                            indsp2(2)=1;                       
%                            valsp2(1)=1;
%                            valsp2(2)=0;
% 
%                        elseif ap_ss(j,inds(2),eta,educ,married,kids)>=agrid(n_agrid)
%                            indsp2(1)=n_agrid;
%                            indsp2(2)=n_agrid;                       
%                            valsp2(1)=1;
%                            valsp2(2)=0;
% 
%                        else
%                            indsp2(1)=ind_aux;
%                            indsp2(2)=ind_aux+1;                       
%                            valsp2(1)=1-((ap_ss(j,inds(2),eta,educ,married,kids)-agrid(indsp2(1)))/(agrid(indsp2(2))-agrid(indsp2(1))));
%                            valsp2(2)=1-valsp2(1);
% 
%                        end
%                        
%                        if j<n_jgrid
%                        
%                            cont1=0;
%                            for etap=1:n_etagrid
%                                for kidsp=1:n_kidsgrid
%                                    cont1=cont1+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married)*(valsp1(1)*V_ss(j+1,indsp1(1),etap,educ,married,kidsp)+valsp1(2)*V_ss(j+1,indsp1(2),etap,educ,married,kidsp));
%                                end
%                            end
% 
%                            cont2=0;
%                            for etap=1:n_etagrid
%                                for kidsp=1:n_kidsgrid
%                                    cont2=cont2+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married)*(valsp2(1)*V_ss(j+1,indsp2(1),etap,educ,married,kidsp)+valsp2(2)*V_ss(j+1,indsp2(2),etap,educ,married,kidsp));
%                                end
%                            end
%                            
%                        else
%                            
%                            cont1=0;
%                            cont2=0;
%                            
%                        end
%                            
% 
%                        V_W(j,a,eta,educ,married,kids)=flow+beta*psi(j)*(vals(1)*cont1+vals(2)*cont2);
   
                   end
               end
           end
       end
   end
   
%  disp(j)
   
end


end