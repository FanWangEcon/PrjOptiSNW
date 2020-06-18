function [V_planner,Phi_mass]=Planner(Phi_true,j,married,kids,welf_checks,ymin,ymax,V_U,V_W,ap_ss,cutoffs)

%% Value function of planner

global theta r agrid epsilon eta_grid SS pi_eta pi_kids pi_unemp n_agrid n_etagrid n_educgrid n_kidsgrid

Phi_norm=zeros(n_agrid,n_etagrid,n_educgrid);

% Find mass of individuals by group (age, marital status, number of
% kids, and income range) by summing over other idiosyncratic states 
for a=1:n_agrid % Assets
   for eta=1:n_etagrid % Productivity
       for educ=1:n_educgrid % Educational level
           
           inc=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ);
           spouse_inc=spousal_income(j,educ,kids,epsilon(j,educ)*theta*exp(eta_grid(eta)),SS(j,educ));
           inc_tot=inc+spouse_inc;
           
           if inc_tot>=ymin && inc_tot<ymax
               Phi_norm(a,eta,educ)=Phi_true(j,a,eta,educ,married,kids);
           end
           
       end
   end
end

% Relative population weight
aux_sum=sum(sum(sum(Phi_norm)));
Phi_mass=aux_sum/sum(sum(sum(sum(sum(sum(Phi_true))))));

% Normalize mass of individuals to 1
if aux_sum>0
    Phi_norm=Phi_norm/aux_sum;
end

Phi_p=zeros(n_agrid,n_etagrid,n_educgrid,n_kidsgrid);

clear aux_sum

% Use policy functions and survival probabilities to get distribution of idiosyncratic states in next period
for a=1:n_agrid % Assets
   for eta=1:n_etagrid % Productivity
       for educ=1:n_educgrid % Educational level

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
                   Phi_p(inds(1),etap,educ,kidsp)=Phi_p(inds(1),etap,educ,kidsp)+Phi_norm(a,eta,educ)*vals(1)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married);
                   Phi_p(inds(2),etap,educ,kidsp)=Phi_p(inds(2),etap,educ,kidsp)+Phi_norm(a,eta,educ)*vals(2)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,married);
               end
           end

       end
   end
end

% Compute value for planner
V_planner=0;

for a=1:n_agrid % Assets in next period
   for eta=1:n_etagrid % Productivity in next period
       for educ=1:n_educgrid % Educational level
           for kids=1:n_kidsgrid % No. of kids in next period
               
               wages=epsilon(j+1,educ)*theta*exp(eta_grid(eta));
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
               
               V_planner=V_planner+Phi_p(a,eta,educ,kids)*( pi_unemp(j+1,wage_ind)*V_U(j+1,a,eta,educ,married,kids,welf_checks+1)+(1-pi_unemp(j+1,wage_ind))*V_W(j+1,a,eta,educ,married,kids,welf_checks+1) );
               
           end
       end
   end
end




end