function [V_planner,Phi_mass]=Planner_grid_search(Phi_true,j,married,kids,welf_checks,ymin,ymax,V_U,V_W,ap_ss,cutoffs,ref_inc_grid,ref_earn_grid,spouse_inc_grid)

%% Value function of planner

global pi_eta pi_kids pi_unemp n_agrid n_etagrid n_educgrid n_kidsgrid

Phi_norm=zeros(n_agrid,n_etagrid,n_educgrid);

% Find mass of individuals by group (age, marital status, number of
% kids, and income range) by summing over other idiosyncratic states 
for eta=1:n_etagrid % Productivity
   for educ=1:n_educgrid % Educational level

       inc_tot=ref_inc_grid(j,1:n_agrid,eta,educ,married,kids)+spouse_inc_grid(j,1:n_agrid,eta,educ,married,kids);

       include=(inc_tot>=ymin & inc_tot<ymax);
       Phi_norm(1:n_agrid,eta,educ)=Phi_true(j,1:n_agrid,eta,educ,married,kids).*include(1:n_agrid);

   end
end

% Relative population weight
aux_sum=sum(sum(sum(Phi_norm)));
Phi_mass=aux_sum/sum(sum(sum(sum(sum(sum(Phi_true))))));

% Normalize mass of individuals to 1 and compute V_planner. Otherwise, V_planner=0 since Phi_p=0
if aux_sum==0
    
    V_planner=0;    
    
elseif aux_sum>0
    
    Phi_norm=Phi_norm/aux_sum;
    
    Phi_p=zeros(n_agrid,n_etagrid,n_educgrid,n_kidsgrid);

    clear aux_sum

    % Use policy functions and survival probabilities to get distribution of idiosyncratic states in next period
    for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
                
               for etap=1:n_etagrid
                   for kidsp=1:n_kidsgrid
                       Phi_p(ap_ss(j,a,eta,educ,married,kids),etap,educ,kidsp)=Phi_p(ap_ss(j,a,eta,educ,married,kids),etap,educ,kidsp)+Phi_norm(a,eta,educ)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                   end
               end
           end
       end
    end


    % Compute value for planner
    V_planner=0;

    for eta=1:n_etagrid % Productivity in next period
       for educ=1:n_educgrid % Educational level
           for kids=1:n_kidsgrid % No. of kids in next period

               wages=ref_earn_grid(j+1,1,eta,educ,married,kids); % Wages/labor earnings are independent of assets (setting equal to 1 without loss of generality). Also independent of marital status and number of children
               
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

               V_planner=V_planner+Phi_p(1:n_agrid,eta,educ,kids)'*( pi_unemp(j+1,wage_ind)*V_U(j+1,1:n_agrid,eta,educ,married,kids,welf_checks+1)'+(1-pi_unemp(j+1,wage_ind))*V_W(j+1,1:n_agrid,eta,educ,married,kids,welf_checks+1)' );

           end
       end
    end
    
else
    disp(aux_sum)
    error('Negative mass! Check code') 
end

end