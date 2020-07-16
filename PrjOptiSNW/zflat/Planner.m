function [V_planner,Phi_mass]=Planner(Phi_true_1,j,married,kids,welf_checks,ymin,ymax,ap_ss,cutoffs,ref_earn_grid,inc_tot_grid,EV,pi_eta,pi_kids,agrid,n_agrid,n_etagrid,n_educgrid,n_kidsgrid)

%% Value function of planner

Phi_norm=zeros(n_agrid,n_etagrid,n_educgrid);

% Find mass of individuals by group (age, marital status, number of
% kids, and income range) by summing over other idiosyncratic states
for eta=1:n_etagrid % Productivity
   for educ=1:n_educgrid % Educational level

       include=(inc_tot_grid(j,:,eta,educ,married,kids)>=ymin & inc_tot_grid(j,:,eta,educ,married,kids)<ymax);
       Phi_norm(:,eta,educ)=Phi_true_1(j,:,eta,educ,married,kids).*include(:)';

   end
end

% Relative population weight
Phi_mass=sum(Phi_norm,'all')/sum(Phi_true_1,'all');

% Normalize mass of individuals to 1 and compute V_planner. Otherwise, V_planner=0 since Phi_p=0
if Phi_mass==0

    V_planner=0;

elseif Phi_mass>0

    Phi_norm=Phi_norm/Phi_mass;

    Phi_p=zeros(n_agrid,n_etagrid,n_educgrid,n_kidsgrid);

    % Use policy functions and survival probabilities to get distribution of idiosyncratic states in next period
    for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level

               if Phi_norm(a,eta,educ)>0

                   if ap_ss(j,a,eta,educ,married,kids)==0
                       inds(1)=1;
                       inds(2)=1;
                       vals(1)=1;
                       vals(2)=0;

                   elseif ap_ss(j,a,eta,educ,married,kids)>=agrid(n_agrid)
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
                           Phi_p(inds(1),etap,educ,kidsp)=Phi_p(inds(1),etap,educ,kidsp)+vals(1)*Phi_norm(a,eta,educ)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                           Phi_p(inds(2),etap,educ,kidsp)=Phi_p(inds(2),etap,educ,kidsp)+vals(2)*Phi_norm(a,eta,educ)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                       end
                   end

               end

           end
       end
    end

    % Compute value for planner
    V_planner=0;

    for a=1:n_agrid % Asset grid
        for eta=1:n_etagrid % Productivity in next period
           for educ=1:n_educgrid % Educational level
               for kids=1:n_kidsgrid % No. of kids in next period

                   if Phi_p(a,eta,educ,kids)>0

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

                       V_planner=V_planner+Phi_p(a,eta,educ,kids)*EV(j+1,a,eta,educ,married,kids,welf_checks+1,wage_ind);

                   end

               end
           end
        end
    end

else
    disp(Phi_mass)
    error('Negative mass! Check code')
end

end
