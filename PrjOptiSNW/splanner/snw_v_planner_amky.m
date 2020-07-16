function [V_planner]=snw_v_planner_amky(...
    Phi_true_1_amk, Phi_mass_amk,...
    a_eta_educ_in_inc_group_idx,...
    j,married,kids,...
    welf_checks,...
    ap_idx_lower_ss,ap_idx_higher_ss,ap_idx_higher_weight_ss,...
    EV,...
    pi_eta,pi_kids,...
    n_etagrid,n_kidsgrid)

%% Value function of planner

V_planner = 0;

% Use policy functions and survival probabilities to get distribution of idiosyncratic states in next period
for it_ctr=1:size(a_eta_educ_in_inc_group_idx,1)
    
    % Get Index
    a = a_eta_educ_in_inc_group_idx(it_ctr,1);
    eta = a_eta_educ_in_inc_group_idx(it_ctr,2);
    educ = a_eta_educ_in_inc_group_idx(it_ctr,3);
            

    inds(1)=ap_idx_lower_ss(j,a,eta,educ,married,kids);
    inds(2)=ap_idx_higher_ss(j,a,eta,educ,married,kids);

    % Linear interpolation
    vals(1)=ap_idx_higher_weight_ss(j,a,eta,educ,married,kids);
    vals(2)=1-vals(1);
    
    condi_mass = Phi_true_1_amk(1,a,eta,educ)/Phi_mass_amk;
    
    for etap=1:n_etagrid
        for kidsp=1:n_kidsgrid
            V_planner=V_planner ...
                + vals(1)*condi_mass*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)*EV(j+1,inds(1),etap,educ,married,kidsp,welf_checks+1)...
                + vals(2)*condi_mass*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)*EV(j+1,inds(2),etap,educ,married,kidsp,welf_checks+1);
            
        end
    end           
end

end
