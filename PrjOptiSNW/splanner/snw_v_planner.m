%% SNW_V_PLANNER Planner Last Year Value Given Income, Marriage, Age, Kids
%    Obtains the mass at each age, marriage status, child count, and income
%    level. And calculates planner value integrating over V_U and V_W,
%    where the state-space fall within the planner's
%    age/marriage/child-count/income-level group, and V_U and V_W are
%    specific to each check. Calculate the planner Value associated with
%    each check level.
%
%    * ST_SOLU_TYPE: 'matlab_minimizer', 'grid_search', 'bisec_vec', three
%    alternative solution methods
%
%    [V_VFI, ap_VFI, cons_VFI] = SNW_VU_VW_CHECKS(ST_SOLU_TYPE) Solves
%    given solution method for value function, policy functions. Using
%    default parameters.
%
%    [V_VFI, ap_VFI, cons_VFI, V_VFI_unemp, ap_VFI_unemp, cons_VFI_unemp] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS) Solves given
%    solution method for value function, policy functions under both
%    employed and unemployed situations, using default parameters. Change
%    model parameters through MP_PARAMS.
%
%    [V_VFI, ap_VFI, cons_VFI, V_VFI_unemp, ap_VFI_unemp, cons_VFI_unemp,
%    V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS) Solves the
%    model for some number of checks for employed and unemployed. The
%    number of checks is set through the 'n_welfchecksgrid' parameters in
%    MP_PARAMS.
%
%    [V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS, V_SS, V_UNEMP)
%    provide V_SS and V_UNEMP solved out elsewhere, and only get two
%    outputs out.
%
%    See also SNW_VU_VW_CHECKS
%

%%
function [varargout]=snw_v_planner(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==8)
        [mp_params, mp_controls, Phi_true, ap_ss, ...
            V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] = varargin{:};
    else
        error('Need to provide 8 parameter inputs');
    end

else
    close all;

    st_solu_type = 'bisec_vec';
    mp_params = snw_mp_param('default_base');
    mp_controls = snw_mp_control('default_test');

    mp_params('TR') = 100/58056;
    mp_params('xi') = 0.5;
    mp_params('b') = 0;
    mp_params('n_welfchecksgrid') = 51;
    n_incgrid = 201;
    mp_params('n_incgrid') = n_incgrid;
    mp_params('inc_grid') = linspace(0, 7, n_incgrid)'; % 7 refers to 7*58056=406392 dollars in 2012USD

    mp_controls('bl_timer') = true;

    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_vu_vw') = false;

    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_vu_vw_verbose') = false;

    mp_controls('bl_print_v_planner') = true;
    mp_controls('bl_print_v_planner_verbose') = true;

    % Calculate Values for All Checks
    [~, ap_ss, ~, ...
        ~, ~, ~, ...
        V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] = ...
        snw_vu_vw_checks(st_solu_type, mp_params, mp_controls);

    % Distribution
    [Phi_true] = snw_ds_main(mp_params, mp_controls);

end

%% Reset All globals
% % Parameters used in this code directly
% global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% % Used in find_a_working function
% global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option throw_in_ocean
% global pi_eta pi_kids pi_unemp

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r' , 'jret'});
[theta,  r, jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids'});
[pi_eta, pi_kids] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'n_welfchecksgrid', 'scaleconvertor'});
[n_welfchecksgrid, scaleconvertor] = params_group{:};

params_group = values(mp_params, {'pi_unemp', 'n_incgrid', 'inc_grid'});
[pi_unemp, n_incgrid, inc_grid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_v_planner', 'bl_print_v_planner_verbose'});
[bl_print_v_planner, bl_print_v_planner_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tm_planner_all = tic;
end

%% Joint Mass 

Phi_true_1=Phi_true/sum(Phi_true,'all');

%% Precompute
cutoffs=snw_wage_cutoffs(Phi_true, ...
    theta, epsilon, eta_H_grid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid, jret);

if (bl_print_v_planner) tm_planner_inc_compute = tic; end

% Income and earning grids used to speed up the code
ref_earn_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
inc_tot_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
inc_tot_ygroup_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Ap lower and upper index
ap_idx_lower_ss = NaN(size(ap_ss));
ap_idx_higher_ss = NaN(size(ap_ss));
ap_idx_higher_weight_ss = NaN(size(ap_ss));

% EV Store, expected check values
EVU=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
EVC=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);

% Loop Over States and Pre-Store
for j=1:n_jgrid
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid
                        
                        % 1. Gather income Information
                        % [inc,earn]=individual_income(j,a,eta,educ);
                        % spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                        [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        
                        ref_earn_grid(j,a,eta,educ,married,kids)=earn;                        
                        int_tot = inc+spouse_inc;
                        inc_tot_grid(j,a,eta,educ,married,kids)=int_tot;
                        
                        % 2. Which Income Group does total income belong to?
                        inc_group_include_idx = NaN;
                        for inc_group=1:n_incgrid
                            if inc_group<n_incgrid
                                ymin = inc_grid(inc_group);
                                ymax = inc_grid(inc_group+1);
                            elseif inc_group==n_incgrid
                                ymin = inc_grid(inc_group);
                                ymax = 10E30;
                            end
                            if (int_tot>=ymin && int_tot<ymax)
                                inc_group_include_idx = inc_group;
                            end
                        end
                        if (Phi_true_1(j,a,eta,educ,married,kids)>0)
                            inc_tot_ygroup_grid(j,a,eta,educ,married,kids)=inc_group_include_idx;
                        else
                            % inc group as 0, not a valid group if no mass
                            inc_tot_ygroup_grid(j,a,eta,educ,married,kids)=0;
                        end
                        
                        % 3. Belongs to which wage group?
                        wages = ref_earn_grid(j,a,eta,educ,married,kids);
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
                        
                        % 4. Expected Value of Checks
                        for welf_checks=1:n_welfchecksgrid
                            EVU(j,a,eta,educ,married,kids,welf_checks)=pi_unemp(j,wage_ind)*V_U_allchecks(j,a,eta,educ,married,kids,welf_checks)+(1-pi_unemp(j,wage_ind))*V_W_allchecks(j,a,eta,educ,married,kids,welf_checks);
                            EVC(j,a,eta,educ,married,kids,welf_checks)=pi_unemp(j,wage_ind)*C_U_allchecks(j,a,eta,educ,married,kids,welf_checks)+(1-pi_unemp(j,wage_ind))*C_W_allchecks(j,a,eta,educ,married,kids,welf_checks);
                        end
                        
                        % 5. Precalculate E(ap_higher) and E(ap_lower)
                        if ap_ss(j,a,eta,educ,married,kids)==0
                            inds(1)=1;
                            inds(2)=1;
                            vals(1)=1;
                        elseif ap_ss(j,a,eta,educ,married,kids)>=agrid(n_agrid)
                            inds(1)=n_agrid;
                            inds(2)=n_agrid;
                            vals(1)=1;
                        else                            
                            ind_aux=find(agrid<=ap_ss(j,a,eta,educ,married,kids),1,'last');                            
                            inds(1)=ind_aux;
                            inds(2)=ind_aux+1;                            
                            % Linear interpolation
                            vals(1)=1-((ap_ss(j,a,eta,educ,married,kids)-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                        end                        
                        ap_idx_lower_ss(j,a,eta,educ,married,kids) = inds(1);
                        ap_idx_higher_ss(j,a,eta,educ,married,kids) = inds(2);
                        ap_idx_higher_weight_ss(j,a,eta,educ,married,kids) = vals(1);
                        
                    end
                end
            end
        end
    end
end

if (bl_print_v_planner)
    tm_planner_inc_end = toc(tm_planner_inc_compute);
    disp(strcat(['SNW_V_PLANNER: pre-computed incomes '...
        ', time-cost:' num2str(tm_planner_inc_end)]));
end

%% Solve planner's problem

C_planner=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid);
V_planner=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid);
Phi_mass=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

for c_or_v=[1,2]
    
    if (bl_print_v_planner)
        disp(strcat(['SNW_V_PLANNER: c(1) or v(2):' num2str(c_or_v)]));
    end
    
    if (c_or_v == 1)
        EV = EVU;
    elseif(c_or_v == 2)
        EV = EVC;
    end
    
    for j=1:(n_jgrid-1) % Age
        
        if (bl_print_v_planner) tm_planner_age = tic; end
        
        for married=1:n_marriedgrid % Marital status
            for kids=1:n_kidsgrid % Number of kids
                
                % income groups that are reached by this age/marry/kids group
                Phi_true_1_amk = Phi_true_1(j,:,:,:,married,kids);
                inc_tot_ygroup_amk = inc_tot_ygroup_grid(j,:,:,:,married,kids);
                inc_tot_ygroup_amk_unique = unique(inc_tot_ygroup_amk(:))';
                
                % loop over income groups that are specific to this amk
                for inc_group=inc_tot_ygroup_amk_unique
                    
                    % zero wage on in the income group, but no mass
                    if(inc_group~=0)
                        
                        % Mass at this amk point
                        Phi_true_1_amk_in_inc_group = (inc_tot_ygroup_amk == inc_group);
                        Phi_mass_amk = sum(Phi_true_1_amk(Phi_true_1_amk_in_inc_group), 'all');
                        Phi_mass(j,married,kids,inc_group) = Phi_mass_amk;
                        
                        % Mass for inc
                        Phi_true_1_amk_inc = zeros(size(Phi_true_1_amk));
                        Phi_true_1_amk_inc(Phi_true_1_amk_in_inc_group) = Phi_true_1_amk(Phi_true_1_amk_in_inc_group);
                        
                        % sets of (a, eta, edu) in this income group.
                        [~,a_idx,eta_idx,educ_idx] = ind2sub(size(inc_tot_ygroup_amk), find(Phi_true_1_amk_in_inc_group));
                        % do not need to loop over all a, edu, educ, only
                        % the following combinations
                        a_eta_educ_in_inc_group_idx = [a_idx,eta_idx,educ_idx];
                        
                        for welf_checks=0:(n_welfchecksgrid-1)
                            
                            % Evaluate for Planner
                            [V_or_C] = snw_v_planner_amky(...
                                Phi_true_1_amk_inc, Phi_mass_amk,...
                                a_eta_educ_in_inc_group_idx,...
                                j,married,kids,welf_checks,...
                                ap_idx_lower_ss,ap_idx_higher_ss,ap_idx_higher_weight_ss,...
                                EV,...
                                pi_eta,pi_kids,n_etagrid,n_kidsgrid);
                            
                            % Store Results
                            if (c_or_v == 1)
                                V_planner(j,married,kids,welf_checks+1,inc_group) = V_or_C;
                            elseif(c_or_v == 2)
                                C_planner(j,married,kids,welf_checks+1,inc_group) = V_or_C;
                            end
                            
                        end
                    end
                end
            end
        end
        
        if (bl_print_v_planner)
            tm_planner_age_end = toc(tm_planner_age);
            disp(strcat(['SNW_V_PLANNER: Finished Age Group:' ...
                num2str(j) ' of ' num2str(n_jgrid-1) ...
                ', time-this-age:' num2str(tm_planner_age_end)]));
        end
        
    end
end

%% Output for computing optimal allocation
Output=zeros((n_jgrid-1)*n_marriedgrid*n_kidsgrid*n_welfchecksgrid*n_incgrid,9);
counter=0;
inc_grid = [inc_grid; 10E30];
for inc_group=1:n_incgrid
    for welf_checks=0:(n_welfchecksgrid-1) % Number of welfare checks
        for kids=1:n_kidsgrid % Number of kids
            for married=1:n_marriedgrid % Marital status
                for j=1:(n_jgrid-1) % Age
                    
                    counter=counter+1;
                    
                    Output(counter,1)=17+j;
                    Output(counter,2)=married-1;
                    Output(counter,3)=kids-1;
                    Output(counter,4)=welf_checks;
                    
                    Output(counter,5)=inc_grid(inc_group);
%                     Output(counter,6)=inc_grid(inc_group+1)*58056;
                    Output(counter,6)=Phi_mass(j,married,kids,inc_group);
                    Output(counter,7)=psi(j);
                    %                     Output(counter,9)=V_planner(j,married,kids,welf_checks+1,inc_group);
                    %                     Output(counter,10)=C_planner(j,married,kids,welf_checks+1,inc_group);
                    
                end
            end
        end
    end
end

Output(:,8) = V_planner(:);
Output(:,9) = C_planner(:);

%% Timing and Profiling End
if (bl_timer)
    tm_planner_all_end = toc(tm_planner_all);
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_V_PLANNER", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time cost=' num2str(tm_planner_all_end)] ...
        ], ";");
    disp(st_complete_vu_vw_checks);
end

%% Print
if (bl_print_v_planner_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('V_planner') = V_planner;
    mp_outcomes('C_planner') = C_planner;
    mp_outcomes('Phi_mass') = Phi_mass;
    mp_outcomes('Output') = Output;
    ff_container_map_display(mp_outcomes);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_planner;
    elseif (it_k==2)
        ob_out_cur = C_planner;
    elseif (it_k==3)
        ob_out_cur = Phi_mass;
    elseif (it_k==4)
        ob_out_cur = Output;
    end
    varargout{it_k} = ob_out_cur;
end

end
