%% SNW_HH_PRECOMPUTE Pre-computes Reusable Household Inc and Other Stats
%    For each check, for the planner, often need to re-compute household
%    head, household spouse income, and other statistics. Pre-computes
%    these matrixes to save time.
%
%    * CL_ST_COMPUTE_LIST: a list of strings of what are pre-computed
%
%    [MP_PRECOMPUTE_RES] = SNW_HH_PRECOMPUTE(MP_PARAMS, MP_CONTROLS,
%    CL_ST_COMPUTE_LIST) given the list of results to to pre-compute, show
%    a container map of results.
%
%    [MP_PRECOMPUTE_RES] = SNW_HH_PRECOMPUTE(MP_PARAMS, MP_CONTROLS,
%    CL_ST_COMPUTE_LIST, AP_SS, PHI_TRUE) if ref_earn_wageind_grid is
%    included in CL_ST_COMPUTE_LIST, also include probability mass
%    PHI_TRUE. At the same time, include optimal savings choice AP_SS,
%    assuming we have 'ap_idx_lower_ss', 'ap_idx_higher_ss',
%    'ap_idx_higher_weight_ss' in CL_ST_COMPUTE_LIST.
%

%%
function [mp_precompute_res]=snw_hh_precompute(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==3)
        [mp_params, mp_controls, cl_st_precompute_list] = varargin{:};
    elseif (length(varargin)==5)
        [mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true] = varargin{:};
    else
        error('Need to provide 3/5 parameter inputs');
    end

else

    close all;

    mp_params = snw_mp_param('default_small');
    mp_controls = snw_mp_control('default_test');

    % set Unemployment Related Variables
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)

    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;

%     cl_st_precompute_list = {'inc', 'spouse_inc', 'a',...
%         'ref_earn_grid', 'inc_tot_grid', 'inc_tot_ygroup_grid'};
    cl_st_precompute_list = {'inc', 'inc_unemp', ...
        'spouse_inc', 'spouse_inc_unemp', ...
        'a', 'ref_earn_grid', 'ref_earn_wageind_grid', ...
        'inc_tot_grid', 'inc_tot_ygroup_grid',...
        'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss',...
        'ar_z_ctr_amz'};

    [~,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);

end

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r' , 'jret'});
[theta,  r, jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp', 'n_incgrid', 'inc_grid'});
[pi_unemp, n_incgrid, inc_grid] = params_group{:};

params_group = values(mp_params, {'n_welfchecksgrid', 'xi','b'});
[n_welfchecksgrid, xi, b] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_precompute', 'bl_print_precompute_verbose'});
[bl_print_precompute, bl_print_precompute_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Initialize

% Used in Check VU VW calculations
if(sum(strcmp(cl_st_precompute_list, 'inc')))
    mn_inc = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'inc_unemp')))
    mn_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'spouse_inc')))
    mn_spouse_inc = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'spouse_inc_unemp')))
    mn_spouse_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'a')))
    mn_a = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end

% Might be used in planner calculation
if(sum(strcmp(cl_st_precompute_list, 'ref_earn_grid')))
    ref_earn_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'ref_earn_wageind_grid')))
    ref_earn_wageind_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'inc_tot_grid')))
    inc_tot_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'inc_tot_ygroup_grid')))
    inc_tot_ygroup_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end

% Ap lower and upper choice index and weight from regulat VFI problem
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_ss')))
    ap_idx_lower_ss = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_higher_ss')))
    ap_idx_higher_ss = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_weight_ss')))
    ap_idx_lower_weight_ss = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end

% Z Counter for FOC based EV Interpolation
if(sum(strcmp(cl_st_precompute_list, 'ar_z_ctr_amz')))
    ar_z_ctr_amz = NaN(n_agrid*n_kidsgrid*n_marriedgrid*n_educgrid*n_etagrid,1);
end

%% Cutoffs
cutoffs = snw_wage_cutoffs(Phi_true, theta, epsilon, eta_H_grid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid, jret);

%% Pre-Store z-Counter
if(sum(strcmp(cl_st_precompute_list, 'ar_z_ctr_amz')))
    % Within each age, a, kids, married, edu, eta are the same
    mn_z_ctr = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    for a=1:n_agrid % Assets
        it_z_ctr = 0;
        % order here is reverse of order of loops at L181, so that it_z_ctr
        % match up with column order after reshape(mn_ev_ap_z, n_agrid, [])
        for kids=1:n_kidsgrid % Number of kids
            for married=1:n_marriedgrid % Marital status
                for educ=1:n_educgrid % Educational level
                    for eta=1:n_etagrid % Productivity

                        % Non-asset position Counter
                        it_z_ctr = it_z_ctr + 1;
                        mn_z_ctr(a, eta, educ, married, kids) = it_z_ctr;

                    end
                end
            end
        end
    end
    ar_z_ctr_amz = mn_z_ctr(:);
end

%% Pre-Store a eta educ married kids sequence
% Loop Over States and Pre-Store
for j=1:n_jgrid

    if (bl_print_precompute) tm_precompute_age = tic; end

    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid

                        % 1. Gather income Information under employment
                        [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        int_tot = inc+(married-1)*spouse_inc*exp(eta_S_grid(eta));

                        if(sum(strcmp(cl_st_precompute_list, 'a')))
                            mn_a(j,a,eta,educ,married,kids) = agrid(a);
                        end
                        if(sum(strcmp(cl_st_precompute_list, 'inc')))
                            mn_inc(j,a,eta,educ,married,kids) = inc;
                        end
                        if(sum(strcmp(cl_st_precompute_list, 'spouse_inc')))
                            mn_spouse_inc(j,a,eta,educ,married,kids) = (married-1)*spouse_inc*exp(eta_S_grid(eta));
                        end
                        if(sum(strcmp(cl_st_precompute_list, 'ref_earn_grid')) || ...
                            sum(strcmp(cl_st_precompute_list, 'ref_earn_wageind_grid')))
                            ref_earn_grid(j,a,eta,educ,married,kids) = earn;
                        end
                        if(sum(strcmp(cl_st_precompute_list, 'inc_tot_grid')))
                            inc_tot_grid(j,a,eta,educ,married,kids) = int_tot;
                        end

                        % 2. Which Income Group does total income belong to?
                        if(sum(strcmp(cl_st_precompute_list, 'inc_tot_ygroup_grid')))
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
                            if (Phi_true(j,a,eta,educ,married,kids)>0)
                                inc_tot_ygroup_grid(j,a,eta,educ,married,kids)=inc_group_include_idx;
                            else
                                % inc group as 0, not a valid group if no mass
                                inc_tot_ygroup_grid(j,a,eta,educ,married,kids)=0;
                            end
                        end

                        % 3. Unemployment Income Information
                        if (sum(strcmp(cl_st_precompute_list, 'inc_unemp')) || ...
                            sum(strcmp(cl_st_precompute_list, 'spouse_inc_unemp')))
                            [inc_umemp,earn_unemp]=snw_hh_individual_income(j,a,eta,educ,...
                                theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
                                xi,b);
                            spouse_inc_unemp=snw_hh_spousal_income(j,educ,kids,earn_unemp,SS(j,educ), jret);

                            if(sum(strcmp(cl_st_precompute_list, 'inc_unemp')))
                                mn_inc_unemp(j,a,eta,educ,married,kids) = inc_umemp;
                            end
                            if(sum(strcmp(cl_st_precompute_list, 'spouse_inc_unemp')))
                                mn_spouse_inc_unemp(j,a,eta,educ,married,kids) = (married-1)*spouse_inc_unemp*exp(eta_S_grid(eta));
                            end

                        end

                        % 4. Belongs to which wage group?
                        if(sum(strcmp(cl_st_precompute_list, 'ref_earn_wageind_grid')))
                            if earn<=cutoffs(1)
                                wage_ind=1;
                            elseif earn>cutoffs(1) && earn<=cutoffs(2)
                                wage_ind=2;
                            elseif earn>cutoffs(2) && earn<=cutoffs(3)
                                wage_ind=3;
                            elseif earn>cutoffs(3) && earn<=cutoffs(4)
                                wage_ind=4;
                            elseif earn>cutoffs(4)
                                wage_ind=5;
                            end
                            ref_earn_wageind_grid(j,a,eta,educ,married,kids)=wage_ind;
                        end

                        % 5. Precalculate E(ap_higher) and E(ap_lower)
                        if (sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_ss')) || ...
                            sum(strcmp(cl_st_precompute_list, 'ap_idx_higher_ss')) || ...
                            sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_weight_ss')))

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

                            if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_ss')))
                                ap_idx_lower_ss(j,a,eta,educ,married,kids) = inds(1);
                            end
                            if(sum(strcmp(cl_st_precompute_list, 'ap_idx_higher_ss')))
                                ap_idx_higher_ss(j,a,eta,educ,married,kids) = inds(2);
                            end
                            if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_weight_ss')))
                                ap_idx_lower_weight_ss(j,a,eta,educ,married,kids) = vals(1);
                            end
                        end

                    end
                end
            end
        end
    end

    if (bl_print_precompute)
        tm_precompute_age_end = toc(tm_precompute_age);
        disp(strcat(['SNW_HH_PRECOMPUTE: Finished Age Group:' ...
            num2str(j) ' of ' num2str(n_jgrid-1) ...
            ', time-this-age:' num2str(tm_precompute_age_end)]));
    end

end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_HH_PRECOMPUTE", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time cost=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete_vu_vw_checks);
end

%% Return
mp_precompute_res = containers.Map('KeyType', 'char', 'ValueType', 'any');
if(sum(strcmp(cl_st_precompute_list, 'inc')))
    ar_inc_amz = mn_inc(:);
    mp_precompute_res('ar_inc_amz') = ar_inc_amz;
end
if(sum(strcmp(cl_st_precompute_list, 'inc_unemp')))
    ar_inc_unemp_amz = mn_inc_unemp(:);
    mp_precompute_res('ar_inc_unemp_amz') = ar_inc_unemp_amz;
end
if(sum(strcmp(cl_st_precompute_list, 'spouse_inc')))
    ar_spouse_inc_amz = mn_spouse_inc(:);
    mp_precompute_res('ar_spouse_inc_amz') = ar_spouse_inc_amz;
end
if(sum(strcmp(cl_st_precompute_list, 'spouse_inc_unemp')))
    ar_spouse_inc_unemp_amz = mn_spouse_inc_unemp(:);
    mp_precompute_res('ar_spouse_inc_unemp_amz') = ar_spouse_inc_unemp_amz;
end
if(sum(strcmp(cl_st_precompute_list, 'a')))
    ar_a_amz = mn_a(:);
    mp_precompute_res('ar_a_amz') = ar_a_amz;
end

% Might be used in planner calculation
if(sum(strcmp(cl_st_precompute_list, 'ref_earn_grid')))
    mp_precompute_res('ref_earn_grid') = ref_earn_grid;
end
if(sum(strcmp(cl_st_precompute_list, 'ref_earn_wageind_grid')))
    mp_precompute_res('ref_earn_wageind_grid') = ref_earn_wageind_grid;
end
if(sum(strcmp(cl_st_precompute_list, 'inc_tot_grid')))
    mp_precompute_res('inc_tot_grid') = inc_tot_grid;
end
if(sum(strcmp(cl_st_precompute_list, 'inc_tot_ygroup_grid')))
    mp_precompute_res('inc_tot_ygroup_grid') = inc_tot_ygroup_grid;
end

% Ap lower and upper choice index and weight from regulat VFI problem
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_ss')))
    mp_precompute_res('ap_idx_lower_ss') = ap_idx_lower_ss;
end
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_higher_ss')))
    mp_precompute_res('ap_idx_higher_ss') = ap_idx_higher_ss;
end
if(sum(strcmp(cl_st_precompute_list, 'ap_idx_lower_weight_ss')))
    mp_precompute_res('ap_idx_lower_weight_ss') = ap_idx_lower_weight_ss;
end

% Z Counter for FOC based EV Interpolation
if(sum(strcmp(cl_st_precompute_list, 'ar_z_ctr_amz')))
    mp_precompute_res('ar_z_ctr_amz') = ar_z_ctr_amz;
end

%% Print
if (bl_print_precompute_verbose)
    ff_container_map_display(mp_precompute_res, 9, 9);
end

end
