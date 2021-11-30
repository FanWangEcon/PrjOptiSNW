%% SNW_calibrate_2009 UI Support for Lost Wages Calibration
%    The ratio of UI benefits to wages and salary is 2.1 percent in 2009.
%    xi in [0,1] governs the duration of unemployment shock for those
%    unemployed. This equals to 0.532 in 2009 (xi=0 no wages earned).
%
%    We solve for total wage earnings from unemployed and employed in 2009,
%    for employed, same as under steady-state. For unemployed, they lose
%    (1-xi) share of the wage they would otherwise have earned. Our
%    unemployment probability in 2009 is conditional on age and edu groups
%    (SNW_UNEMP_2008.m) computed based on rectiilnear restriction.
%
%    We know total UI amount (multiply its share of total "Wages and
%    salary" by total "wages and salary". We know how much wage was lost
%    due to xi. The ratio of these two levels is b, which is the
%    parameter that is the share of lost-wage recovered. Note that this is
%    based on exogenous wage earnings, so we do not have to worry about
%    endogenous changes to savings. We will solve for the steady-state
%    distribution, which generates mass of people by age, education,
%    marital status, kids count, etc.
%
%    [FL_B_CALIBRATED_BY_UI_SHARE, MP_STATS_WAGE_UI_SPENDING,
%    MN_EARN_TOT_WGTED, MN_EARN_UNEMP_WGTED, MN_EARN_UNEMP_TOT_WGTED,
%    MN_EARN_UNEMP_WEIGHTED_WGTED] = SNW_CALIBRATE_2009_B(MP_PARAMS,
%    MP_CONTROLS, FL_RATIO_UI_BENEFITS_TO_WAGE) where
%    FL_RATIO_UI_BENEFITS_TO_WAGE is the ratio of UI benefits to teh total
%    wages and salary bill, including wage for the household head and
%    pouse, considering unemployment probability by states and loss of
%    income due to unemployment duration.
%
%    See also SNWX_V0808_JAEEMK, SNW_V0808_JAEEMK, SNWX_UNEMP_2008,
%    SNW_UNEMP_2008, SNW_UNEMP_2008
%

%%
function [varargout]=snw_calibrate_2009_b(varargin)

%% Default and Parse
if (~isempty(varargin))

    fl_ratio_ui_benefits_to_wage = 0.021;
    if (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==3)
        [mp_params, mp_controls, fl_ratio_ui_benefits_to_wage] = varargin{:};
    else
        error('Need to provide 2/3 parameter inputs');
    end

else
    clc;
    close all;

    % Solve the VFI Problem and get Value Function
    mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
    mp_more_inputs('fl_ss_non_college') = 0.225;
    mp_more_inputs('fl_ss_college') = 0.271;
    mp_more_inputs('fl_scaleconvertor') = 54831;
    st_param_group = 'default_small';
%     st_param_group = 'default_dense';
    mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
    mp_controls = snw_mp_control('default_test');

    % no b, solving for b, b set to 0 when solving for wages
    % set Unemployment Related Variables
    xi=0.532; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    % TR parameter does not matter below
    TR=100/54831; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values

    mp_params('xi') = xi;
    mp_params('TR') = TR;

    % Solve for Unemployment Values
    mp_controls('bl_print_calibrate_2009') = true;
    mp_controls('bl_print_calibrate_2009_verbose') = true;

    fl_ratio_ui_benefits_to_wage = 0.021;

end

%% Set b = 0
% solving for b, b=0 to get wages without b
fl_b_zero = 0;
mp_params('b') = 0;

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

params_group = values(mp_params, {'pi_unemp_2009_edu_age'});
[pi_unemp_2009_edu_age] = params_group{:};

params_group = values(mp_params, {'xi','b'});
[xi, b] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_calibrate_2009', 'bl_print_calibrate_2009_verbose'});
[bl_print_calibrate_2009, bl_print_calibrate_2009_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Solve for distribution
% the savings distirbution will not impact later calculations, since
% calcualtions based on wages only, not savings. But will compute the
% distributions to get the joint distribution of age, children count,
% marital status, etc.

[V_ss, ap_ss, cons_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
[Phi_true, Phi_adj_ss] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss);


%% Solve 08 Policy and Value given 09 Unemployment shock
% 09 has two states, employed or not (extensive), and if unemployed (intensive duration)

% mn_earn_tot = wage earning household head + household-head-spouse (no interest no SS)
mn_earn_tot = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% mn_earn_unemp = household head only wage earning under unemployment
mn_earn_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% mn_earn_unemp_tot = wage earning household head + household-head-spouse (no interest no SS, no UI benefits) under unemployment
mn_earn_unemp_tot = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% mn_earn_unemp_weighted = weighted sum of mn_earn_tot and mn_earn_unemp_tot
mn_earn_unemp_weighted = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Loop Over States and Pre-Store
for j=1:n_jgrid
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid

                        % 1. Get unemployment probability conditional on age and edu
                        % columns are edu: 1st column high school; 2nd col college
                        % rows are age groups: 18--24; 25-54; 55 to 65
                        if (j<=7)
                            it_age_grp = 1;
                        elseif (j<=37)
                            it_age_grp = 2;
                        else
                            it_age_grp = 3;
                        end
                        fl_unemp_edu_age = pi_unemp_2009_edu_age(it_age_grp, educ);

                        % 2. Income and Wages employed
                        % inc = SS + wages + interest earnings + bequest
                        % earn = wages
                        [inc, earn]=snw_hh_individual_income(j,a,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        % total household income and earnings
                        earn_tot = earn+(married-1)*spouse_inc*exp(eta_S_grid(eta));

                        % 3. Income unemployed
                        % inc_umemp = SS + wages*(xi+b*(1-xi)) + interest earnings + bequest
                        % earn_unemp = wages*(xi+b*(1-xi)), NOT multiplied by (xi+b*(1-xi))
                        [inc_umemp,earn_unemp]=snw_hh_individual_income(j,a,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
                            xi,b);
                        % We assume spousal income NOT impcated by (xi+b*(1-xi))
                        spouse_inc_unemp=snw_hh_spousal_income(j,educ,kids,earn_unemp,SS(j,educ), jret);
                        % Total household income and earnings under unemployment
                        % under earn_unemp_tot, note UI = 0, b = 0
                        % xi=0 means no wage earnings in 2009
                        earn_unemp_tot = earn_unemp*(xi)+(married-1)*spouse_inc_unemp*exp(eta_S_grid(eta));

                        % 4. Collect results to mn Matrixes
                        mn_earn_tot(j,a,eta,educ,married,kids) = earn_tot*(1 - fl_unemp_edu_age);
                        mn_earn_unemp(j,a,eta,educ,married,kids) = (earn_unemp*(xi))*fl_unemp_edu_age;
                        mn_earn_unemp_tot(j,a,eta,educ,married,kids) = earn_unemp_tot*fl_unemp_edu_age;
                        mn_earn_unemp_weighted(j,a,eta,educ,married,kids) = ...
                            fl_unemp_edu_age*earn_unemp_tot + (1 - fl_unemp_edu_age)*earn_tot;

                    end
                end
            end
        end
    end
end

%% Compute weighted values, weight with mass by states
mn_earn_tot_wgted = Phi_true.*mn_earn_tot;
mn_earn_unemp_wgted = Phi_true.*mn_earn_unemp;
mn_earn_unemp_tot_wgted = Phi_true.*mn_earn_unemp_tot;
mn_earn_unemp_weighted_wgted = Phi_true.*mn_earn_unemp_weighted;

%% Compute b: The Ratio of UI Benefits to Wage is 2.1 percent in 2009
% 1. we know the total wages, household head and household spouse
fl_total_wage = sum(mn_earn_unemp_weighted_wgted,'all');

% 2. what is the total b spending?
% fl_ratio_ui_benefits_to_wage = 0.021;
fl_total_b_spending = fl_total_wage*fl_ratio_ui_benefits_to_wage;

% 3. Total wage earning by household head unemployed (with unemployment duration xi)
fl_total_wage_unemp_hhhead = sum(mn_earn_unemp_wgted, 'all');

% 4. Wages lost for unemployed household head
fl_total_wage_unemp_hhhead_lost = (fl_total_wage_unemp_hhhead/xi)*(1-xi);

% 5. Compute the b, share of lost income recovered parameter.
fl_b_calibrated_by_ui_share = fl_total_b_spending/fl_total_wage_unemp_hhhead_lost;

% Gather statistics
mp_stats_wage_ui_spending = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_stats_wage_ui_spending('fl_total_wage') = fl_total_wage;
mp_stats_wage_ui_spending('fl_total_b_spending') = fl_total_b_spending;
mp_stats_wage_ui_spending('fl_total_wage_unemp_hhhead') = fl_total_wage_unemp_hhhead;
mp_stats_wage_ui_spending('fl_total_wage_unemp_hhhead_lost') = fl_total_wage_unemp_hhhead_lost;
mp_stats_wage_ui_spending('fl_b_calibrated_by_ui_share') = fl_b_calibrated_by_ui_share;

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_calibrate_2009", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_calibrate_2009)
    ff_container_map_display(mp_stats_wage_ui_spending, 9, 9);
    if (bl_print_calibrate_2009_verbose)
        mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
        mp_outcomes('mn_earn_tot_wgted') = mn_earn_tot_wgted;
        mp_outcomes('mn_earn_unemp_wgted') = mn_earn_unemp_wgted;
        mp_outcomes('mn_earn_unemp_tot_wgted') = mn_earn_unemp_tot_wgted;
        mp_outcomes('mn_earn_unemp_weighted_wgted') = mn_earn_unemp_weighted_wgted;
        ff_container_map_display(mp_outcomes, 9, 9);
    end
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = fl_b_calibrated_by_ui_share;
    elseif (it_k==2)
        ob_out_cur = mp_stats_wage_ui_spending;
    elseif (it_k==3)
        ob_out_cur = mn_earn_tot_wgted;
    elseif (it_k==4)
        ob_out_cur = mn_earn_unemp_wgted;
    elseif (it_k==5)
        ob_out_cur = mn_earn_unemp_tot_wgted;
    elseif (it_k==6)
        ob_out_cur = mn_earn_unemp_weighted_wgted;
    end
    varargout{it_k} = ob_out_cur;
end

end
