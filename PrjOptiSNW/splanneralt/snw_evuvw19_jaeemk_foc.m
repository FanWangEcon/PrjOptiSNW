%% SNW_EVUVW19_JAEEMK_FOC Solves 2019 and 2007 Expected Values Given 2020 and 2008 EVs/Vs
%    Despite the name, this problem supports solving the 2019 looking into
%    2020 as well as the 2007 looking into 2008 problems. The idea is that
%    the planner only has information from 2019 and from 2007, and must
%    allocate using those information. Stimulus, however, is given in 2020
%    and in 2008. So the planner needs to consider expected values in
%    consumption or welfare given the transition probabilities of states in
%    2007 to 2008 and in 2019 to 2020. The SNW_EVUVW19_JMKY file then
%    aggregates the full state-space results to just JMKY state-space,
%    which is the extend of information available to the planner.
%
%    Given 2020 JAEEMK (age, endogenous savings, education, income shock,
%    marital status, kids count), what is the expected value for the
%    planner given 2020 JAEEMK and transition between 2019 to 2020 JAEEMK
%    given some stimulus check assignment based on 2019 information?
%    (Stimulus amount set by WELF_CHECKS). This is similar to
%    SNW_EVUVW19_JAEEMK, except the solution here, under
%    SNW_EVUVW19_JAEEMK_FOC, relies on First Order Conditions, and are
%    hence faster.
%
%    Additionally, given 2008 policy and value functions (given expectation
%    of 2009 crisis unemployment shocks), call SNW_V08_JAEEMK to solve for
%    value and consumption given stimulus checks. And then integrate over
%    08 JAEEMK states given 07 JAEEMK states (age, endogenous savings,
%    education, income shock, marital status, kids count). The stimulus
%    will be provisioned based on 07 JAEEMK states. Note that
%    SNW_EVUVW19_JAEEMK does not solve the 07/08 problem.
%
%    [EV19_JAEEMK, EC19_JAEEMK, EV20_JAEEMK, EC20_JAEEMK] =
%    SNW_EVUVW19_JAEEMK_FOC(WELF_CHECKS, ST_SOLU_TYPE, MP_PARAMS,
%    MP_CONTROLS, V_SS_2020, AP_SS, CONS_SS_2020, V_UNEMP_2020,
%    CONS_UNEMP_2020, MP_PRECOMPUTE_RES) solves the 2019 biden or trump
%    check problem given 2020 (SS) policy and value and 2020 unemployment
%    policy value and consumtion. The SNWX_EVUVW19_JAEEMK_FOC vignette
%    tests these inputs.
%
%    [EV07_JAEEMK, EC07_JAEEMK, EV08_JAEEMK, EC08_JAEEMK] =
%    SNW_EVUVW19_JAEEMK_FOC(WELF_CHECKS, ST_SOLU_TYPE, MP_PARAMS,
%    MP_CONTROLS, AP_SS, V_2008, CONS_2008, MP_PRECOMPUTE_RES) solves the
%    2007 bush check problem given the steady-state savings function and
%    value and consutmpion in 2008. The 2008 Value and Consumption function
%    should be solved given 2009 expected unemployment shocks due to the
%    Great Recession.  the SNWX_EVUVW07_JAEEMK_FOC vignette tests these
%    inputs.
%
%    See also SNW_EVUVW19_JYMK, SNW_EVUVW19_JAEEMK_FOC,
%    SNWX_EVUVW19_JAEEMK_FOC, SNWX_EVUVW08_JAEEMK_FOC, SNW_EVUVW19_JAEEMK,
%    SNW_EVUVW20_JAEEMK, SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jaeemk_foc(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==5)
        [welf_checks, st_solu_type, mp_params, mp_controls, ...
            spt_mat_path] = varargin{:};
    elseif (length(varargin)==7)
        [welf_checks, mp_params, mp_controls, ...
            ap_ss, V_2008, cons_2008, ...
            mp_precompute_res] = varargin{:};
    elseif (length(varargin)==10)
        [welf_checks, st_solu_type, mp_params, mp_controls, ...
            V_ss_2020, ap_ss, cons_ss_2020, ...
            V_unemp_2020, cons_unemp_2020, ...
            mp_precompute_res] = varargin{:};
    else
        error('Need to provide 10 parameter inputs');
    end

else
    clc;
    close all;

    % The Number of Checks to Provide
    welf_checks = 2;

    % The Bush or the Biden/Trump Stimulus Problems
    st_biden_or_trump = 'bushchck';
%     st_biden_or_trump = 'bidenchk';

    % Printing Controls
    mp_controls = snw_mp_control('default_test');
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_precompute') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_v08p08_jaeemk') = false;
    mp_controls('bl_print_v08p08_jaeemk_verbose') = false;
    mp_controls('bl_print_v08_jaeemk') = true;
    mp_controls('bl_print_v08_jaeemk_verbose') = true;

    if (strcmp(st_biden_or_trump, 'bushchck'))
        %% Prepare inputs for the Bush Stimulus Problem

        % 1. generate MP_PARAMS specific to 2008 stimulus
        % Use non-default values for Bush Stimulus
        mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
        mp_more_inputs('fl_ss_non_college') = 0.225;
        mp_more_inputs('fl_ss_college') = 0.271;
        fl_p50_hh_income_07 = 54831;
        mp_more_inputs('fl_scaleconvertor') = fl_p50_hh_income_07;
        st_param_group = 'default_small';
    %     st_param_group = 'default_dense';
        mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
        mp_params('st_biden_or_trump') = st_biden_or_trump;

        % 2. Solve value steady state (2009 employed)
        [V_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
        V_emp_2009 = V_ss;

        % 3. Solve value unemployed 2009
        mp_params('xi') = 0.532;
        mp_params('b') = 0.37992;
        mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
        mp_params('TR') = 100/fl_p50_hh_income_07;
        [V_unemp_2009] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);

        % 4. Value and Optimal choice in 2008
        [V_2008, ap_2008, cons_2008, ev_empshk_2009] = ...
            snw_v08p08_jaeemk(mp_params, mp_controls, V_emp_2009, V_unemp_2009);

        % 5. matrixes to pre-compute
        % Only using the SNW_A4CHK_WRK_BISEC_VEC function, no unemployment
        % related matrixes needed Also don't need REF_EARN_WAGEIND_GRID,
        % become unemployment not conditional on wage in 2009.
        cl_st_precompute_list = {'a', ...
            'inc', 'inc_unemp', 'spouse_inc',...
            'ar_z_ctr_amz'};

    else
        %% Prepare inputs for the Biden/Trump Stimulus Problem

        % 1. generate MP_PARAMS specific to 2008 stimulus, set solution type
    %     st_solu_type = 'matlab_minimizer';
        st_solu_type = 'bisec_vec';
    %     st_solu_type = 'grid_search';

        % Solve the VFI Problem and get Value Function
        mp_params = snw_mp_param('default_tiny', false, 'tauchen', false, 8, 8);
    %     mp_params = snw_mp_param('default_dense');
    %     mp_params = snw_mp_param('default_moredense');

        % 2. Solve value steady state
        [V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);

        % 3a. 2020 V and C same as V_SS and cons_ss if tax the same, o/w resolve
        mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');
        if (mp_params('a2_covidyr') == mp_params('a2'))
            V_ss_2020 = V_ss;
            cons_ss_2020 = cons_ss;
        else
            % change xi and b to for people without unemployment shock
            % solving for employed but 2020 tax results
            mp_params('xi') = 1;
            mp_params('b') = 0;
            [V_ss_2020,~,cons_ss_2020,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
            mp_params('xi') = xi;
            mp_params('b') = b;
        end

        % 3b. 2020 V and C unemployed given MIT shock
        mp_params('xi') = 0.5;
        mp_params('b') = 0;
        fl_p50_hh_income_19 = 62502;
        mp_params('TR') = 100/fl_p50_hh_income_19;
        mp_params('st_biden_or_trump') = st_biden_or_trump;
        [V_unemp_2020,~,cons_unemp_2020,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);

        % 4. matrixes to pre-compute
        cl_st_precompute_list = {'a', ...
            'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid',...
            'ar_z_ctr_amz'};
    end

    % Shared: Steady-State distribution
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);

    % Shared: precompute, get Matrixes
    % note, the mp_params inputs are based on unemployed in 2020 (MIT) or unemployed in 2009 (Expected)
    % note, however, for the 2008/9 problem, only will use inc, inc_unemp, spouse_inc
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);

end

%% Parse Model Parameters
params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'bl_store_shock_trans', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, bl_store_shock_trans, cl_mt_pi_jem_kidseta, psi] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'st_biden_or_trump'});
[st_biden_or_trump] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_evuvw19_jaeemk', 'bl_print_evuvw19_jaeemk_verbose'});
[bl_print_evuvw19_jaeemk, bl_print_evuvw19_jaeemk_verbose] = params_group{:};

%% Solve evuvw19 Given Current Check
if (strcmp(st_biden_or_trump, 'bushchck'))
    if exist('spt_mat_path','var')
        [V_2008_check, C_2008_check] = snw_v08_jaeemk(...
            welf_checks, mp_params, mp_controls, spt_mat_path);
    else
        [V_2008_check, C_2008_check] = snw_v08_jaeemk(...
            welf_checks, ...
            mp_params, mp_controls, ...
            V_2008, cons_2008, ...
            mp_precompute_res);
    end
    ev20_jaeemk = V_2008_check;
    ec20_jaeemk = C_2008_check;
else
    if exist('spt_mat_path','var')
        [ev20_jaeemk, ec20_jaeemk] = snw_evuvw20_jaeemk(...
            welf_checks, ...
            st_solu_type, mp_params, mp_controls, ...
            spt_mat_path);
    else
        [ev20_jaeemk, ec20_jaeemk] = snw_evuvw20_jaeemk(...
            welf_checks, ...
            st_solu_type, mp_params, mp_controls, ...
            V_ss_2020, cons_ss_2020, ...
            V_unemp_2020, cons_unemp_2020, ...
            mp_precompute_res);
    end
end

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end


%% Parse Pre-Computes
if exist('spt_mat_path','var')
    load(spt_mat_path, 'ar_z_ctr_amz');
else
    params_group = values(mp_precompute_res, {'ar_z_ctr_amz'});
    [ar_z_ctr_amz] = params_group{:};
end

%% Planner val2019(J,A,E,E,M,K) = E(val2020(J+1,A',E',E,M,K')) given trans.

ev19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ec19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

if exist('spt_mat_path','var')
    load(spt_mat_path, 'ap_ss');
end

for j=1:n_jgrid-1 % Age

    % A1. mn_z_ctr and mn_aprime
    % array states/shocks all
    ar_aprime_amz = reshape(ap_ss(j,:,:,:,:,:), [], 1);

    % B1. Solve For EV(ap,z) = EV(ap,zp|z)f(zp|z) for all possible ap points
    mn_ev_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_ec_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    for educ=1:n_educgrid % Educational level
        for married=1:n_marriedgrid % Marital status

            % B1a. Get P(S'|S), S = [eta x kids] by [eta x kids] transition matrix
            if (bl_store_shock_trans)
                mt_pi_jem_kidseta = cl_mt_pi_jem_kidseta{j,educ,married};
            else
                mt_pi_jem_kidseta = kron(pi_kids(:,:,j,educ,married), pi_eta);
            end

            % B1b. Get age/edu/marry submatrix
            mn_ev20_jem = permute(ev20_jaeemk(j+1,:,:,educ,married,:), [2,3,6,1,4,5]);
            mn_ec20_jem = permute(ec20_jaeemk(j+1,:,:,educ,married,:), [2,3,6,1,4,5]);

            % B1c. 2D rows are savings states, columns are [eta x kids]
            mt_ev20_jem = reshape(mn_ev20_jem, n_agrid, []);
            mt_ec20_jem = reshape(mn_ec20_jem, n_agrid, []);

            % B1d. EV = V([a] by [eta x kids]) x TRANS([eta x kids] by [eta x kids])
            mt_ev19_jem = mt_ev20_jem*mt_pi_jem_kidseta';
            mt_ec19_jem = mt_ec20_jem*mt_pi_jem_kidseta';

            % B1e. Reshape Back and Store
            mn_ev19_jem = reshape(mt_ev19_jem, [n_agrid, n_etagrid, 1, 1, n_kidsgrid]);
            mn_ec19_jem = reshape(mt_ec19_jem, [n_agrid, n_etagrid, 1, 1, n_kidsgrid]);
            mn_ev_ap_z(:, :, educ, married, :) = mn_ev19_jem;
            mn_ec_ap_z(:, :, educ, married, :) = mn_ec19_jem;

        end
    end

    % C1. z specific EV Slope: EV(ap,z)/d(ap)
    mn_deri_dev_dap = diff(mn_ev_ap_z, 1)./diff(agrid);
    mn_deri_dec_dap = diff(mn_ec_ap_z, 1)./diff(agrid);
    % C2. ND dimensional Array to 2D dimensional Array:
    mt_ev_ap_z = reshape(mn_ev_ap_z, n_agrid, []);
    mt_ec_ap_z = reshape(mn_ec_ap_z, n_agrid, []);
    mt_deri_dev_dap = reshape(mn_deri_dev_dap, n_agrid-1, []);
    mt_deri_dec_dap = reshape(mn_deri_dec_dap, n_agrid-1, []);

    % D, Evaluate
    [ar_ev_aprime_z] = ffi_vec_v_ap(...
        ar_aprime_amz, agrid', ...
        ar_z_ctr_amz, ...
        mt_ev_ap_z, mt_deri_dev_dap);
    [ar_ec_aprime_z] = ffi_vec_v_ap(...
        ar_aprime_amz, agrid', ...
        ar_z_ctr_amz, ...
        mt_ec_ap_z, mt_deri_dec_dap);

    % G. Record Results
    mn_ev_aprime_z = reshape(ar_ev_aprime_z, size(mn_ev_ap_z));
    mn_ec_aprime_z = reshape(ar_ec_aprime_z, size(mn_ev_ap_z));

    % H. To main store:
    ev19_jaeemk(j,:,:,:,:,:) = mn_ev_aprime_z;
    ec19_jaeemk(j,:,:,:,:,:) = mn_ec_aprime_z;

end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_EVUVW19_JAEEMK_FOC", ...
         ['st_biden_or_trump=' char(mp_params('st_biden_or_trump'))], ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_evuvw19_jaeemk_verbose)

    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if (strcmp(st_biden_or_trump, 'bushchck'))
        st_ev_tyr = '07';
        st_ev_tp1 = '08';
    else
        st_ev_tyr = '19';
        st_ev_tp1 = '20';
    end
    % ev07_jaeemk or ev19_jaeemk
    mp_outcomes(['ev', st_ev_tyr, '_jaeemk']) = ev19_jaeemk;
    % ec07_jaeemk or ec19_jaeemk
    mp_outcomes(['ec', st_ev_tyr, '_jaeemk']) = ec19_jaeemk;
    % ev08_jaeemk or ev20_jaeemk
    mp_outcomes(['ev', st_ev_tp1, '_jaeemk']) = ev20_jaeemk;
    % ec08_jaeemk or ec20_jaeemk
    mp_outcomes(['ec', st_ev_tp1, '_jaeemk']) = ec20_jaeemk;
    ff_container_map_display(mp_outcomes, 9, 9);

end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev19_jaeemk;
    elseif (it_k==2)
        ob_out_cur = ec19_jaeemk;
    elseif (it_k==3)
        ob_out_cur = ev20_jaeemk;
    elseif (it_k==4)
        ob_out_cur = ec20_jaeemk;
    end
    varargout{it_k} = ob_out_cur;
end

end

% Utility given choices
function [ar_ev_aprime_z] = ffi_vec_v_ap(...
    ar_aprime, ar_a, ar_z_ctr_amz, ...
    mt_ev_ap_z, mt_deri_dev_dap)

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_it_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_it_ap_near_lower_idx(ar_it_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% the marginal effects of additional asset is determined by the slope
ar_deri_lin_idx = sub2ind(size(mt_deri_dev_dap), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_ev_lin_idx = sub2ind(size(mt_ev_ap_z), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_deri_dev_dap = mt_deri_dev_dap(ar_deri_lin_idx);
ar_ev_ap_lower_idx = mt_ev_ap_z(ar_ev_lin_idx);
clear ar_ev_lin_idx

% Ev(a_lower_idx,z) + slope*(fl_aprime - fl_a_lower)
ar_ev_aprime_z = ar_ev_ap_lower_idx + (ar_aprime - ar_a(ar_it_ap_near_lower_idx)').*ar_deri_dev_dap;
end
