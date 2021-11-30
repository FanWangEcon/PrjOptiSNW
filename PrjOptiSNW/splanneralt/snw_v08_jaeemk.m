%% SNW_V08_JAEEMK Solve for V and C in 2008 given Stimulus Checks
%    This is simialr to SNW_EVUVW20_JAEEMK, but for the 2008 Bush stimulus.
%    SNW_V08P08_JAEEMK already solved for optimal policy and value functions in 2008,
%    given expected unemployment shock in 2009.
%
%    In this function, given some stimulus amount, we use
%    SNW_A4CHK_WRK_BISEC_VEC to compute the updated optimal V and C in 2008
%    given the stimulus amount, based on the alues for V and C without
%    stimulus computed by SNW_V08P08_JAEEMK.
%
%    Note that SNW_A4CHK_WRK_BISEC_VEC computes the adjustment in the
%    savings state that would be equivalent to the increase in stimulus
%    amount (which is not a state variable) to current resources, this is
%    faster than resolving 2008 optimal V and C at specific stimulus check
%    amount levels.
%
%    Note SNW_EVUVW20_JAEEMK has EVUVW, but here, we only have V08, because
%    in the 2020 problem, households receive checks ex-post of the COVID
%    MIT shocks in 2008 and the EVUVW is the weighted average in V between
%    the MIT unemployed and non-shock employed state. In 2008, however,
%    there are no shocks yet. The Bush stimulus is provided ex-ante of the
%    shock realization. The 2009 shocks due to the great recession is not a
%    MIT shock, but expected shock. The effect of the 2009 shock on
%    consumption, savings is solved by SNW_V08P08_JAEEMK. The expectation
%    over shock, in another word, for the SNW_V08_JAEEMK is already
%    included in EV' in 2008 for 2009.
%
%    [V_2008_CHECK, C_2008_CHECK] = SNW_V08_JAEEMK(WELF_CHECKS,
%    MP_PARAMS, MP_CONTROLS, V_2008, CONS_2008, MP_PRECOMPUTE_RES) provide
%    V_SS and V_UNEMP solved out elsewhere, and only get two outputs out.
%
%    See also SNW_V08_JAEEMK, SNWX_EVUVW08_JAEEMK, SNW_HH_PRECOMPUTE,
%    SNW_A4CHK_WRK_BISEC_VEC, SNW_EVUVW20_JAEEMK
%

%%
function [varargout]=snw_v08_jaeemk(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==4)
        [welf_checks, mp_params, mp_controls, spt_mat_path] = varargin{:};
    elseif (length(varargin)==6)
        [welf_checks, mp_params, mp_controls, V_2008, cons_2008, mp_precompute_res] = varargin{:};
    else
        error('Need to provide 6 parameter inputs');
    end

else
    clc;
    close all;

    % 1. Paramters
    % The Number of Checks to Provide
    welf_checks = 2;

    mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
    mp_more_inputs('fl_ss_non_college') = 0.225;
    mp_more_inputs('fl_ss_college') = 0.271;
    fl_p50_hh_income_07 = 54831;
    mp_more_inputs('fl_scaleconvertor') = fl_p50_hh_income_07;
    st_param_group = 'default_small';
%     st_param_group = 'default_dense';
    mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);

    mp_controls = snw_mp_control('default_test');
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_v08p08_jaeemk') = false;
    mp_controls('bl_print_v08p08_jaeemk_verbose') = false;
    mp_controls('bl_print_v08_jaeemk') = true;
    mp_controls('bl_print_v08_jaeemk_verbose') = true;

    % 2. Solve value steady state (2009 employed)
    [V_VFI_ss, ap_VFI_ss, cons_VFI_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    V_emp_2009 = V_VFI_ss;
    % Solve for probability mass, needed for pre-compute
    [Phi_true] = snw_ds_main_vec(mp_params, mp_controls, ap_VFI_ss, cons_VFI_ss, mp_valpol_more_ss);

    % 3. Solve value unemployed 2009
    % Set Unemployment Related Variables
    mp_params('xi') = 0.532;
    mp_params('b') = 0.37992;
    mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
    mp_params('TR') = 100/fl_p50_hh_income_07; % Value of a stimulus check (can receive multiple checks). TO DO: Update with alternative values;
    [V_unemp_2009] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_VFI_ss);

    % 4. Value and Optimal choice in 2008
    [V_2008, ap_2008, cons_2008, ev_empshk_2009] = ...
        snw_v08p08_jaeemk(mp_params, mp_controls, V_emp_2009, V_unemp_2009);

    % 5. pre-compute
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid'};
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_VFI_ss, Phi_true);

end

%% Parse Pre-Computes
if ~exist('spt_mat_path','var')
    params_group = values(mp_precompute_res, {'ar_a_amz'});
    [ar_a_amz] = params_group{:};
    params_group = values(mp_precompute_res, {'ar_inc_amz'});
    [ar_inc_amz] = params_group{:};
    params_group = values(mp_precompute_res, {'ar_spouse_inc_amz'});
    [ar_spouse_inc_amz] = params_group{:};
end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp'});
[pi_unemp] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_v08_jaeemk', 'bl_print_v08_jaeemk_verbose'});
[bl_print_v08_jaeemk, bl_print_v08_jaeemk_verbose] = params_group{:};

%% Compute Check V_2008_check and check V_U Values
if (bl_print_v08_jaeemk)
    disp(['Solve for V_2008_check for ', num2str(welf_checks), ' stimulus checks'])
end

if exist('spt_mat_path','var')
    [V_2008_check,C_2008_check]=snw_a4chk_wrk_bisec_vec(welf_checks, mp_params, mp_controls, spt_mat_path);

else
    [V_2008_check,C_2008_check]=snw_a4chk_wrk_bisec_vec( ...
        welf_checks, V_2008, cons_2008, mp_params, mp_controls, ...
        ar_a_amz, ar_inc_amz, ar_spouse_inc_amz);
    clear V_2008 cons_2008 ar_inc_amz ar_spouse_inc_amz

end

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_V08_JAEEMK", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['timeEUEC=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_v08_jaeemk_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('V_2008_check') = V_2008_check;
    mp_outcomes('C_2008_check') = C_2008_check;
    ff_container_map_display(mp_outcomes, 9, 9);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_2008_check;
    elseif (it_k==2)
        ob_out_cur = C_2008_check;
    end
    varargout{it_k} = ob_out_cur;
end

end
