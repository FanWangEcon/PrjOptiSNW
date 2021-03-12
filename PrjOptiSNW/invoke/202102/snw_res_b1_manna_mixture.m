%% Solve for MPC At Each Check Level over Beta and Alpha Types
%{
Given calibration parameters from the mixture calibration from
snw_calibrate_tax_normgdp_multitypes.m, solve at 2 beta levels and 2
education levels separately, given the same tax/normalizing parameters. 

Solve for the MPC of Biden checks (3rd Round).

DO NOT run other files that will shift snw_hh_spousal_income.m if running
this file, otherwise will over-ride low vs high education income equation
specifications.

%}

clc;
clear all;

% The two mixture beta values
ls_fl_beta_val = [0.60, 0.95];
% Solve for HS and College for each type
ls_it_edu_simu_type = [1, 2];
% Biden check (3rd check)
st_biden_or_trump = 'bidenchk';
st_biden_or_trump = 'trumpchk';

% Solve for beta and edu combinations
for fl_beta_val = ls_fl_beta_val
    for it_edu_simu_type = ls_it_edu_simu_type

        %% A1. Computing Specifications
        % 1a. Parfor controls
        bl_parfor = true;
        it_workers = 8;
        % 1b. Export Controls
        bl_export = true;
        % 1c. Solution Type
        st_solu_type = 'bisec_vec';

        %% A2. Parameter Specifications by Education Groups.
        % 2. Set Up Parameters
        mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
        if (it_edu_simu_type == 1)
            mp_more_inputs('st_edu_simu_type') = 'low';
            st_param_group_suffix = '_e1lm2';
        elseif (it_edu_simu_type == 2)
            mp_more_inputs('st_edu_simu_type') = 'high';
            st_param_group_suffix = '_e2hm2';
        end
        
        % param group name
        % st_param_group_base = 'default_tiny';
        % st_param_group_base = 'default_small53';
        % st_param_group_base = 'default_dense';
        % st_param_group_base = 'default_docdense';
        % st_param_group_base = 'default_moredense_a65zh133zs5';
        st_param_group_base = 'default_moredense_a65zh266zs5';
        st_param_group = [st_param_group_base st_param_group_suffix];
        % get parametesr
        mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
        
        %% B1. Mixture Calibrated Parameters
        % see
        % Output202101/log/cali_tax_normgdp_multitypes_default_docdense.out
        mp_params('a2') = 1.6996;
        mp_params('a2_covidyr_manna_heaven') = 1.6996;
        mp_params('theta') = 0.5577;
                
        %% B2. Unemployment Shock and Benefits
        xi=0;
        b=1;
        mp_params('xi') = xi;
        mp_params('b') = b;
            
        %% B3. Welfare Check Value And Numbers
        % The number of welfare checks to consider and the value of each checks
        TR=100/58056;
        % Biden check: The current proposal by Congress is to give $1,400 to adults, $2,800 to married couples, and $1,400 per child
        n_welfchecksgrid = 169;
        mp_params('TR') = TR;
        mp_params('n_welfchecksgrid') = n_welfchecksgrid;

        %% B4. Tax in 2020
        mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');

        %% C. Income Grid Solution Precision
        % Let's use the current proposal for the third round by Congress.
        % The excel file shows how much single and married people with and
        % without children receive given income. It also shows when the
        % stimulus checks will phase out to zero (reaching 0 at an income
        % of $100k for singles without children; reaching 0 at an income of
        % $150k for singles with children (independent of how many
        % children); and reaching 0 at an income of $200k for married
        % couples (independent of how many children they have). In each
        % case, it's a linear slope from when it starts phasing out until
        % it reaches zero. The magnitude of the slope is a little bit
        % different compared to the common slope of 5 dollar stimulus
        % reduction for each 100 dollar income increase under the Trump
        % policy.
        % fl_max_phaseout = 238000;
%         fl_max_phaseout = 318000; % 7*58056 is 406392
        fl_max_phaseout = 200000;
        fl_multiple = 58056;
        it_bin_dollar_before_phaseout = 500;
        it_bin_dollar_after_phaseout = 5000;
        fl_thres = fl_max_phaseout/fl_multiple;
        inc_grid1 = linspace(0,fl_thres,(fl_max_phaseout)/it_bin_dollar_before_phaseout);
        inc_grid2 = linspace(fl_thres, 7, (7*fl_multiple-fl_max_phaseout)/it_bin_dollar_after_phaseout);
        inc_grid=sort(unique([inc_grid1 inc_grid2]'));

        mp_params('n_incgrid') = length(inc_grid);
        mp_params('inc_grid') = inc_grid;

        % 3. Controls
        mp_controls = snw_mp_control('default_test');

        %% D. Parameter overrid for gamma beta and r
        mp_params('beta') = fl_beta_val;
        st_file_name_suffix = ['_bt' num2str(round(fl_beta_val*100))];

        %% E. Display Control Parameters
        mp_controls('bl_print_vfi') = true;
        mp_controls('bl_print_vfi_verbose') = false;
        mp_controls('bl_print_ds') = true;
        mp_controls('bl_print_ds_verbose') = true;
        mp_controls('bl_print_precompute') = true;
        mp_controls('bl_print_precompute_verbose') = false;
        mp_controls('bl_print_a4chk') = false;
        mp_controls('bl_print_a4chk_verbose') = false;
        mp_controls('bl_print_evuvw20_jaeemk') = false;
        mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
        mp_controls('bl_print_evuvw19_jaeemk') = false;
        mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
        mp_controls('bl_print_evuvw19_jmky') = false;
        mp_controls('bl_print_evuvw19_jmky_verbose') = false;
        mp_controls('bl_print_evuvw19_jmky_mass') = false;
        mp_controls('bl_print_evuvw19_jmky_mass_verbose') = false;

        %% F1. Output Save nmae
        snm_suffix = ['_b1_xi0_manna_' num2str(n_welfchecksgrid-1) st_file_name_suffix];

        %% F2. Start log
        mp_paths = snw_mp_path('fan');
        spt_simu_outputs_log = mp_paths('spt_simu_outputs_log');

        snm_invoke_suffix = strrep(mp_params('mp_params_name'), 'default_', '');
        snm_file = ['snwx_v_planner_' char(snm_invoke_suffix) char(snm_suffix)];
        spn_log = fullfile(mp_paths('spt_simu_outputs_log'), [snm_file '.log']);

        diary(spn_log);

        %% F3. Log Show Parameters for Simulation
        ff_container_map_display(mp_params);
        ff_container_map_display(mp_controls);

        %% G. Run Checks Programs
        bl_load_mat = true;
        snw_evuvw19_jmky_allchecks(mp_params, mp_controls, ...
            st_biden_or_trump, st_solu_type, ...
            bl_parfor, it_workers, ...
            bl_export, bl_load_mat, snm_suffix);

        %% H. Stop Log
        diary off;
    end
end