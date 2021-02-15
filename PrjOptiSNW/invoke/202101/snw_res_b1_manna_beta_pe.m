%% Adjust Gamma and Beta
%{
Run comparative statics exercises with different beta values, but assuming
all consumers have the same beta. We do this in a "partial equilibrium"
manner where we keep all other model parameters, including endogenous
objects such as the income tax rate required to balance the government
budget, constant. We will then make some graphs with for example the
beta-value on the horizontal axis and the economy-wide average MPC (out of
a $1,200 check) on the vertical axis to understand how the average MPC
varies with beta.

Solve the model under the following beta-values: 
0.50, 0.60, 0.7, 0.8, 0.9, 0.95, 0.99. 
Then we can consider other values if 0.5 still isn't
sufficient to match the empirical MPC estimates.

see invoke/2020/snw_res_b0_xi0p25_manna.m for more detailed comments on
each section below.
%}

clc;
clear all;

% ls_fl_beta_val = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99];
% ls_fl_beta_val = [0.50, 0.60, 0.70, 0.80, 0.90];
ls_fl_beta_val = [0.90];

for fl_beta_val = ls_fl_beta_val   
    %% A1. Computing Specifications
    % 1a. Parfor controls
    bl_parfor = true;
    it_workers = 8;    
    % 1b. Export Controls
    bl_export = true;    
    % 1c. Solution Type
    st_solu_type = 'bisec_vec';    
    % 2. Set Up Parameters
    mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', false, 'tauchen', false, 8, 8);
    
    %% B1. Unemployment Shock and Benefits
    xi=0;
    b=1;  
    mp_params('xi') = xi;
    mp_params('b') = b;
    
    %% B2. Welfare Check Value And Numbers
    % The number of welfare checks to consider and the value of each checks
    TR=100/58056;
    n_welfchecksgrid = 89;    
    mp_params('TR') = TR;
    mp_params('n_welfchecksgrid') = n_welfchecksgrid;
    
    %% B3. Tax in 2020
    mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');
    
    %% C. Income Grid Solution Precision
    fl_max_phaseout = 238000;
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
    snw_evuvw19_jmky_allchecks(mp_params, mp_controls, st_solu_type, ...
        bl_parfor, it_workers, ...
        bl_export, bl_load_mat, snm_suffix);
    
    %% H. Stop Log
    diary off;
end