%% Adjust Gamma and Beta
%{
Following AER insight comments, conduct exercises to vary beta, gamma and
r. Each time, solve for 89 check increments, which allows for double adults
as well as for doubling kids check amounts.
Based on December 21st chat, run these exercises, holding:

beta=0.971162552785405, gamma=2, and r=0.04, which were their original/current
values, solve and simulate 6 times:

1. gamma to 3
2. gamma to 1.05
3. beta to 0.95
4. beta to 0.99
5. r = 0.02
6. r = 0.01

see invoke/2020/snw_res_b0_xi0p25_manna.m for more detailed comments on
each section below.

for beta, 
%}

clc;
clear all;

% for it_run_type=[1,2,3,4,5,6]
for it_run_type=[3,4]
    %% A1. Beta, Gamma and R Adjustments
    gamma_override = NaN;
    rho_override = NaN;
    r_override = NaN;
    if (it_run_type == 1)
        gamma_override = 3;
    elseif (it_run_type == 2)
        gamma_override = 1.05;
    elseif (it_run_type == 3)
        rho_override = 0.90;
    elseif (it_run_type == 4)
        rho_override = 0.95;
    elseif (it_run_type == 5)
        r_override = 0.02;
    elseif (it_run_type == 6)
        r_override = 0.01;
    else
        error('it_run_type has to be 1 through 6');
    end
    
    %% A2. Computing Specifications
    % 1a. Parfor controls
    bl_parfor = true;
    it_workers = 8;
    % bl_parfor = false;
    % it_workers = 1;
    
    % 1b. Export Controls
    % bl_export = false;
    bl_export = true;
    
    % 1c. Solution Type
    st_solu_type = 'bisec_vec';
    
    % 2. Set Up Parameters
    % mp_params = snw_mp_param('default_tiny', false, 'tauchen', false, 8, 8);
    % mp_params = snw_mp_param('default_small', false, 'tauchen', false, 8, 8);
    % mp_params = snw_mp_param('default_base', false, 'tauchen', false, 8, 8);
    mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', false, 'tauchen', false, 8, 8);
    
    %% B1. Unemployment Shock and Benefits
    % xi=0, meaning that wage loss given covid shock is 100 percent
    % b=1, this is assuming there are full unemployment benefits
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
    st_file_name_suffix = ''
    if (~isnan(gamma_override))
        mp_params('gamma') = gamma_override;
        st_file_name_suffix = [st_file_name_suffix '_gm' num2str(round(gamma_override*100))];
    end
    if (~isnan(rho_override))
        mp_params('beta') = rho_override;
        st_file_name_suffix = [st_file_name_suffix '_bt' num2str(round(rho_override*100))];
    end
    if (~isnan(r_override))
        mp_params('r') = r_override;
        st_file_name_suffix = [st_file_name_suffix '_rr' num2str(round(r_override*100))];
    end
    
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