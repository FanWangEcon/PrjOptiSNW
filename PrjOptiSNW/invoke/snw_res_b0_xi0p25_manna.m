%% Zero Unemployment Benefits
% Mana from Heaven (no COVID check Tax Adjustments

clc;
clear all;

%% A. Computing Specifications
% On precision, with 12 cores and 190 GB of memory can run it_worker = 12
% for 245 checks, including initial distributional computing time etc
% should take between 18 hours to 24 hours. If Distribution and VFI already
% stored, takes less time closer to the lower bound. If not, requires up to
% closer to 24 hours.

% Store or not, parallel or not, parameters to use
% 1a. Parfor controls
bl_parfor = true;
it_workers = 12;
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
mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', false, 'tauchen', false, 8, 8);

%% B1. Unemployment Shock and Benefits
% xi=0.25, meaning that wage loss given covid shock is 75 percent
% b=0, this is assuming there are no unemployment benefits

xi=0.25;
b=0;

mp_params('xi') = xi;
mp_params('b') = b;

%% B2. Welfare Check Value And Numbers
% The number of welfare checks to consider and the value of each checks

TR=100/58056;
n_welfchecksgrid = 245;

mp_params('TR') = TR;
mp_params('n_welfchecksgrid') = n_welfchecksgrid;

%% B3. Tax in 2020
% Can either have the same tax in all years, tax calibrated to pre-covid
% level. Or use COVID specific tax. In particular, pay for all covid check
% costs in one year (2020) by adjusting one year tax rate. Found via
% SNW_FIND_TAX_RATE

mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');

%% C. Income Grid Solution Precision
% How much 2019 income groups to generate
% 500 dollar precision before full threshold, 5000 dollar pecision
% afterwards

% 4 refers to 4*58056=232224 dollars in 2012USD
% max 7 refers to 7*58056=406392 dollars in 2012USD
% all phase out = (4400/5)*100 + 150000 = 238000
% if 500 dollar interval, need 476 inc groups before 238000
% if have 85 percent of points betwen 238000,

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

%% D. Display Control Parameters
% 6. Display Controls
% Solve for Unemployment Values
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

%% E1. Output Save nmae
% output will be saved with the name of the param group, plus suffix
% indicating tax, b, xi parameters that can be specified via snw_suffix

snm_suffix = ['_b0_xi0p25_manna_' num2str(n_welfchecksgrid-1)];

%% E2. Start log
mp_paths = snw_mp_path('fan');
spt_simu_outputs_log = mp_paths('spt_simu_outputs_log');

snm_invoke_suffix = strrep(mp_params('mp_params_name'), 'default_', '');
snm_file = ['snwx_v_planner_' char(snm_invoke_suffix) char(snm_suffix)];
spn_log = fullfile(mp_paths('spt_simu_outputs_log'), [snm_file '.log']);

diary(spn_log);

%% E3. Log Show Parameters for Simulation

ff_container_map_display(mp_params);
ff_container_map_display(mp_controls);

%% F. Run Checks Programs
bl_load_mat = true;
snw_evuvw19_jmky_allchecks(mp_params, mp_controls, st_solu_type, ...
    bl_parfor, it_workers, ...
    bl_export, bl_load_mat, snm_suffix);

%% G. Stop Log

diary off;