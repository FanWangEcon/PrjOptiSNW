%% 2019 Age, Income, Kids, Marry EV and EC All Checks
% This is the example vignette for function:  <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/splanneralt/snw_evuvw20_jaeemk.m 
% *snw_evuvw20_jaeemk*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.* 2019 integrated over VU and VW
%% Test SNW_EVUVW19_JMKY_ALLCHECKS Parameters
% Save a result that is low in memory cost so that it can be loaded quickly 
% for various allocation tests. Turn off Various Printing Controls. Call function 
% with wide income bins to reduce memory storage and retrievel costs

clear all;
% Start mp contorls
mp_controls = snw_mp_control('default_test');
% Solve for Unemployment Values
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
mp_controls('bl_print_precompute') = false;
mp_controls('bl_print_precompute_verbose') = false;
mp_controls('bl_print_a4chk') = false;
mp_controls('bl_print_a4chk_verbose') = false;
mp_controls('bl_print_evuvw20_jaeemk') = false;
mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw19_jaeemk') = false;
mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw19_jmky') = false;
mp_controls('bl_print_evuvw19_jmky_verbose') = false;
%% 
% Dense default, and unemployment parameters:

% default dense load
mp_params = snw_mp_param('default_docdense');
% Unemployment
xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
TR=100/58056; % Value of a wezlfare check (can receive multiple checks). TO DO: Update with alternative values
mp_params('xi') = xi;
mp_params('b') = b;
mp_params('TR') = TR;
% Check Count: 89 checks to allow for both the first and the second round
n_welfchecksgrid = 3;
mp_params('n_welfchecksgrid') = n_welfchecksgrid;
%% 
% Income bins:

% Income Grid
% 4 refers to 4*58056=232224 dollars in 2012USD
% max 7 refers to 7*58056=406392 dollars in 2012USD
% all phase out = (4400/5)*100 + 150000 = 238000
% if 500 dollar interval, need 476 inc groups before 238000  
% if have 85 percent of points betwen 238000, 
fl_max_phaseout = 238000;
fl_multiple = 58056;
it_bin_dollar_before_phaseout = 5000;
it_bin_dollar_after_phaseout = 25000;
fl_thres = fl_max_phaseout/fl_multiple;
inc_grid1 = linspace(0,fl_thres,(fl_max_phaseout)/it_bin_dollar_before_phaseout);
inc_grid2 = linspace(fl_thres, 7, (7*fl_multiple-fl_max_phaseout)/it_bin_dollar_after_phaseout); 
inc_grid=sort(unique([inc_grid1 inc_grid2]'));
mp_params('n_incgrid') = length(inc_grid);
mp_params('inc_grid') = inc_grid;
%% SNW_EVUVW19_JMKY_ALLCHECKS Low Storage Invoke
% The simulation here (dense) requires less than 10 GB of memory with 8 workers 
% (8 threads needed), simulating over 88 checks takes with 8 workers 

st_solu_type = 'bisec_vec';
bl_parfor = false;
it_workers = 1;
bl_export = false; 
snm_suffix = ['_ybin' num2str(it_bin_dollar_before_phaseout)];
[ev19_jmky_allchecks, ec19_jmky_allchecks, output] = ...
    snw_evuvw19_jmky_allchecks(mp_params, mp_controls, ...
    st_solu_type, bl_parfor, it_workers, bl_export, snm_suffix);