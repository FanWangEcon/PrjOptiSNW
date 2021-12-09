%% 2007 Age, Income, Kids, Marry EV and EC All Checks (Bush Checks)
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw19_jmky_allchecks.m 
% *snw_evuvw19_jmky_allchecks*> from the <https://fanwangecon.github.io/PrjOptiSNW/ 
% *PrjOptiSNW Package*>*.* 2019 integrated over VU and VW
% 
% The function <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw19_jmky_allchecks.m 
% *snw_evuvw19_jmky_allchecks*> was initiall designed to handle the COVID problem, 
% the revised version of the program handles both the 2007/8/9 Bush stimulus check 
% problem, and the 2019/20/21 Trump and Biden stimulus check problems.
% 
% The key features of the Bush stimulus checks are: i) determined based on 2007 
% information; ii) checks received in 2008, when the great recession has not arrived 
% yet, but all expect it to in 2009; iii) the Great Recession hits in 2009, putting 
% some people, based on education and age, into unemployment state with a shared 
% unemployment duration and lost income and also UI benefits (calibrated to match 
% overall UI share of wages); iv) the economy returns to steady-state in 2010.
%% Test SNW_EVUVW19_JMKY_ALLCHECKS Parameters for Bush Checks
% Save a result that is low in memory cost so that it can be loaded quickly 
% for various allocation tests. Turn off Various Printing Controls. Call function 
% with wide income bins to reduce memory storage and retrievel costs

clear all;
% Start mp contorls
mp_controls = snw_mp_control('default_test');
% Solve for Unemployment Values
mp_controls('bl_timer') = true;
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = true;
mp_controls('bl_print_precompute') = false;
mp_controls('bl_print_precompute_verbose') = false;
mp_controls('bl_print_a4chk') = false;
mp_controls('bl_print_a4chk_verbose') = false;
mp_controls('bl_print_v08p08_jaeemk') = false;
mp_controls('bl_print_v08p08_jaeemk_verbose') = false;
mp_controls('bl_print_v08_jaeemk') = false;
mp_controls('bl_print_v08_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw20_jaeemk') = false;
mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw19_jaeemk') = false;
mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw19_jmky') = false;
mp_controls('bl_print_evuvw19_jmky_verbose') = false;
%% 
% Dense default, and unemployment parameters:

% default dense load
% 1. generate MP_PARAMS specific to 2008 stimulus
% Use non-default values for Bush Stimulus
mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
mp_more_inputs('fl_ss_non_college') = 0.225;
mp_more_inputs('fl_ss_college') = 0.271;
fl_p50_hh_income_07 = 54831;
mp_more_inputs('fl_scaleconvertor') = fl_p50_hh_income_07;
% st_param_group = 'default_small';
% st_param_group = 'default_dense';
st_param_group = 'default_docdense';
mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
% mp_params = snw_mp_param('default_docdense')
mp_params('beta') = 0.95;
fl_scaleconvertor = 62502;
mp_more_inputs('fl_scaleconvertor') = fl_scaleconvertor;
% Unemployment
mp_params('xi') = 0.532;
mp_params('b') = 0.37992;
mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
mp_params('TR') = 100/fl_p50_hh_income_07;
% Check Count: 89 checks to allow for both the first and the second round
n_welfchecksgrid = 3;
mp_params('n_welfchecksgrid') = n_welfchecksgrid;
mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
%% 
% Income bins:

% Income Grid
% 4 refers to 4*58056=232224 dollars in 2012USD
% max 7 refers to 7*58056=406392 dollars in 2012USD
% all phase out = (4400/5)*100 + 150000 = 238000
% if 500 dollar interval, need 476 inc groups before 238000  
% if have 85 percent of points betwen 238000, 
fl_max_phaseout = 238000;
fl_multiple = fl_scaleconvertor;
it_bin_dollar_before_phaseout = 5000;
it_bin_dollar_after_phaseout = 25000;
fl_thres = fl_max_phaseout/fl_multiple;
inc_grid1 = linspace(0,fl_thres,(fl_max_phaseout)/it_bin_dollar_before_phaseout);
inc_grid2 = linspace(fl_thres, 7, (7*fl_multiple-fl_max_phaseout)/it_bin_dollar_after_phaseout); 
inc_grid=sort(unique([inc_grid1 inc_grid2]'));
mp_params('n_incgrid') = length(inc_grid);
mp_params('inc_grid') = inc_grid;
%% SNW_EVUVW19_JMKY_ALLCHECKS Low Storage Invoke for Bush Checks
% The simulation here (dense) requires less than 10 GB of memory with 8 workers 
% (8 threads needed), simulating over 88 checks takes with 8 workers 

st_biden_or_trump = 'bushchck';
st_solu_type = 'bisec_vec';
bl_parfor = false;
it_workers = 1;
bl_export = false; 
bl_load_mat = false;
snm_suffix = ['_test_ybin' num2str(it_bin_dollar_before_phaseout)];
[ev19_jmky_allchecks, ec19_jmky_allchecks, output] = ...
    snw_evuvw19_jmky_allchecks(mp_params, mp_controls, ...
    st_biden_or_trump, st_solu_type, ...
    bl_parfor, it_workers, ...
    bl_export, bl_load_mat, snm_suffix);