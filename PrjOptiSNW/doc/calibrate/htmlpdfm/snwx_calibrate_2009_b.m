%% UI Benefit Unemployment Lost Wage Recovery Parameter b Calibration
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/calibrate/snw_calibrate_2009_b.m 
% *snw_calibrate_2009_b*> from the <https://fanwangecon.github.io/PrjOptiSNW/ 
% *PrjOptiSNW Package*>*.* 
% 
% The ratio of UI benefits to wages and salary is 2.1 percent in 2009. $\xi\in[0,1]$ 
% governs the duration of unemployment shock for those unemployed. This equals 
% to 0.532 in 2009 ($\xi$ = 0 no wages earned).
% 
% We solve for total wage earnings from unemployed and employed in 2009, for 
% employed, same as under steady-state. For unemployed, they lose ($1 - \xi$) 
% share of the wage they would otherwise have earned. Our unemployment probability 
% in 2009 is conditional on age and edu groups (SNW_UNEMP_2008.m) computed based 
% on rectiilnear restriction. 
% 
% We know total UI amount (multiply its share of total "Wages and salary" by 
% total "wages and salary". We know how much wage was lost due to $\xi$. The ratio 
% of these two levels is b, which is the parameter that is the share of lost-wage 
% recovered. Note that this is based on exogenous wage earnings, so we do not 
% have to worry about endogenous changes to savings. We will solve for the steady-state 
% distribution, which generates mass of people by age, education, marital status, 
% kids count, etc.
%% Calibrate b with 2.1% UI Benefits to Wages Ratio and $\xi=0.532$
% Using various default parameters, including the default unemployment in 2009 
% matrix, and the default $\xi=0.532$ parameter, compute b.

% Solve parameters
mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
mp_more_inputs('fl_ss_non_college') = 0.225;
mp_more_inputs('fl_ss_college') = 0.271;
mp_more_inputs('fl_scaleconvertor') = 54831;
% st_param_group = 'default_small';
st_param_group = 'default_docdense';
mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
% Controls
mp_controls = snw_mp_control('default_test');

% no b, solving for b, b set to 0 when solving for wages
xi=0.532; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
mp_params('xi') = xi;

% Solve for Unemployment Values
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
mp_controls('bl_print_calibrate_2009') = true;
mp_controls('bl_print_calibrate_2009_verbose') = false;

% 2.1% UI Benefits to Wages and Salary Ratio
fl_ratio_ui_benefits_to_wage = 0.021;

% Solve
[fl_b_calibrated_by_ui_share, ...
    mp_stats_wage_ui_spending, ...
    mn_earn_tot_wgted, mn_earn_unemp_wgted, ...
    mn_earn_unemp_tot_wgted, mn_earn_unemp_weighted_wgted] = ...
    snw_calibrate_2009_b(mp_params, mp_controls, ...
    fl_ratio_ui_benefits_to_wage);
%% Calibrate b with 5.68% UI Benefits to Wages Ratio and $\xi=0.651$
% Change the benefit share and $\xi$ parameter to COVID values. The $b$ we find 
% below is not what should be used for COVID, the unemployment probability is 
% based on 2009 crisis still. That is hard-coded into the   <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/calibrate/snw_calibrate_2009_b.m 
% *snw_calibrate_2009_b*> function.

% Solve parameters
mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
mp_more_inputs('fl_ss_non_college') = 0.225;
mp_more_inputs('fl_ss_college') = 0.271;
mp_more_inputs('fl_scaleconvertor') = 54831;
% st_param_group = 'default_small';
st_param_group = 'default_dense';
mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
% Controls
mp_controls = snw_mp_control('default_test');

% no b, solving for b, b set to 0 when solving for wages
xi=0.651; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
mp_params('xi') = xi;

% Solve for Unemployment Values
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
mp_controls('bl_print_calibrate_2009') = true;
mp_controls('bl_print_calibrate_2009_verbose') = false;

% 2.1% UI Benefits to Wages and Salary Ratio
fl_ratio_ui_benefits_to_wage = 0.0568;

% Solve
[fl_b_calibrated_by_ui_share, ...
    mp_stats_wage_ui_spending, ...
    mn_earn_tot_wgted, mn_earn_unemp_wgted, ...
    mn_earn_unemp_tot_wgted, mn_earn_unemp_weighted_wgted] = ...
    snw_calibrate_2009_b(mp_params, mp_controls, ...
    fl_ratio_ui_benefits_to_wage);