%% Compute for Equilibrium Tax
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/fgov/snw_find_tax_rate.m 
% *snw_find_tax_rate*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*,> this function solves for equilibrium tax rate.
%% Parameter Controls

clear all;    
mp_params = snw_mp_param('default_docdense');
% mp_params = snw_mp_param('default_dense');
% mp_params = snw_mp_param('default_base');
% mp_params = snw_mp_param('default_small');
mp_params('beta') = 0.95;
xi=0.651; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)       
mp_params('xi') = xi;
mp_params('b') = b;
mp_controls = snw_mp_control('default_test');
%% 
% Parameters for COVID related Costs:

% Average check per household, given COVID actual policy payment schedule
% And given distribution. The number is from averaging over the actual
% allocations given distribution.
Covid_checks_per_capita = 18.7255856*100/62502;
% Covid_checks_per_capita = 0;
% which tax parameter to change a2 is the deafult, a0 shifts max tax rate
bl_adjust_a0 = false;
bl_load_existing = false;
%% 
% Graph Controls etc:

mp_controls('bl_timer') = true;
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
mp_controls('bl_print_find_tax_rate') = true;
mp_controls('bl_print_find_tax_rate_verbose') = true;
%% Solve for New Tax Rate
% Solve for Equilibrium Tax rate that clears government costs and income. In 
% the extreme bounding exercise, we assume the government will pay COVID costs 
% all in one year. This is to test whether an extreme tax scenario will lead to 
% changes in allocation results.
% 
% Given the checks that the government hands out and the taxes imposed, individual 
% resources post-tax are different in 2020. Households' savings decisions in 2020 
% vary with taxes and checks. However, the policy function post 2020 shifts back 
% to thte previous non-COVID world's policy function because the COVID shock is 
% an one period surprise shock.

a2 = snw_find_tax_rate(mp_params, mp_controls, Covid_checks_per_capita, bl_adjust_a0, bl_load_existing);
%% 
%% 
%%