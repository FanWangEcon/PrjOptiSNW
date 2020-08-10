%% Model Parameters
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/params/snw_mp_param.m 
% *snw_mp_param*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.* This function sets and gets different parameters
%% Documentation Run Parameters Docdense
% Parameters used for documentation vig.

mp_params = snw_mp_param('default_docdense', true, 100, 6);
%% Parameters Used for Test Simulation
% Rather than solving for all ages between 18 to 100, this solves for age groups, 
% and has limited shocks and asset levels.

mp_params = snw_mp_param('default_small', true, 100, 6);
%% Parameters Used for Paper Simulations
% Using 266 household head income shocks. Requires 150GB of memory.

% mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', true, 100, 6);
%%