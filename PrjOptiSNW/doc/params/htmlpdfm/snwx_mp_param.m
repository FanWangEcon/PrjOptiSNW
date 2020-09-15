%% Model Parameters
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/params/snw_mp_param.m 
% *snw_mp_param*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.* This function sets and gets different parameters.
%% Parameters Used for Test Simulation
% Rather than solving for all ages between 18 to 100, this solves for age groups, 
% and has limited shocks and asset levels. Used for testing.

mp_params = snw_mp_param('default_small', true, 100, 6);
%% Documentation Run Parameters Docdense
% Parameters used for documentation vig. "docdense" uses less shocks than the 
% version of the model used to implement the allocation problems in the <https://cepr.org/active/publications/discussion_papers/dp.php?dpno=15283 
% Nygaard, Sorensen and Wang (2020)>.

mp_params = snw_mp_param('default_docdense', true, 100, 6);
%% Parameters Used for Paper Simulations
% Full version of parameters used in <https://cepr.org/active/publications/discussion_papers/dp.php?dpno=15283 
% Nygaard, Sorensen and Wang (2020)>. This is not printed to save space.

% mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', true, 100, 6);
%%