%% Distribution with One Period Policy Shift
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.*  
%% One-period Deviation from Steady-State given Alternative Policy Function
% In addition to solving for distribution given one policy function, <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*> can also solve for the distributional shift from "steady-state" 
% with a one-period policy shift. 
% 
% If a 6th parameter, PHI_ADJ_BASE, is provided to <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*>*,* solve for next-period forward distribution conditional 
% on PHI_ADJ_BASE, using the policy function provided to <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*> as the 3rd and 4th parameters.
% 
% When PHI_ADJ_BASE is provided, if the AP_SS, CONS_SS policy functions inputs 
% are from the same problem that generated PHI_ADJ_BASE, output PHI_ADJ will be 
% identical to PHI_ADJ_BASE. However, if AP_SS, CONS_SS are different policy functions 
% from those that induced PHI_ADJ_BASE,  PHI_ADJ output will be different from 
% PHI_ADJ_BASE input. 
% 
% This allows for obtaining the distributional impact of a one period policy, 
% allowing for deviation from "steady-state" distribution. This is used to solve 
% for the distribution after one-period MIT shock, given stimulus checks provided 
% in that period.
% 
% This is used to model the distributional effects of CARES Act, the two rounds 
% of Trump Stimulus Checks, on household asset distribution when then receive 
% the Biden stimulus checks from the the American Recovery Act. In effect, we 
% have two MIT shock periods.
%% Solve for "Steady-State" Policy and Value Functions
% Steady-state policy and value functions

% mp_params = snw_mp_param('default_dense');
mp_params = snw_mp_param('default_docdense');
% mp_params = snw_mp_param('default_moredense_a65zh133zs5_e2m2');
mp_controls = snw_mp_control('default_test');
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
[V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
%% Solve for "Steady-State" Distribution
% Solve for steady-state distributions, using steady-state policy functions.

[Phi_true_ss,Phi_adj_ss,A_agg_ss,Y_inc_agg_ss,~,mp_dsvfi_results_ss] = ...
    snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss);
% [Phi_true,Phi_adj] = snw_ds_main(mp_params, mp_controls);
Phi_true_ss = Phi_true_ss/sum(Phi_true_ss(:));
%% 
% Show distributional results.

mp_cl_mt_xyz_of_s = mp_dsvfi_results_ss('mp_cl_mt_xyz_of_s');
disp(mp_cl_mt_xyz_of_s('tb_outcomes'));
%% Solve for Policy Function Under Trump Stimulus
% Same continuation value as prior (steady-state continuation), but now solve 
% for new policy (one round) due to Trump stimulus. Same tax rate in covid and 
% other years, manna-from-heaven. This calls the <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec_stimulus.m 
% *snw_vfi_main_bisec_vec_stimulus*> function, which provides the stimulus checks 
% as a function of income and family status.

mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');
[~,ap_trumpchecks,cons_trumpchecks, mp_valpol_more_trumpchecks] = ...
    snw_vfi_main_bisec_vec_stimulus(mp_params, mp_controls, V_ss);
%% Solve for Updated Distribution given Trump Stimulus 
% Fixing mass at their steady-state distribution, policy functions shift to 
% the Trump stimulus policies, resolve for one-period forward distribution. The 
% distributional code is almost identical, except uses steady-state distribution 
% as the "base" distribution via parameter PHI_ADJ_SS.

[Phi_true_trumpchecks,Phi_adj_trumpchecks,...
    A_agg_trumpchecks,Y_inc_agg_trumpchecks,~,mp_dsvfi_results_trumpchecks] = snw_ds_main_vec(...
    mp_params, mp_controls, ...
    ap_trumpchecks, cons_trumpchecks, ...
    mp_valpol_more_trumpchecks, ...
    Phi_adj_ss);
Phi_true_trumpchecks = Phi_true_trumpchecks/sum(Phi_true_trumpchecks(:));
%% 
% Show distributional results.

mp_cl_mt_xyz_of_s = mp_dsvfi_results_trumpchecks('mp_cl_mt_xyz_of_s');
disp(mp_cl_mt_xyz_of_s('tb_outcomes'));
%% Debug Check, SNW_DS_MAIN_VEC with Steady State Policies
% This is to confirm that code is working properly. If we use steady-state policy 
% functions and also provide as a sixth parameter the steady-state distribution, 
% PHI_ADJ_SS, to <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*>, we should get back the same distribution, PHI_TRUE_SS_WITH_EXISTDIST_DEBUG, 
% which is the same as PHI_ADJ_SS. See that the distributional outputs at the 
% end of this subsection is the same as the distributional table before the table 
% directly prior.

[Phi_true_ss_with_existdist_debug,~,~,~,~,mp_dsvfi_results_ss_with_existdist_debug] = snw_ds_main_vec(...
    mp_params, mp_controls, ...
    ap_ss, cons_ss, ...
    mp_valpol_more_ss, ...
    Phi_adj_ss);
Phi_true_ss_with_existdist_debug = Phi_true_ss_with_existdist_debug/sum(Phi_true_ss_with_existdist_debug(:));
%% 
% Show distributional results.

mp_cl_mt_xyz_of_s = mp_dsvfi_results_ss_with_existdist_debug('mp_cl_mt_xyz_of_s');
disp(mp_cl_mt_xyz_of_s('tb_outcomes'));