%% Distribution Exact Savings Choices Vectorized
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_vec.m 
% *snw_ds_main_vec*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.* This function solves for vfi and gets distribution induced by policy 
% functions and exogenous distributions.More Dense Simulation.  Vectorized to 
% get distribution, but uses *bisect vec* for VFI.
%% Test SNW_DS_MAIN_VEC Defaults Dense
% Call the function with testing defaults.

mp_params = snw_mp_param('default_docdense');
mp_controls = snw_mp_control('default_test');
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_ds') = true;
mp_controls('bl_print_ds_verbose') = false;
[Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] = snw_ds_main_vec(mp_params, mp_controls);
% [Phi_true,Phi_adj] = snw_ds_main(mp_params, mp_controls);
Phi_true = Phi_true/sum(Phi_true(:));
%% Show All Info in mp_dsvfi_results More Dense

mp_cl_mt_xyz_of_s = mp_dsvfi_results('mp_cl_mt_xyz_of_s');
disp(mp_cl_mt_xyz_of_s('tb_outcomes'))