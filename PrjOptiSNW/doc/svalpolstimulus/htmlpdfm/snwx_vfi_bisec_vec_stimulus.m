%% Life Cycle Dynamic Programming under with CARES Act Stimulus Checks
% This is the example vignette for function:  <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec_stimulus.m 
% *snw_vfi_main_bisec_vec_stimulus*> from the <https://fanwangecon.github.io/PrjOptiSNW/ 
% *PrjOptiSNW Package*>*.* This function solves for policy function using Exact 
% Vectorized Solution. Value in 2020 with surprise COVID unemployment Shock, with 
% non-covid year Value as the continuation function, and provides households with 
% stimulus checks specified in the 1st and 2nd round under actual Trump admin 
% policies. The file focuses on the change in value function, asset choice, and 
% consumption choice given a one period unemployment shock (that does not reappear 
% in the future again). Solving this provides the distribution needed for the 
% Biden checks, American Rescue Plan, problem.
%% Test SNW_VFI_MAIN_BISEC_VEC_STIMULUS
% Solve the Regular Value and Also the Unemployment Value.
% 
% First, solve for value without unemployment issue (use the vectorized code 
% that was previously tested):

mp_params = snw_mp_param('default_docdense');
mp_controls = snw_mp_control('default_test');
[V_VFI_ss,ap_VFI_ss,cons_VFI_ss,mp_valpol_more_ss] = ...
    snw_vfi_main_bisec_vec(mp_params, mp_controls);
%% 
% Second, solve for the unemployment value, use the exact-bisec result code, 
% call the snw_vfi_main_bisec_vec.m function with a third input of existing value. 
% xi is the share of income lost during covid year given surprise covid shock, 
% b is the share of income loss that is covered by unemployment insurance. xi=0.5 
% and b=0 means will lose 50 percent of income given COVID shocks, and the loss 
% will not be covered at all by unemployment insurance. Calling the <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec_stimulus.m 
% *snw_vfi_main_bisec_vec_stimulus*> means households will receive positive amounts 
% of stimulus given household structure (marital status and children count), as 
% well as their total household income level. 

mp_params('xi') = 0.5;
mp_params('b') = 0;
mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');
[V_VFI_wthtrumpchk,ap_VFI_wthtrumpchecks,cons_VFI_wthtrumpchk,mp_valpol_more_wthtrumpchk] = ...
    snw_vfi_main_bisec_vec_stimulus(mp_params, mp_controls, V_VFI_ss);
%% 
% Difference Between Value and Choices In Unemployment and Future Periods

V_VFI_wthtrumpchk_drop = V_VFI_ss - V_VFI_wthtrumpchk;
ap_VFI_wthtrumpchk_drop = ap_VFI_ss - ap_VFI_wthtrumpchecks;
cons_VFI_wthtrumpchk_drop = cons_VFI_ss - cons_VFI_wthtrumpchk;
%% Define Parameter Frames
% Define the matrix dimensions names and dimension vector values. Policy and 
% Value Functions share the same ND dimensional structure.

% Grids:
age_grid = 18:100;
agrid = mp_params('agrid')';
eta_H_grid = mp_params('eta_H_grid')';
eta_S_grid = mp_params('eta_S_grid')';
ar_st_eta_HS_grid = string(cellstr([num2str(eta_H_grid', 'hz=%3.2f;'), num2str(eta_S_grid', 'wz=%3.2f')]));
edu_grid = [0,1];
marry_grid = [0,1];
kids_grid = (1:1:mp_params('n_kidsgrid'))';
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cl_mp_datasetdesc = {};
cl_mp_datasetdesc{1} = containers.Map({'name', 'labval'}, {'age', age_grid});
cl_mp_datasetdesc{2} = containers.Map({'name', 'labval'}, {'savings', agrid});
cl_mp_datasetdesc{3} = containers.Map({'name', 'labval'}, {'eta', 1:length(eta_H_grid)});
cl_mp_datasetdesc{4} = containers.Map({'name', 'labval'}, {'edu', edu_grid});
cl_mp_datasetdesc{5} = containers.Map({'name', 'labval'}, {'marry', marry_grid});
cl_mp_datasetdesc{6} = containers.Map({'name', 'labval'}, {'kids', kids_grid});
%% Analyze Savings and Shocks
% First, analyze Savings Levels and Shocks, Aggregate Over All Others, and do 
% various other calculations.

% Generate some Data
mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_support_graph('cl_st_xtitle') = {'Savings States, a'};
mp_support_graph('st_legend_loc') = 'eastoutside';
mp_support_graph('bl_graph_logy') = true; % do not log
mp_support_graph('it_legend_select') = 15; % how many shock legends to show
mp_support_graph('cl_colors') = 'jet';
%% 
% MEAN(VAL(A,Z) - VAL(A,Z|CARESActChecks)), MEAN(AP(A,Z) - AP(A,Z|CARESActChecks)), 
% MEAN(C(A,Z) - C(A,Z|CARESActChecks))
% 
% Tabulate value and policies along savings and shocks:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [1,4,5,6,3,2];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(v(A,Z) - v(A,Z|CARESActChecks))", V_VFI_wthtrumpchk_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(AP(A,Z) - AP(A,Z|CARESActChecks))", ap_VFI_wthtrumpchk_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
%% 
% 

% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(C(A,Z) - C(A,Z|CARESActChecks))", cons_VFI_wthtrumpchk_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(val(a,z) - val(a,z|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(val(a,z) - val(a,z|CARESActChecks))'};
ff_graph_grid((tb_az_v{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(a,z) - ap(a,z|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(a,z) - ap(a,z|CARESActChecks))'};
ff_graph_grid((tb_az_ap{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(C(a,z) - C(a,z|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(C(a,z) - C(a,z|CARESActChecks))'};
ff_graph_grid((tb_az_c{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% Analyze Kids and Marriage and Age
% Aggregating over education, savings, and shocks, what are the differential 
% effects of Marriage and Age.

% Generate some Data
mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
ar_row_grid = [...
    "k0M0", "K1M0", "K2M0", "K3M0", "K4M0", ...
    "k0M1", "K1M1", "K2M1", "K3M1", "K4M1"];
mp_support_graph('cl_st_xtitle') = {'Age'};
mp_support_graph('st_legend_loc') = 'best';
mp_support_graph('bl_graph_logy') = true; % do not log
mp_support_graph('st_rounding') = '6.2f'; % format shock legend
mp_support_graph('cl_scatter_shapes') = {...
    'o', 'd' ,'s', 'x', '*', ...
    'o', 'd', 's', 'x', '*'};
mp_support_graph('cl_colors') = {...
    'red', 'red', 'red', 'red', 'red'...
    'blue', 'blue', 'blue', 'blue', 'blue'};
%% 
% MEAN(V(KM,J) - V(KM,J | CARESActChecks)), MEAN(ap(KM,J) - ap(KM,J | CARESActChecks)), 
% MEAN(c(KM,J) - c(KM,J | CARESActChecks))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,4,1,6,5];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(V(KM,J) - V(KM,J | CARESActChecks))", V_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(ap(KM,J) - ap(KM,J | CARESActChecks))", ap_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% 

% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(c(KM,J) - c(KM,J | CARESActChecks))", cons_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(v(KM,J) - v(KM,J|CARESActChecks), a=age, z=kids+marry'};
mp_support_graph('cl_st_ytitle') = {'MEAN(v(KM,J) - v(KM,J|CARESActChecks)'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(KM,J) - ap(KM,J|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(KM,J) - ap(KM,J|CARESActChecks))'};
ff_graph_grid((tb_az_ap{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(c(KM,J) - c(KM,J|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(c(KM,J) - c(KM,J|CARESActChecks))'};
ff_graph_grid((tb_az_c{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% Analyze Education and Marriage and Age
% Aggregating over education, savings, and shocks, what are the differential 
% effects of Marriage and Age.

% Generate some Data
mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
ar_row_grid = ["E0M0", "E1M0", "E0M1", "E1M1"];
mp_support_graph('cl_st_xtitle') = {'Age'};
mp_support_graph('st_legend_loc') = 'best';
mp_support_graph('bl_graph_logy') = true; % do not log
mp_support_graph('st_rounding') = '6.2f'; % format shock legend
mp_support_graph('cl_scatter_shapes') = {'*', 'p', '*','p' };
mp_support_graph('cl_colors') = {'red', 'red', 'blue', 'blue'};
%% 
% MEAN(v(EKM,J) - v(EKM,J|CARESActChecks)), MEAN(ap(EM,J) - ap(EM,J|CARESActChecks)), 
% MEAN(c(EM,J) - c(EM,J|CARESActChecks))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,6,1,4,5];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(v(EM,J) - v(EM,J|CARESActChecks))", V_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(ap(EM,J) - ap(EM,J|CARESActChecks))", ap_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(c(EM,J) - c(EM,J|CARESActChecks))", cons_VFI_wthtrumpchk_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(v(EM,J) - v(EM,J|CARESActChecks)), a=age, z=kids+marry'};
mp_support_graph('cl_st_ytitle') = {'MEAN(v(EM,J) - v(EM,J|CARESActChecks))'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(EM,J) - ap(EM,J|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(EM,J) - ap(EM,J|CARESActChecks))'};
ff_graph_grid((tb_az_ap{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(c(EM,J) - c(EM,J|CARESActChecks)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(c(EM,J) - c(EM,J|CARESActChecks))'};
ff_graph_grid((tb_az_c{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);