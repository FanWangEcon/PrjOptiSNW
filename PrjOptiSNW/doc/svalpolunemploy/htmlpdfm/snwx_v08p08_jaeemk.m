%% Life Cycle Dynamic Programming under Great Recession Unemployment Shock
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_v08p08_jaeemk.m 
% *snw_v08p08_jaeemk*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*>*.* Solving the dynamic programming problem conditional on having an 
% one period unemployment shock that is expected with known unemployment probability. 
% Unemployment probability is a function of the realized state-space next year, 
% specifically, it is determined by age and education. Bush 2008 checks were received 
% by households in expectation of forth-coming unemployment shocks, ex-ante the 
% realization of shocks. During COVID, the shocks were received ex-post the realization 
% of shocks. In both cases, stimulus checks were determined based on ex-ante information.
% 
% Due to expected shock, households consume less and save more in 2008 than 
% under steady-state, as shown below. Value/welfare overall is lower in 2008 than 
% under steady-state.
%% Test SNW_V08P08_JAEEMK
% First, solve for value without unemployment issue (use the vectorized code 
% that was previously tested). This is the steady state results, but also the 
% results in 2009 without unemployment.

mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
mp_more_inputs('fl_ss_non_college') = 0.225;
mp_more_inputs('fl_ss_college') = 0.271;
mp_more_inputs('fl_scaleconvertor') = 54831;
% st_param_group = 'default_small';
% st_param_group = 'default_dense';
st_param_group = 'default_docdense';
mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
mp_controls = snw_mp_control('default_test');
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
[V_VFI_ss, ap_VFI_ss, cons_VFI_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
V_emp_2009 = V_VFI_ss;
%% 
% Second, solve for the unemployment value, use the exact-bisec result code, 
% call the snw_vfi_main_bisec_vec.m function with a third input of existing value. 
% xi is the share of income lost during covid year given surprise covid shock, 
% b is the share of income loss that is covered by unemployment insurance. If 
% xi=0.5 and b=0 means will lose 50 percent of income given 2009 great recession 
% shocks, and the loss will not be covered at all by unemployment insurance.

mp_params('xi') = 0.532;
mp_params('b') = 0.37992;
mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
[V_unemp_2009] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_VFI_ss);
%% 
% Third, solve for 2008 policy and value funtion given employed and unemployed 
% value function in 2009, 

[V_2008, ap_2008, cons_2008, ev_empshk_2009] = ...
    snw_v08p08_jaeemk(mp_params, mp_controls, V_emp_2009, V_unemp_2009);
%% 
% Difference Between Value and Choices In steady state and in 2008, given expected 
% unemployment (one-period) shock due to the great recession,  <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_v08p08_jaeemk.m 
% *snw_v08p08_jaeemk*>.

V_VFI_unemp_drop = V_VFI_ss - V_2008;
ap_VFI_unemp_drop = ap_VFI_ss - ap_2008;
cons_VFI_unemp_drop = cons_VFI_ss - cons_2008;
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
% MEAN(VAL(A,Z) - VAL(A,Z, 08wthEV09unemshk)), MEAN(AP(A,Z) - AP(A,Z, 08wthEV09unemshk)), 
% MEAN(C(A,Z) - C(A,Z, 08wthEV09unemshk))
% 
% Tabulate value and policies along savings and shocks:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [1,4,5,6,3,2];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(v(A,Z) - v(A,Z, 08wthEV09unemshk))", V_VFI_unemp_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(AP(A,Z) - AP(A,Z, 08wthEV09unemshk))", ap_VFI_unemp_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
%% 
% 

% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(C(A,Z) - C(A,Z, 08wthEV09unemshk))", cons_VFI_unemp_drop, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(val(a,z) - val(a,z, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(val(a,z) - val(a,z, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_v{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(a,z) - ap(a,z, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(a,z) - ap(a,z, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_ap{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(C(a,z) - C(a,z, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(C(a,z) - C(a,z, 08wthEV09unemshk))'};
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
% MEAN(V(KM,J) - V(KM,J | unemp)), MEAN(ap(KM,J) - ap(KM,J | unemp)), MEAN(c(KM,J) 
% - c(KM,J | unemp))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,4,1,6,5];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(V(KM,J) - V(KM,J | unemp))", V_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(ap(KM,J) - ap(KM,J | unemp))", ap_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% 

% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(c(KM,J) - c(KM,J | unemp))", cons_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(v(KM,J) - v(KM,J, 08wthEV09unemshk), a=age, z=kids+marry'};
mp_support_graph('cl_st_ytitle') = {'MEAN(v(KM,J) - v(KM,J, 08wthEV09unemshk)'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(KM,J) - ap(KM,J, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(KM,J) - ap(KM,J, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_ap{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(c(KM,J) - c(KM,J, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(c(KM,J) - c(KM,J, 08wthEV09unemshk))'};
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
% MEAN(v(EKM,J) - v(EKM,J, 08wthEV09unemshk)), MEAN(ap(EM,J, steady) - ap(EM,J, 
% 08wthEV09unemshk)), MEAN(c(EM,J, steady) - c(EM,J, 08wthEV09unemshk))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,6,1,4,5];
% Value Function
tb_az_v = ff_summ_nd_array("MEAN(v(EM,J, steady) - v(EM,J, 08wthEV09unemshk))", V_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Aprime Choice
tb_az_ap = ff_summ_nd_array("MEAN(ap(EM,J, steady) - ap(EM,J, 08wthEV09unemshk))", ap_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Consumption Choices
tb_az_c = ff_summ_nd_array("MEAN(c(EM,J, steady) - c(EM,J, 08wthEV09unemshk))", cons_VFI_unemp_drop, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(v(EM,J, steady) - v(EM,J, 08wthEV09unemshk)), a=age, z=kids+marry'};
mp_support_graph('cl_st_ytitle') = {'MEAN(v(EM,J, steady) - v(EM,J, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Savings Choices Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(ap(EM,J, steady) - ap(EM,J, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(ap(EM,J, steady) - ap(EM,J, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_ap{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption Change:

mp_support_graph('cl_st_graph_title') = {'MEAN(c(EM,J, steady) - c(EM,J, 08wthEV09unemshk)), a=x, z=color'};
mp_support_graph('cl_st_ytitle') = {'MEAN(c(EM,J, steady) - c(EM,J, 08wthEV09unemshk))'};
ff_graph_grid((tb_az_c{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);