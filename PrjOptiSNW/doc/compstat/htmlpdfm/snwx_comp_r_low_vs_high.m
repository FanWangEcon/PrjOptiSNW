%% Value and Consumption Low vs Higher Interest Rates Results Comparison
% This is the example vignette for function: <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m 
% *snw_vfi_main_bisec_vec*> from the <https://fanwangecon.github.io/PrjOptiSNW/ 
% *PrjOptiSNW Package*>*.* This function solves for the V(states) for individuals 
% at lower and higher savings interest rate. Note that welfare improves for all 
% when interest rate goes up for savings in a model where borrowing is not allowed. 
% However, an increase change in the interest rate generates both an income effect 
% (higher resources) and also changes relative price of consumption today vs tomorrow. 
% The change in income increases incentive to consume, the change in relative 
% price depresses incentives to consume today. The combined effect of rising interest 
% rate on savings on consumption/savings differs by the state-space, households 
% might overall consume more or less depending on their state-space.
%% Solve Model at 4 Percent Interest Rate
% Solve the benchmark model at 4 percent savings interest rate.

% mp_params = snw_mp_param('default_dense');
mp_params = snw_mp_param('default_docdense');
mp_params('beta') = 0.95;
fl_higher_r = 0.04;
fl_lower_r = 0.02;
mp_params('r') = fl_higher_r;
mp_controls = snw_mp_control('default_test');
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_timer') = true;
[V_ss_r04,~,cons_ss_r04,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
%% Solve Model at 2 Percent Interest Rate
% Solve the benchmark model at 2 percent savings interest rate.

mp_params('r') = fl_lower_r;
mp_controls = snw_mp_control('default_test');
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_timer') = true;
[V_ss_r02,~,cons_ss_r02,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
%% Generate Interest Rate Comparison Matrixes
% Take the difference between 4 percent and 2 percent savings interest rate 
% results. When interest rates are higher, greater incentive to save, but leads 
% to heterogeneous responses by income and other characteristics, note that this 
% changes both relative prices as well as total resource/budget, so there is both 
% income and price effects that differ in magnitudes depending on the individual's 
% statespace. Welfare does improve for all higher higher r. Welfare is converted 
% to units in fixed life-time consumption.

gamma = mp_params('gamma');
mn_V_gain_r = snw_hh_welfare(V_ss_r04, gamma) - snw_hh_welfare(V_ss_r02, gamma);
mn_C_gain_r = cons_ss_r04 - cons_ss_r02;
fl_r_gap = fl_higher_r - fl_lower_r;
%% Dense Param Results Define Frames
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
%% Analyze Difference in V and C with Higher and Lower Savings Interest Rate
% The difference between V and C with higher and lower $r$.

% Generate some Data
mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_support_graph('cl_st_xtitle') = {'Savings States, a'};
mp_support_graph('st_legend_loc') = 'eastoutside';
mp_support_graph('bl_graph_logy') = true; % do not log
mp_support_graph('it_legend_select') = 21; % how many shock legends to show
mp_support_graph('cl_colors') = 'jet';
%% 
% MEAN(MN_V_GAIN(A,Z))
% 
% Tabulate value and policies along savings and shocks:

% Set
ar_permute = [1,4,5,6,3,2];
% Value Function
st_title = ['MEAN(MN_V_Gain(A,Z)), r_gap=' num2str(fl_r_gap) ];
tb_az_v = ff_summ_nd_array(st_title, mn_V_gain_r, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
% Consumption
st_title = ['MEAN(MN_C_Gain(A,Z)), r_gap=' num2str(fl_r_gap) ];
tb_az_c = ff_summ_nd_array(st_title, mn_C_gain_r, true, ["mean"], 4, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values:

st_title = ['MEAN(MN\_V\_GAIN(A,Z)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_V\_GAIN(a,z))'};
ff_graph_grid((tb_az_v{1:end, 3:end})', ar_st_eta_HS_grid, agrid, mp_support_graph);
%% 
% Graph Mean Consumption:

st_title = ['MEAN(MN\_C\_GAIN(A,Z)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_C\_GAIN(a,z))'};
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
% MEAN(VAL(KM,J)), MEAN(AP(KM,J)), MEAN(C(KM,J))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,4,1,6,5];
% Value Function
st_title = ['MEAN(MN_V_Gain(KM,J)), r_gap=' num2str(fl_r_gap) ];
tb_az_v = ff_summ_nd_array(st_title, mn_V_gain_r, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Consumption Function
st_title = ['MEAN(MN_C_Gain(KM,J)), r_gap=' num2str(fl_r_gap) ];
tb_az_c = ff_summ_nd_array(st_title, mn_C_gain_r, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values:

st_title = ['MEAN(MN\_V\_GAIN(KM,J)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_V\_GAIN(KM,J))'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption:

st_title = ['MEAN(MN\_C\_GAIN(KM,J)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_C\_GAIN(KM,J))'};
ff_graph_grid((tb_az_c{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% Analyze Education and Marriage
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
% MEAN(VAL(EM,J)), MEAN(AP(EM,J)), MEAN(C(EM,J))
% 
% Tabulate value and policies:

% Set
% NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ar_permute = [2,3,6,1,4,5];
% Value Function
st_title = ['MEAN(MN_V_Gain(EM,J)), r_gap=' num2str(fl_r_gap) ];
tb_az_v = ff_summ_nd_array(st_title, mn_V_gain_r, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
% Consumption
st_title = ['MEAN(MN_C_Gain(EM,J)), r_gap=' num2str(fl_r_gap) ];
tb_az_c = ff_summ_nd_array(st_title, mn_C_gain_r, true, ["mean"], 3, 1, cl_mp_datasetdesc, ar_permute);
%% 
% Graph Mean Values:

st_title = ['MEAN(MN\_V\_GAIN(EM,J)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_V\_GAIN(EM,J))'};
ff_graph_grid((tb_az_v{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);
%% 
% Graph Mean Consumption:

st_title = ['MEAN(MN\_C\_GAIN(EM,J)), r\_gap=' num2str(fl_r_gap)  ''];
mp_support_graph('cl_st_graph_title') = {st_title};
mp_support_graph('cl_st_ytitle') = {'MEAN(MN\_C\_GAIN(EM,J))'};
ff_graph_grid((tb_az_c{1:end, 4:end}), ar_row_grid, age_grid, mp_support_graph);