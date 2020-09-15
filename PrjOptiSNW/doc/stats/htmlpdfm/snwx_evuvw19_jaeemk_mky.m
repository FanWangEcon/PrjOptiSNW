%% 2019 Full States MPC and Distributional Statistics by Marital, Kids, and Income Groups.
% In the file here, we consider marital, kids and income groups, and summarize 
% various statistics for each bin. 
%% Test SNW_EVUVW19_JAEEMK Defaults Dense
% VFI and Distribution
% 
% Call the function with defaults.

clear all;
st_solu_type = 'bisec_vec';
bl_save_csv = false;

% Solve the VFI Problem and get Value Function
% mp_params = snw_mp_param('default_dense');
% mp_params = snw_mp_param('default_docdense');
mp_params = snw_mp_param('default_moredense_a65zh133zs5_e2m2');
mp_controls = snw_mp_control('default_test');

% set Unemployment Related Variables
xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values

mp_params('xi') = xi;
mp_params('b') = b;
mp_params('TR') = TR;

% Solve for Unemployment Values
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = true;
mp_controls('bl_print_ds_verbose') = true;    
mp_controls('bl_print_precompute') = false;
mp_controls('bl_print_precompute_verbose') = false;
mp_controls('bl_print_a4chk') = false;
mp_controls('bl_print_a4chk_verbose') = false;
mp_controls('bl_print_evuvw20_jaeemk') = false;
mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
mp_controls('bl_print_evuvw19_jaeemk') = false;
mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;

% Solve the Model to get V working and unemployed
[V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
inc_VFI = mp_valpol_more_ss('inc_VFI');
spouse_inc_VFI = mp_valpol_more_ss('spouse_inc_VFI');
total_inc_VFI = inc_VFI + spouse_inc_VFI;
% tax during covid year
mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');
% Solve unemployment
[V_unemp,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
[Phi_true, Phi_adj, A_agg, Y_inc_agg, ~, mp_dsvfi_results] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
% Get Matrixes
cl_st_precompute_list = {'a', 'ar_z_ctr_amz', ...
    'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid',...
    'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss'};
mp_controls('bl_print_precompute_verbose') = false;
[mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);
%% Solve for 2019 Evuvw With 0 and 1 Checks

% Call Function
welf_checks = 0;
[ev19_jaeemk_check0, ec19_jaeemk_check0, ev20_jaeemk_check0, ec20_jaeemk_check0] = snw_evuvw19_jaeemk_foc(...
    welf_checks, st_solu_type, mp_params, mp_controls, ...
    V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, mp_precompute_res);
% Call Function
welf_checks = 1;
[ev19_jaeemk_check2, ec19_jaeemk_check2, ev20_jaeemk_check2, ec20_jaeemk_check2] = snw_evuvw19_jaeemk_foc(...
    welf_checks, st_solu_type, mp_params, mp_controls, ...
    V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, mp_precompute_res);
%% 
% Differences between Checks in Expected Value and Expected Consumption

mn_V_U_gain_check = ev19_jaeemk_check2 - ev19_jaeemk_check0;
mn_MPC_C_gain_share_check = (ec19_jaeemk_check2 - ec19_jaeemk_check0)./(welf_checks*mp_params('TR'));
%% Additional Variables
% Create additional Staet-Spac Arrays

% (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% Children Array
ar_kids = (1:mp_params('n_kidsgrid')) - 1;
mn_kids = zeros(1,1,1,1,1,length(ar_kids));
mn_kids(1,1,1,1,1,:) = ar_kids;
kids_ss = repmat(mn_kids, [mp_params('n_jgrid'), mp_params('n_agrid'), mp_params('n_etagrid'), ...
    mp_params('n_educgrid'), mp_params('n_marriedgrid'), 1]);
% Marital Status Arrays
ar_marital = (1:mp_params('n_marriedgrid')) - 1;
mn_marital = zeros(1,1,1,1,length(ar_marital),1);
mn_marital(1,1,1,1,:,1) = ar_marital;
marital_ss = repmat(mn_marital, [mp_params('n_jgrid'), mp_params('n_agrid'), mp_params('n_etagrid'), ...
    mp_params('n_educgrid'), 1, mp_params('n_kidsgrid')]);
% Educational Status Arrays
ar_educ = (1:mp_params('n_educgrid')) - 1;
mn_educ = zeros(1,1,1,length(ar_educ),1,1);
mn_educ(1,1,1,:,1,1) = ar_educ;
educ_ss = repmat(mn_educ, [mp_params('n_jgrid'), mp_params('n_agrid'), mp_params('n_etagrid'), ...
    1, mp_params('n_marriedgrid'), mp_params('n_kidsgrid')]);
% Age Array
ar_age = (1:mp_params('n_jgrid')) + 18;
mn_age = zeros(length(ar_age),1,1,1,1,1);
mn_age(:,1,1,1,1,1) = ar_age;
age_ss = repmat(mn_age, [1, mp_params('n_agrid'), mp_params('n_etagrid'), ...
    mp_params('n_educgrid'), mp_params('n_marriedgrid'), mp_params('n_kidsgrid')]);
%% Adjust to Probability Mass Function

Phi_true_1 = Phi_true./sum(Phi_true,'all');
%% Age Bounds

% 1 = 18
min_age = 1
% retirement, 46+18=64, the year prior to retirement year.
max_age = 46;
%% Scale Statistics to Thousands of Dollars

a_ss = mp_dsvfi_results('a_ss')*58.056;
ap_ss = mp_dsvfi_results('ap_ss')*58.056;
c_ss = mp_dsvfi_results('cons_ss')*58.056;
n_ss = mp_dsvfi_results('n_ss');
% household head + spousal (realized) income
y_all = mp_dsvfi_results('y_all_ss')*58.056;
y_head_inc = mp_dsvfi_results('y_head_inc_ss')*58.056;
y_spouse_inc = mp_dsvfi_results('y_spouse_inc_ss')*58.056;

yshr_wage = mp_dsvfi_results('yshr_wage_ss');
yshr_SS = mp_dsvfi_results('yshr_SS_ss');
yshr_nttxss = mp_dsvfi_results('yshr_nttxss_ss');
%% Distributional Statistics Overall All Ages

% construct input data
marital_grp = marital_ss(min_age:82, :, :, : ,: ,:);
y_all_grp = y_all(min_age:82, :, :, : ,: ,:);
age_ss_grp = age_ss(min_age:82, :, :, : ,: ,:);
educ_ss_grp = educ_ss(min_age:82, :, :, : ,: ,:);
a_ss_grp = a_ss(min_age:82, :, :, : ,: ,:);
ap_ss_grp = ap_ss(min_age:82, :, :, : ,: ,:);
mn_MPC_C_gain_share_check_grp = mn_MPC_C_gain_share_check(min_age:82, :, :, :, :, :);        
Phi_true_grp = Phi_true_1(min_age:82, :, :, : ,: ,:);
c_ss_grp = c_ss(min_age:82, :, :, : ,: ,:);
y_head_inc_grp = y_head_inc(min_age:82, :, :, : ,: ,:);
y_spouse_inc_grp = y_spouse_inc(min_age:82, :, :, : ,: ,:);
yshr_nttxss_grp = yshr_nttxss(min_age:82, :, :, : ,: ,:);    

mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('married') = {marital_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_all') = {y_all_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('age_ss') = {age_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('educ_ss') = {educ_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('a_ss') = {a_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('MPC') = {mn_MPC_C_gain_share_check_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('Mass') = {Phi_true_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('c_ss') = {c_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_head_inc') = {y_head_inc_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_spouse') = {y_spouse_inc_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('yshr_nttxss') = {yshr_nttxss_grp(:), zeros(1)};

mp_cl_ar_xyz_of_s('ar_st_y_name') = ["married", "y_all", "age_ss", "educ_ss", "a_ss", "ap_ss", "MPC", "Mass", "c_ss", "y_head_inc", "y_spouse", "yshr_nttxss"];

% controls
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('ar_fl_percentiles') = [0.01 10 25 50 75 90 99.99];
mp_support('bl_display_final') = true;
mp_support('bl_display_detail') = false;    
mp_support('bl_display_drvm2outcomes') = false;  
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;

% Call Function     
mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);
tb_dist_stats_all = mp_cl_mt_xyz_of_s('tb_outcomes');
%% Distributional Statistics Overall 18 to 64
% Statistics overall distributionally for 18 to 64 year olds.

% construct input data
marital_grp = marital_ss(min_age:max_age, :, :, : ,: ,:);
y_all_grp = y_all(min_age:max_age, :, :, : ,: ,:);
age_ss_grp = age_ss(min_age:max_age, :, :, : ,: ,:);
educ_ss_grp = educ_ss(min_age:max_age, :, :, : ,: ,:);
a_ss_grp = a_ss(min_age:max_age, :, :, : ,: ,:);
ap_ss_grp = ap_ss(min_age:max_age, :, :, : ,: ,:);
mn_MPC_C_gain_share_check_grp = mn_MPC_C_gain_share_check(min_age:max_age, :, :, :, :, :);        
Phi_true_grp = Phi_true_1(min_age:max_age, :, :, : ,: ,:);
c_ss_grp = c_ss(min_age:max_age, :, :, : ,: ,:);
y_head_inc_grp = y_head_inc(min_age:max_age, :, :, : ,: ,:);
y_spouse_inc_grp = y_spouse_inc(min_age:max_age, :, :, : ,: ,:);
yshr_nttxss_grp = yshr_nttxss(min_age:max_age, :, :, : ,: ,:);    

mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('married') = {marital_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_all') = {y_all_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('age_ss') = {age_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('educ_ss') = {educ_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('a_ss') = {a_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('MPC') = {mn_MPC_C_gain_share_check_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('Mass') = {Phi_true_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('c_ss') = {c_ss_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_head_inc') = {y_head_inc_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('y_spouse') = {y_spouse_inc_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('yshr_nttxss') = {yshr_nttxss_grp(:), zeros(1)};

mp_cl_ar_xyz_of_s('ar_st_y_name') = ["married", "y_all", "age_ss", "educ_ss", "a_ss", "ap_ss", "MPC", "Mass", "c_ss", "y_head_inc", "y_spouse", "yshr_nttxss"];

% controls
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('ar_fl_percentiles') = [0.01 10 25 50 75 90 99.99];
mp_support('bl_display_final') = true;
mp_support('bl_display_detail') = false;    
mp_support('bl_display_drvm2outcomes') = false;  
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;

% Call Function     
mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);
tb_dist_stats_all_18to64 = mp_cl_mt_xyz_of_s('tb_outcomes');
%% Distributional Statistics By Kids Count
% Various statistics, including MPC (of the first check) by Children Count

it_row_ctr = 0;
for it_ctr=1:mp_params('n_kidsgrid')
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    display(['kids =' num2str(ar_kids(it_ctr))]);
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    
    % construct input data
    marital_grp = marital_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    y_all_grp = y_all(min_age:max_age, :, :, : ,: ,it_ctr);
    age_ss_grp = age_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    educ_ss_grp = educ_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    a_ss_grp = a_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    ap_ss_grp = ap_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    mn_MPC_C_gain_share_check_grp = mn_MPC_C_gain_share_check(min_age:max_age, :, :, :, :, it_ctr);        
    Phi_true_grp = Phi_true_1(min_age:max_age, :, :, : ,: ,it_ctr);
    c_ss_grp = c_ss(min_age:max_age, :, :, : ,: ,it_ctr);
    y_head_inc_grp = y_head_inc(min_age:max_age, :, :, : ,: ,it_ctr);
    y_spouse_inc_grp = y_spouse_inc(min_age:max_age, :, :, : ,: ,it_ctr);
    yshr_nttxss_grp = yshr_nttxss(min_age:max_age, :, :, : ,: ,it_ctr);    
    
    mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
    mp_cl_ar_xyz_of_s('married') = {marital_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('y_all') = {y_all_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('age_ss') = {age_ss_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('educ_ss') = {educ_ss_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('a_ss') = {a_ss_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('MPC') = {mn_MPC_C_gain_share_check_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('Mass') = {Phi_true_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('c_ss') = {c_ss_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('y_head_inc') = {y_head_inc_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('y_spouse') = {y_spouse_inc_grp(:), zeros(1)};
    mp_cl_ar_xyz_of_s('yshr_nttxss') = {yshr_nttxss_grp(:), zeros(1)};
    
    mp_cl_ar_xyz_of_s('ar_st_y_name') = ["married", "y_all", "age_ss", "educ_ss", "a_ss", "ap_ss", "MPC", "Mass", "c_ss", "y_head_inc", "y_spouse", "yshr_nttxss"];
    
    % controls
    mp_support = containers.Map('KeyType','char', 'ValueType','any');
    mp_support('ar_fl_percentiles') = [0.01 10 25 50 75 90 99.99];
    mp_support('bl_display_final') = true;
    mp_support('bl_display_detail') = false;    
    mp_support('bl_display_drvm2outcomes') = false;  
    mp_support('bl_display_drvstats') = false;
    mp_support('bl_display_drvm2covcor') = false;

    % Call Function     
    mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);

    it_kids = ar_kids(it_ctr);
    
    tb_dist_stats = mp_cl_mt_xyz_of_s('tb_outcomes');
    
    fl_married_mean = tb_dist_stats{"married", "mean"};
    
    fl_age_mean = tb_dist_stats{"age_ss", "mean"};
    fl_age_p50 = tb_dist_stats{"age_ss", "p50"};
   	
    fl_educ_mean = tb_dist_stats{"educ_ss", "mean"};
    
    fl_a_mean = tb_dist_stats{"a_ss", "mean"};
    fl_a_p50 = tb_dist_stats{"a_ss", "p50"};
    
    fl_ap_mean = tb_dist_stats{"ap_ss", "mean"};
    fl_ap_p50 = tb_dist_stats{"ap_ss", "p50"};
    
    fl_y_all_mean = tb_dist_stats{"y_all", "mean"};
    fl_y_all_p50 = tb_dist_stats{"y_all", "p50"};
    
    fl_mpc_mean = tb_dist_stats{"MPC", "mean"};
    fl_mpc_p50 = tb_dist_stats{"MPC", "p50"};
    
    fl_mass = tb_dist_stats{"Mass", "unweighted_sum"};
    
    fl_c_ss_mean = tb_dist_stats{"c_ss", "mean"};
    fl_c_ss_p50 = tb_dist_stats{"c_ss", "p50"};
    
    fl_y_head_inc_mean = tb_dist_stats{"y_head_inc", "mean"};
    fl_y_spouse_mean = tb_dist_stats{"y_spouse", "mean"};
                            
    ar_store_stats = [it_kids, fl_married_mean, ...
       	fl_age_mean, fl_age_p50, fl_educ_mean, ...
       	fl_a_mean, fl_a_p50, fl_ap_mean, fl_ap_p50, ...
       	fl_y_all_mean, fl_y_all_p50, ...
       	fl_mpc_mean, fl_mpc_p50, ...
        fl_mass, ...
       	fl_c_ss_mean, fl_c_ss_p50, ...
       	fl_y_head_inc_mean, fl_y_spouse_mean];
    
    it_row_ctr = it_row_ctr + 1;
    
    if (it_row_ctr>1)
        mt_store_stats_by_k = [mt_store_stats_by_k;ar_store_stats];
    else
        mt_store_stats_by_k = [ar_store_stats];
    end    
end
%% Distributional Statistics By Marital Status and Kids Count
% Various statistics, including MPC (of the first check) by Marital Status and 
% Kids COunt

it_row_ctr = 0;
for it_marry_ctr=1:mp_params('n_marriedgrid')
    
    display(['']);
    display(['']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    display(['Marital =' num2str(ar_marital(it_marry_ctr))]);
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    
    for it_kids_ctr=1:mp_params('n_kidsgrid')
        display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
        display(['Marital =' num2str(ar_marital(it_marry_ctr)) ' and kids =' num2str(ar_kids(it_kids_ctr))]);
        display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
        
        % construct input data
        y_all_grp = y_all(min_age:max_age, :, :, : ,it_marry_ctr ,it_ctr);
        age_ss_grp = age_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        educ_ss_grp = educ_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        a_ss_grp = a_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        ap_ss_grp = ap_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        mn_MPC_C_gain_share_check_grp = mn_MPC_C_gain_share_check(min_age:max_age, :, :, :, it_marry_ctr, it_kids_ctr);
        Phi_true_grp = Phi_true_1(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        c_ss_grp = c_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        y_head_inc_grp = y_head_inc(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        y_spouse_inc_grp = y_spouse_inc(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        yshr_nttxss_grp = yshr_nttxss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        
        mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
        mp_cl_ar_xyz_of_s('y_all') = {y_all_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('age_ss') = {age_ss_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('educ_ss') = {educ_ss_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('a_ss') = {a_ss_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('MPC') = {mn_MPC_C_gain_share_check_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('Mass') = {Phi_true_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('c_ss') = {c_ss_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_head_inc') = {y_head_inc_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_spouse') = {y_spouse_inc_grp(:), zeros(1)};
        mp_cl_ar_xyz_of_s('yshr_nttxss') = {yshr_nttxss_grp(:), zeros(1)};
        
        mp_cl_ar_xyz_of_s('ar_st_y_name') = ["y_all", "age_ss", "educ_ss", "a_ss", "ap_ss", "MPC", "Mass", "c_ss", "y_head_inc", "y_spouse", "yshr_nttxss"];
        
        % controls
        mp_support = containers.Map('KeyType','char', 'ValueType','any');
        mp_support('ar_fl_percentiles') = [0.01 10 25 50 75 90 99.99];
        mp_support('bl_display_final') = true;
        mp_support('bl_display_detail') = false;
        mp_support('bl_display_drvm2outcomes') = false;
        mp_support('bl_display_drvstats') = false;
        mp_support('bl_display_drvm2covcor') = false;
        
        % Call Function
        mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);
            
        it_marital = ar_marital(it_marry_ctr);
        it_kids = ar_kids(it_kids_ctr);
        
        tb_dist_stats = mp_cl_mt_xyz_of_s('tb_outcomes');
        fl_age_mean = tb_dist_stats{"age_ss", "mean"};
        fl_age_p50 = tb_dist_stats{"age_ss", "p50"};
       	
        fl_educ_mean = tb_dist_stats{"educ_ss", "mean"};
        
        fl_a_mean = tb_dist_stats{"a_ss", "mean"};
        fl_a_p50 = tb_dist_stats{"a_ss", "p50"};
        
        fl_ap_mean = tb_dist_stats{"ap_ss", "mean"};
        fl_ap_p50 = tb_dist_stats{"ap_ss", "p50"};
        
        fl_y_all_mean = tb_dist_stats{"y_all", "mean"};
        fl_y_all_p50 = tb_dist_stats{"y_all", "p50"};
        
        fl_mpc_mean = tb_dist_stats{"MPC", "mean"};
        fl_mpc_p50 = tb_dist_stats{"MPC", "p50"};
        
        fl_mass = tb_dist_stats{"Mass", "unweighted_sum"};
        
        fl_c_ss_mean = tb_dist_stats{"c_ss", "mean"};
        fl_c_ss_p50 = tb_dist_stats{"c_ss", "p50"};
        
        fl_y_head_inc_mean = tb_dist_stats{"y_head_inc", "mean"};
        fl_y_spouse_mean = tb_dist_stats{"y_spouse", "mean"};
                                
        ar_store_stats = [it_marital, it_kids, ...
           	fl_age_mean, fl_age_p50, fl_educ_mean, ...
           	fl_a_mean, fl_a_p50, fl_ap_mean, fl_ap_p50, ...
           	fl_y_all_mean, fl_y_all_p50, ...
           	fl_mpc_mean, fl_mpc_p50, ...
            fl_mass, ...
           	fl_c_ss_mean, fl_c_ss_p50, ...
           	fl_y_head_inc_mean, fl_y_spouse_mean];
        
        it_row_ctr = it_row_ctr + 1;
        
        if (it_row_ctr>1)
            mt_store_stats_by_mk = [mt_store_stats_by_mk;ar_store_stats];
        else
            mt_store_stats_by_mk = [ar_store_stats];
        end
    end
end
%% Distributional Statistics By Marital Status, Kids Count and Income Bins
% Various statistics, including MPC (of the first check) by Marital Status and 
% Kids COunt and income bins

it_row_ctr = 0;
for it_marry_ctr=1:mp_params('n_marriedgrid')
    
    display(['']);
    display(['']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    display(['Marital =' num2str(ar_marital(it_marry_ctr))]);
    display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    display(['-------------------------------']);
    display(['-------------------------------']);
    
    for it_kids_ctr=1:mp_params('n_kidsgrid')
        display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
        display(['Marital =' num2str(ar_marital(it_marry_ctr)) ' and kids =' num2str(ar_kids(it_kids_ctr))]);
        display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
        
        % construct input data
        y_all_grp = y_all(min_age:max_age, :, :, : ,it_marry_ctr ,it_ctr);
        age_ss_grp = age_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        educ_ss_grp = educ_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        a_ss_grp = a_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        ap_ss_grp = ap_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        mn_MPC_C_gain_share_check_grp = mn_MPC_C_gain_share_check(min_age:max_age, :, :, :, it_marry_ctr, it_kids_ctr);
        Phi_true_grp = Phi_true_1(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        c_ss_grp = c_ss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        y_head_inc_grp = y_head_inc(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        y_spouse_inc_grp = y_spouse_inc(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        yshr_nttxss_grp = yshr_nttxss(min_age:max_age, :, :, : ,it_marry_ctr, it_kids_ctr);
        
        % Income Bins
        ar_y_all = y_all_grp(:);
        ar_age_ss = age_ss_grp(:);
        ar_educ_ss = educ_ss_grp(:);
        ar_a_ss = a_ss_grp(:);
        ar_ap_ss = ap_ss_grp(:);
        ar_mn_MPC_C_gain_share_check = mn_MPC_C_gain_share_check_grp(:);
        ar_Phi_true = Phi_true_grp(:);
        ar_c_ss = c_ss_grp(:);
        ar_y_head_inc = y_head_inc_grp(:);
        ar_y_spouse_inc = y_spouse_inc_grp(:);
        ar_yshr_nttxss = yshr_nttxss_grp(:);
        
        % income bins loop
        for it_y_all_ctr=1:6
            
            % Current y group index
            % y is in thousands of dollars
            y_all_start = (it_y_all_ctr-1)*20;
            if (it_y_all_ctr == 6)
                y_all_end = max(ar_y_all);
            else
                y_all_end = it_y_all_ctr*20;                
            end
            
            display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
            display(['Marital =' num2str(ar_marital(it_marry_ctr)) ', kids =' num2str(ar_kids(it_kids_ctr)) ', ybin =' num2str(y_all_start) ' to ' num2str(y_all_end)]);
            display(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
            
            ar_y_idx = (ar_y_all >= y_all_start & ar_y_all <y_all_end);
            
            ar_mky_y_all = ar_y_all(ar_y_idx);
            ar_mky_age_ss = ar_age_ss(ar_y_idx);
            ar_mky_educ_ss = ar_educ_ss(ar_y_idx);
            ar_mky_a_ss = ar_a_ss(ar_y_idx);
            ar_mky_ap_ss = ar_ap_ss(ar_y_idx);
            ar_mky_mn_MPC_C_gain_share_check = ar_mn_MPC_C_gain_share_check(ar_y_idx);
            ar_mky_Phi_true = ar_Phi_true(ar_y_idx);
            ar_mky_c_ss = ar_c_ss(ar_y_idx);
            ar_mky_y_head_inc = ar_y_head_inc(ar_y_idx);
            ar_mky_y_spouse_inc = ar_y_spouse_inc(ar_y_idx);
            ar_mky_yshr_nttxss = ar_yshr_nttxss(ar_y_idx);
            
            mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
            mp_cl_ar_xyz_of_s('y_all') = {ar_mky_y_all(:), zeros(1)};
            mp_cl_ar_xyz_of_s('age_ss') = {ar_mky_age_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('educ_ss') = {ar_mky_educ_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('a_ss') = {ar_mky_a_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('ap_ss') = {ar_mky_ap_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('MPC') = {ar_mky_mn_MPC_C_gain_share_check(:), zeros(1)};
            mp_cl_ar_xyz_of_s('Mass') = {ar_mky_Phi_true(:), zeros(1)};
            mp_cl_ar_xyz_of_s('c_ss') = {ar_mky_c_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('y_head_inc') = {ar_mky_y_head_inc(:), zeros(1)};
            mp_cl_ar_xyz_of_s('y_spouse') = {ar_mky_y_spouse_inc(:), zeros(1)};
            mp_cl_ar_xyz_of_s('yshr_nttxss') = {ar_mky_yshr_nttxss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('ar_st_y_name') = ["y_all", "age_ss", "educ_ss", "a_ss", "ap_ss", "MPC", "Mass", "c_ss", "y_head_inc", "y_spouse", "yshr_nttxss"];
            
            % controls
            mp_support = containers.Map('KeyType','char', 'ValueType','any');
            mp_support('ar_fl_percentiles') = [0.01 10 25 50 75 90 99.99];
            mp_support('bl_display_final') = true;
            mp_support('bl_display_detail') = false;
            mp_support('bl_display_drvm2outcomes') = false;
            mp_support('bl_display_drvstats') = false;
            mp_support('bl_display_drvm2covcor') = false;
            
            % Call Function
            mp_cl_mt_xyz_of_s = ff_simu_stats(ar_mky_Phi_true(:)/sum(ar_mky_Phi_true,'all'), mp_cl_ar_xyz_of_s, mp_support);            
            
            it_marital = ar_marital(it_marry_ctr);
            it_kids = ar_kids(it_kids_ctr);
            fl_y_all_start = y_all_start;
            fl_y_all_end = y_all_end;
            
            tb_dist_stats = mp_cl_mt_xyz_of_s('tb_outcomes');
            fl_age_mean = tb_dist_stats{"age_ss", "mean"};
            fl_age_p50 = tb_dist_stats{"age_ss", "p50"};
           	
            fl_educ_mean = tb_dist_stats{"educ_ss", "mean"};
            
            fl_a_mean = tb_dist_stats{"a_ss", "mean"};
            fl_a_p50 = tb_dist_stats{"a_ss", "p50"};
            
            fl_ap_mean = tb_dist_stats{"ap_ss", "mean"};
            fl_ap_p50 = tb_dist_stats{"ap_ss", "p50"};
            
            fl_y_all_mean = tb_dist_stats{"y_all", "mean"};
            fl_y_all_p50 = tb_dist_stats{"y_all", "p50"};
            
            fl_mpc_mean = tb_dist_stats{"MPC", "mean"};
            fl_mpc_p50 = tb_dist_stats{"MPC", "p50"};
            
            fl_mass = tb_dist_stats{"Mass", "unweighted_sum"};
            
            fl_c_ss_mean = tb_dist_stats{"c_ss", "mean"};
            fl_c_ss_p50 = tb_dist_stats{"c_ss", "p50"};
            
            fl_y_head_inc_mean = tb_dist_stats{"y_head_inc", "mean"};
            fl_y_spouse_mean = tb_dist_stats{"y_spouse", "mean"};
                                    
            ar_store_stats = [it_marital, it_kids, fl_y_all_start, fl_y_all_end, ...
               	fl_age_mean, fl_age_p50, fl_educ_mean, ...
               	fl_a_mean, fl_a_p50, fl_ap_mean, fl_ap_p50, ...
               	fl_y_all_mean, fl_y_all_p50, ...
               	fl_mpc_mean, fl_mpc_p50, ...
                fl_mass, ...
               	fl_c_ss_mean, fl_c_ss_p50, ...
               	fl_y_head_inc_mean, fl_y_spouse_mean];
            
            it_row_ctr = it_row_ctr + 1;
            
            if (it_row_ctr>1)
                mt_store_stats_by_mky = [mt_store_stats_by_mky;ar_store_stats];
            else
                mt_store_stats_by_mky = [ar_store_stats];
            end
            
        end       
    end
end
%% Store Aggregate To File
% Store Several Files:
%% 
% # Overall Aggregate Statistics All Distribution
% # Aggregate Statistics Only for 18 to 64 year olds
% # Group Statistics by Kids
% # Group Statistics by Marital + Kids
% # Group Statistics by Marital + Kids + Income Bins

if (bl_save_csv)
    % All Stats All Ages
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_dist_stats_all, [spt_simu_results_csv 'stats_all_allages.csv'], 'WriteRowNames', true);
    % All Stats 18 to 64 Year old
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_dist_stats_all_18to64, [spt_simu_results_csv 'stats_all_18t64.csv'], 'WriteRowNames', true);
    % Group by K: Kids only
    tb_store_stats_by_k = array2table(mt_store_stats_by_k, 'VariableNames', ...
        {'kids', 'married_mean' ...
        'age_mean', 'age_p50', 'educ_mean', ...
        'a_mean', 'a_p50', 'ap_mean', 'ap_p50', ...
        'y_all_mean', 'y_all_p50', ...
        'mpc_mean', 'mpc_p50', ...
        'mass',...
        'c_ss_mean', 'c_ss_p50', ...
        'y_head_inc_mean', 'y_spouse_mean'});
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_store_stats_by_k, [spt_simu_results_csv 'stats_by_kids.csv']);
    % Group by MK: marry + kids only
    tb_store_stats_by_mk = array2table(mt_store_stats_by_mk, 'VariableNames', ...
        {'marital', 'kids', ...
        'age_mean', 'age_p50', 'educ_mean', ...
        'a_mean', 'a_p50', 'ap_mean', 'ap_p50', ...
        'y_all_mean', 'y_all_p50', ...
        'mpc_mean', 'mpc_p50', ...
        'mass',...
        'c_ss_mean', 'c_ss_p50', ...
        'y_head_inc_mean', 'y_spouse_mean'});
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_store_stats_by_mk, [spt_simu_results_csv 'stats_by_marital_kids.csv']);
    % Group by MKY
    tb_store_stats_by_mky = array2table(mt_store_stats_by_mky, 'VariableNames', ...
        {'marital', 'kids', 'y_all_start', 'y_all_end', ...
        'age_mean', 'age_p50', 'educ_mean', ...
        'a_mean', 'a_p50', 'ap_mean', 'ap_p50', ...
        'y_all_mean', 'y_all_p50', ...
        'mpc_mean', 'mpc_p50', ...
        'mass',...
        'c_ss_mean', 'c_ss_p50', ...
        'y_head_inc_mean', 'y_spouse_mean'});
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_store_stats_by_mky, [spt_simu_results_csv 'stats_by_marital_kids_20kincbins.csv']);
end
%% Store Key Stats to Compare to Key US Distributional Statistics
% Earning, income and Wealth.
% 
% Income = interest earnings + Social Security + labor income + spousal income. 
% This is equal to y_all.
% 
% Earnings = labor income + spousal income. 

% Income Variable
if (min(abs(total_inc_VFI*58.056 - y_all), [], 'all')>0)
    error('someothing is wrong, total_inc_VFI should be equal to y_all');
end
income = y_all;
% Earning variable
% earn*fl_earn_ratio generated earn_VFI
earning = (mp_valpol_more_ss('earn_VFI') + spouse_inc_VFI)*58.056;
% Wealth Varaible 
wealth = a_ss;
%% 
% Generate Key Statistics for these three variables only, distributional Statistics 
% Overall All Ages:

% construct input data
income_grp = income(min_age:82, :, :, : ,: ,:);
earning_grp = earning(min_age:82, :, :, : ,: ,:);
wealth_grp = wealth(min_age:82, :, :, : ,: ,:);
Phi_true_grp = Phi_true_1(min_age:82, :, :, : ,: ,:);

mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('earning') = {earning_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('income') = {income_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('wealth') = {wealth_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('earninglog') = {log(earning_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('incomelog') = {log(income_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('wealthlog') = {log(wealth_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('ar_st_y_name') = ["earning", "income", "wealth", "earninglog", "incomelog", "wealthlog"];

% controls
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('ar_fl_percentiles') = [20 30 40 60 50 80 90 95 99];
mp_support('bl_display_final') = true;
mp_support('bl_display_detail') = false;    
mp_support('bl_display_drvm2outcomes') = false;  
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;

% Call Function     
mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);
tb_dist_stats_all = mp_cl_mt_xyz_of_s('tb_outcomes');
% Select columns
tb_dist_stats_all_save = tb_dist_stats_all(1:3,:);
ar_st_columns = ["coefofvar", "gini", "varianceoflog", ...
				 "p99p50ratio", "p90p50ratio", "meantomedian", "p50p30ratio", ...
				 "fracP0toP20", "fracP20toP40", "fracP40toP60", "fracP60toP80", "fracP80toP100", ...
				 "fracP90toP95", "fracP95toP99", "fracP99toP100"];
				 
varianceoflog = tb_dist_stats_all{4:6,"sd"}.^2;

p99p50ratio = tb_dist_stats_all_save{:,"p99"}./tb_dist_stats_all_save{:,"p50"};
p90p50ratio = tb_dist_stats_all_save{:,"p90"}./tb_dist_stats_all_save{:,"p50"};
meantomedian = tb_dist_stats_all_save{:,"mean"}./tb_dist_stats_all_save{:,"p50"};
p50p30ratio = tb_dist_stats_all_save{:,"p50"}./tb_dist_stats_all_save{:,"p30"};
fracP0toP20 = tb_dist_stats_all_save{:,"fracByP20"};
fracP20toP40 = tb_dist_stats_all_save{:,"fracByP40"} - tb_dist_stats_all_save{:,"fracByP20"};
fracP40toP60 = tb_dist_stats_all_save{:,"fracByP60"} - tb_dist_stats_all_save{:,"fracByP40"};
fracP60toP80 = tb_dist_stats_all_save{:,"fracByP80"} - tb_dist_stats_all_save{:,"fracByP60"};
fracP80toP100 = 1 - tb_dist_stats_all_save{:,"fracByP80"};

fracP90toP95 = tb_dist_stats_all_save{:,"fracByP95"} - tb_dist_stats_all_save{:,"fracByP90"};
fracP95toP99 = tb_dist_stats_all_save{:,"fracByP99"} - tb_dist_stats_all_save{:,"fracByP95"};
fracP99toP100 = 1 - tb_dist_stats_all_save{:,"fracByP99"};

tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, varianceoflog, 'Before', 'gini');
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p99p50ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p90p50ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, meantomedian);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p50p30ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP0toP20);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP20toP40);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP40toP60);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP60toP80);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP80toP100);

tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP90toP95);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP95toP99);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP99toP100);
disp(tb_dist_stats_all_save(:, ar_st_columns));
% Core Stats Table
if (bl_save_csv)
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_dist_stats_all_save(:, ar_st_columns), [spt_simu_results_csv 'stats_all_allages_vrrcore.csv'], 'WriteRowNames', true);
end
%% 
% Statistics overall distributionally for 18 to 64 year olds.

% construct input data
income_grp = income(min_age:max_age, :, :, : ,: ,:);
earning_grp = earning(min_age:max_age, :, :, : ,: ,:);
wealth_grp = wealth(min_age:max_age, :, :, : ,: ,:);
Phi_true_grp = Phi_true_1(min_age:max_age, :, :, : ,: ,:);

mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('income') = {income_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('earning') = {earning_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('wealth') = {wealth_grp(:), zeros(1)};
mp_cl_ar_xyz_of_s('earninglog') = {log(earning_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('incomelog') = {log(income_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('wealthlog') = {log(wealth_grp(:)), zeros(1)};
mp_cl_ar_xyz_of_s('ar_st_y_name') = ["earning", "income", "wealth", "earninglog", "incomelog", "wealthlog"];

% controls
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('ar_fl_percentiles') = [20 30 40 60 50 80 90 95 99];
mp_support('bl_display_final') = true;
mp_support('bl_display_detail') = false;    
mp_support('bl_display_drvm2outcomes') = false;  
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;

% Call Function     
mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_grp(:)/sum(Phi_true_grp,'all'), mp_cl_ar_xyz_of_s, mp_support);
tb_dist_stats_all = mp_cl_mt_xyz_of_s('tb_outcomes');
% Select columns
tb_dist_stats_all_save = tb_dist_stats_all(1:3,:);
ar_st_columns = ["coefofvar", "gini", "varianceoflog", ...
				 "p99p50ratio", "p90p50ratio", "meantomedian", "p50p30ratio", ...
				 "fracP0toP20", "fracP20toP40", "fracP40toP60", "fracP60toP80", "fracP80toP100", ...
				 "fracP90toP95", "fracP95toP99", "fracP99toP100"];
				 
varianceoflog = tb_dist_stats_all{4:6,"sd"}.^2;

p99p50ratio = tb_dist_stats_all_save{:,"p99"}./tb_dist_stats_all_save{:,"p50"};
p90p50ratio = tb_dist_stats_all_save{:,"p90"}./tb_dist_stats_all_save{:,"p50"};
meantomedian = tb_dist_stats_all_save{:,"mean"}./tb_dist_stats_all_save{:,"p50"};
p50p30ratio = tb_dist_stats_all_save{:,"p50"}./tb_dist_stats_all_save{:,"p30"};
fracP0toP20 = tb_dist_stats_all_save{:,"fracByP20"};
fracP20toP40 = tb_dist_stats_all_save{:,"fracByP40"} - tb_dist_stats_all_save{:,"fracByP20"};
fracP40toP60 = tb_dist_stats_all_save{:,"fracByP60"} - tb_dist_stats_all_save{:,"fracByP40"};
fracP60toP80 = tb_dist_stats_all_save{:,"fracByP80"} - tb_dist_stats_all_save{:,"fracByP60"};
fracP80toP100 = 1 - tb_dist_stats_all_save{:,"fracByP80"};

fracP90toP95 = tb_dist_stats_all_save{:,"fracByP95"} - tb_dist_stats_all_save{:,"fracByP90"};
fracP95toP99 = tb_dist_stats_all_save{:,"fracByP99"} - tb_dist_stats_all_save{:,"fracByP95"};
fracP99toP100 = 1 - tb_dist_stats_all_save{:,"fracByP99"};

tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, varianceoflog, 'Before', 'gini');
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p99p50ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p90p50ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, meantomedian);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, p50p30ratio);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP0toP20);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP20toP40);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP40toP60);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP60toP80);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP80toP100);

tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP90toP95);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP95toP99);
tb_dist_stats_all_save = addvars(tb_dist_stats_all_save, fracP99toP100);
disp(tb_dist_stats_all_save(:, ar_st_columns));
% Core Stats Table
if (bl_save_csv)
    mp_path = snw_mp_path('fan');
    spt_simu_results_csv = mp_path('spt_simu_results_csv');
    writetable(tb_dist_stats_all_save(:, ar_st_columns), [spt_simu_results_csv 'stats_all_18t64_vrrcore.csv'], 'WriteRowNames', true);
end