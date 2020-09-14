%% 2019 Age, Income, Kids, Marry MPC given Actual and Optimal Allocations
% Given the simulation results and allocation results, load in allocation CSV 
% files, and analyze not the instantaneous marginal propensity to consume, which 
% is an output of the dynamic programming problem, but the realized marginal propensity 
% to consume given allocations assigned to each individual. This marginal propensity 
% distribution shifts depending on the allocation scheme. 
% 
% Under different allocation schemes, recipients of allocations differ, and 
% how much each receives differs. The marginal propensity to consume given allocations 
% can only be computed for individuals who receive allocations under a particular 
% scheme. 
%% Common Parameters

clc;
clear all;
% Paths
mp_paths = snw_mp_path('fan');
spt_simu_results_csv = mp_paths('spt_simu_results_csv');
spt_b1_a64_mulone_18t64_feasible = fullfile(spt_simu_results_csv, 'b1_a64_mulone_18t64', 'df_alloc_all_feasible_mulone.csv');
% Distributional Controls
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('bl_display_detail') = false;
mp_support('bl_display_final') = false;
mp_support('bl_display_drvm2outcomes') = false;
mp_support('ar_fl_percentiles') = [1:1:99];
ar_st_pcols = strrep(string(strcat('p', num2str(mp_support('ar_fl_percentiles')'))),' ', '')';
ar_st_pcols_subset = strrep(string(strcat('p', num2str([1 5 10 25 50 75 90 99]'))),' ', '')';
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;
%% Load in Actual Allocation File and Test
% For the Average Propensity given actual allocations, can use csv from any 
% of the allocation folders, they all share the same actual allocation MPC column. 
% But each allocation file has its own distribution of APC. 
% 
% Actually the actual allocation column does differ for G4 vs G47 vs feasible. 
% The one that is the closest to reality is the feasible allocation result. 
% 
% There are 17 allocation results presented in Table 1 of of the main text of 
% the paper. Each figure presents the aggregate marginal propensity to consume 
% across households. But we also know the distribution for each, so we can present 
% that. 

% get data
tb_b1_a64_mulone_18t64_feasible = readtable(spt_b1_a64_mulone_18t64_feasible, ...
    'PreserveVariableNames', true, 'TreatAsEmpty',{'NA'});
mt_b1_a64_mulone_18t64_feasible = rmmissing(tb_b1_a64_mulone_18t64_feasible{:,["mass", "apc_actual_1st"]});
% Get columns
apc_actual_1st_feasible_hh_mass = mt_b1_a64_mulone_18t64_feasible(:,1);
apc_actual_1st_feasible_hh_mass = apc_actual_1st_feasible_hh_mass/sum(apc_actual_1st_feasible_hh_mass(:));
apc_actual_1st_feasible = mt_b1_a64_mulone_18t64_feasible(:,2);
% generate distribution
mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('apc_actual_1st_feasible') = {apc_actual_1st_feasible(:), zeros(1)};
mp_cl_ar_xyz_of_s('ar_st_y_name') = ["apc_actual_1st_feasible"];
% Call Function
mp_cl_mt_xyz_of_s_out = ff_simu_stats(apc_actual_1st_feasible_hh_mass(:), mp_cl_ar_xyz_of_s, mp_support);
tb_stats_percentiles_all = mp_cl_mt_xyz_of_s_out('tb_outcomes');
% Modify Table
cl_st_varrownames = tb_stats_percentiles_all.Properties.RowNames;
csv_col_name = string(cl_st_varrownames);
tb_stats_percentiles_all = addvars(tb_stats_percentiles_all, csv_col_name, 'Before', 1);
csv = ["df_alloc_all_feasible_mulone.csv"];
tb_stats_percentiles_all = addvars(tb_stats_percentiles_all, csv, 'Before', 1);
folder = ["b1_a64_mulone_18t64"];
tb_stats_percentiles_all = addvars(tb_stats_percentiles_all, folder, 'Before', 1);
% row name
tb_stats_percentiles_all.Properties.RowNames = ["row=1"];
% select columns
tb_stats_percentiles_all = tb_stats_percentiles_all(:, ["folder" "csv" "csv_col_name" "min" "max" ar_st_pcols]);
%% Allocation Results under Alternative Allocation Structures
% Load in all key optimal allocation results under alternative constraints

% Allocation MPC Column 
svr_mpc_select = 'apc_optiexpmin_1st';
st_common_prefix = 'b1_a64';
% Load in Allocation Results
it_row_ctr = 1;
for it_year_bounds=1
    if (it_year_bounds == 1)
        st_yr_bounds = '_18t64';
    elseif (it_year_bounds == 2)
        st_yr_bounds = '_18t99';
    end
    for it_alloc=[1 2:8 10:11 12:18]
%     for it_alloc=[1 2:6 8 10 12:16 18]
        it_row_ctr = it_row_ctr + 1;
        
        % Construct name components
        if (it_alloc >= 1 && it_alloc <= 8)
            st_age_type = '_feasible';
        end
        
        if ((it_alloc == 10) || (it_alloc >= 12 && it_alloc <= 18))
            st_age_type = '_optimalg4';
        end
        if (it_alloc == 11)
            st_age_type = '_optimalg47';
        end
        
        % Bound types
        if (it_alloc == 1 || it_alloc == 10 || it_alloc == 11)
            st_bound_type = '_mulone';
        end
        if (it_alloc == 2 || it_alloc == 12)
            st_bound_type = '_dbladt';
        end
        if (it_alloc == 3 || it_alloc == 13)
            st_bound_type = '_trbadt';
        end
        if (it_alloc == 4 || it_alloc == 14)
            st_bound_type = '_dblkid';
        end
        if (it_alloc == 5 || it_alloc == 15)
            st_bound_type = '_trbkid';
        end
        if (it_alloc == 6 || it_alloc == 16)
            st_bound_type = '_dblbth';
        end
        if (it_alloc == 7 || it_alloc == 17)
            st_bound_type = '_trbbth';
        end
        if (it_alloc == 8 || it_alloc == 18)
            st_bound_type = '_20k';
        end
        
        % get file path
        srn_folder = [st_common_prefix st_bound_type st_yr_bounds];
        snm_csv_file = ['df_alloc_all' st_age_type st_bound_type '.csv'];
        spt_file = fullfile(spt_simu_results_csv, srn_folder, snm_csv_file);
        st_name = [st_common_prefix st_bound_type st_yr_bounds st_bound_type];
        
        % Get Data
        tb_read = readtable(spt_file, 'PreserveVariableNames', true, 'TreatAsEmpty',{'NA'});
        mt_read = rmmissing(tb_read{:,["mass", svr_mpc_select]});
                
        % Get columns
        ar_mass = mt_read(:,1);
        ar_mass = ar_mass/sum(ar_mass(:));
        ar_apc = mt_read(:,2);
        % generate distribution
        mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
        mp_cl_ar_xyz_of_s(svr_mpc_select) = {ar_apc(:), zeros(1)};
        mp_cl_ar_xyz_of_s('ar_st_y_name') = [string(svr_mpc_select)];
        % Call Function
        mp_cl_mt_xyz_of_s_out = ff_simu_stats(ar_mass(:), mp_cl_ar_xyz_of_s, mp_support);
        tb_stats_percentiles = mp_cl_mt_xyz_of_s_out('tb_outcomes');        
        
        % Modify table for concatenation        
        cl_st_varrownames = tb_stats_percentiles.Properties.RowNames;
        csv_col_name = string(cl_st_varrownames);
        tb_stats_percentiles = addvars(tb_stats_percentiles, csv_col_name, 'Before', 1);
        csv = string(snm_csv_file);
        tb_stats_percentiles = addvars(tb_stats_percentiles, csv, 'Before', 1);
        folder = string(srn_folder);
        tb_stats_percentiles = addvars(tb_stats_percentiles, folder, 'Before', 1);        
        % row name
        tb_stats_percentiles.Properties.RowNames = string(['row=' num2str(it_row_ctr)]);
        tb_stats_percentiles = tb_stats_percentiles(:, ["folder" "csv" "csv_col_name" "min" "max" ar_st_pcols]);
        % stack
        tb_stats_percentiles_all = [tb_stats_percentiles_all;tb_stats_percentiles];
    end
end
%% Save MPC given Allocation to File
% save to csv and select to keep only percentile columns

% write to table
spt_simu_results_csv = [mp_paths('spt_simu_results_csv') '/apc_dist_all_percentile.csv'];
writetable(tb_stats_percentiles_all, spt_simu_results_csv);
% write to table
spt_simu_results_csv = [mp_paths('spt_simu_results_csv') '/apc_dist_all_percentile_subset.csv'];
writetable(tb_stats_percentiles_all(:, ["folder" "csv" "csv_col_name" "min" "max" ar_st_pcols_subset])...
    , spt_simu_results_csv);