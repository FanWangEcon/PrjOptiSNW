%% SNW Calibrate Lock Down Consumption Effects
% Lock-down consumption effects

clc;
close all;
clear all;

%% Gamma 
fl_gamma = 1;

%% Parameters
for it_param_group=4:1:4
    for it_trumpbiden=1:1:1
        
        if (it_trumpbiden == 1)
            st_biden_or_trump = 'trumpchk';
        elseif (it_trumpbiden ==2)
            % Don't need to do this here. 
            st_biden_or_trump = 'bidenchk';
        end
        
        if (it_param_group == 1)
            st_param_group_base = 'default_tiny';
        elseif (it_param_group == 2)
            st_param_group_base = 'default_small53';
        elseif (it_param_group == 3)
            st_param_group_base = 'default_base';
        elseif (it_param_group == 4)
            st_param_group_base = 'default_dense';
        elseif (it_param_group == 5)
            st_param_group_base = 'default_docdense';
        end
        % st_param_group_base = 'default_moredense_a65zh21zs5_e2m2';
        % st_param_group_base = 'default_moredense_a65zh81zs5_e2m2';
        
        % beta and edu types
        ls_fl_beta_val = [0.60, 0.95];
        % ls_fl_beta_val = [0.971162552785405, 0.971162552785405];
        % 1 is low education, 2 is high education
        ls_it_edu_simu_type = [1, 2];
        
        fl_wgt_college_frac = 0.303;
        fl_wgt_low_beta_college =  0.1;
        fl_wgt_low_beta_non_college = 0.4;
        
        %% Path and Log
        mp_paths = snw_mp_path('fan');
        spt_simu_outputs_log = mp_paths('spt_simu_outputs_log');
        
        % start diary
        spn_diary_filename = fullfile(spt_simu_outputs_log, ...
            ['cali_lockdown_c_' st_param_group_base '_' ...
            st_biden_or_trump '_gamma' num2str(fl_gamma*100) '.out'] );
        diary(spn_diary_filename);
        
        %% Control Map
        % Set up print defaults
        % mp_params = containers.Map('KeyType', 'char', 'ValueType', 'any');
        % mp_params = snw_mp_param(st_param_group_base);
        mp_controls = snw_mp_control('default_test');
        mp_controls('bl_timer') = false;
        mp_controls('bl_print_vfi') = false;
        mp_controls('bl_print_vfi_verbose') = false;
        mp_controls('bl_print_ds') = false;
        mp_controls('bl_print_ds_verbose') = false;
        mp_controls('bl_print_ds_aggregation') = true;
        mp_controls('bl_print_ds_aggregation_verbose') = false;
        
        %% Array of Consumption drops to consider
        % Same min and max and grid points
        [fl_invbtlock_min, fl_invbtlock_max, it_invbtlock_points] = deal(0, 1, 20);
        st_grid_type = 'grid_powerspace';
        mp_grid_control = containers.Map('KeyType', 'char', 'ValueType', 'any');
        mp_grid_control('grid_powerspace_power') = 1.5;
        [ar_fl_invbtlock] = ff_saveborr_grid(fl_invbtlock_min, fl_invbtlock_max, it_invbtlock_points, ...
            st_grid_type, mp_grid_control);
        ar_fl_invbtlock = 1-ar_fl_invbtlock;
        
        
        %% Calibration
        disp('Start calibration')
        
        cl_mp_dsvfi_results_wgtbeta = cell(length(ar_fl_invbtlock), 1);
        ar_fl_cons_mean_betaedu_innerwgt = NaN(size(ar_fl_invbtlock));
        
        for it_invbtlock=1:length(ar_fl_invbtlock)
            
            % Parameter Map
            mp_params = containers.Map('KeyType','char', 'ValueType','any');
            mp_params('a2_covidyr') = 1.6996;
            % modify consumption losses today
            mp_params('invbtlock') = ar_fl_invbtlock(it_invbtlock);
            mp_params('gamma') = fl_gamma;
            % solve
            [~, ~, ~, ~, ~, mp_dsvfi_results_betaedu] = ...
                snw_ds_main_vec_multitypes(st_param_group_base, st_biden_or_trump, ...
                mp_params, mp_controls, ...
                ls_fl_beta_val, ls_it_edu_simu_type, ...
                fl_wgt_college_frac, fl_wgt_low_beta_college, fl_wgt_low_beta_non_college);
            
            %% display
            ff_container_map_display(mp_dsvfi_results_betaedu);
            
            % collect
            cl_mp_dsvfi_results_wgtbeta{it_invbtlock} = mp_dsvfi_results_betaedu;
            ar_fl_cons_mean_betaedu_innerwgt(it_invbtlock) = mp_dsvfi_results_betaedu('fl_cons_mean_betaedu_innerwgt');
        end
        
        %% Graph Results
        figure();
        hold on;
        ar_x = 1-ar_fl_invbtlock;
        ar_y = -(1-ar_fl_cons_mean_betaedu_innerwgt/ar_fl_cons_mean_betaedu_innerwgt(1));
        ar_x = ar_x*100;
        ar_y = ar_y*100;
        scatter(ar_x, ar_y, 300, [57 106 177]./255, 'd');
        line = plot(ar_x, ar_y);
        line.Color = [83 81 84]./255;
        line.LineStyle = '--';
        line.LineWidth = 2;
        
        % labeling
        title('Current Utility Lock Down');
        ylabel('Aggregate C Percent Reduction in Lockdown Year');
        xlabel('Lock Down Current Utility Percent/Proportional Reduction');
        grid on;
        
        %% Table Dispaly
        % Generate Table
        tb_show = array2table([ar_x,ar_y]);
        it_num_rows = length(ar_x);
        
        % Generate Row and Column Names
        cl_col_names = {'Lock Down Util Perc Reduce', 'Agg C Percent Reduc Lockdown Yr'};
        cl_row_names = strcat('lockdown_', string((1:it_num_rows)));
        
        tb_show.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
        tb_show.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
        disp(tb_show);
        
        % dairy off
        diary('off');
        
    end
end