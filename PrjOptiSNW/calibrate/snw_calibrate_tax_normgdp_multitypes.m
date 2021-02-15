%% SNW Calibrate Beta and Normalize 
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/PrjOptiSNW/params/snw_ds_main.m 
% *snw_ds_main*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*,> this function calibrates the discount factor and also solves for 
% the normalizing constant. 

%% Calibrate Parameter Controls for SNW Functions
% Set up controls for shock process and tiny/small/dense/densemore
clc;
clear all;

mp_paths = snw_mp_path('fan');
spt_simu_outputs_log = mp_paths('spt_simu_outputs_log');

% Defaults
% st_shock_method = 'tauchen';
% st_param_group_base = 'default_tiny';
% st_param_group_base = 'default_base';
% st_param_group_base = 'default_dense';
st_param_group_base = 'default_docdense';
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

% a2_init = 1.528571486531964;
% a2_init = 1.0849;
% a2_init = 1.2445;
% a2_init = 1.2026;
a2_init = 1.6919;

% theta_init=0.565228521783443;
% theta_init=0.6027;
theta_init=0.55912;

bl_print_calibrate_multitypes = true;
bl_print_calibrate_multitypes_verbose = false;

% start diary
spn_diary_filename = fullfile(spt_simu_outputs_log, ...
    ['cali_tax_normgdp_multitypes_' st_param_group_base '.out'] );
diary(spn_diary_filename);
%% key info
disp(['ls_fl_beta_val=' num2str(ls_fl_beta_val)]);


%% 
% Set up print defaults
% mp_params = containers.Map('KeyType', 'char', 'ValueType', 'any');
% mp_params = snw_mp_param(st_param_group_base);
mp_controls = snw_mp_control('default_test');
mp_controls('bl_timer') = true;
mp_controls('bl_print_vfi') = false;
mp_controls('bl_print_vfi_verbose') = false;
mp_controls('bl_print_ds') = false;
mp_controls('bl_print_ds_verbose') = false;
mp_controls('bl_print_ds_aggregation') = true;
mp_controls('bl_print_ds_aggregation_verbose') = false;

%% Calibrate Routine

%% Calibration
bl_check_theta_only = true;
err=1;
tol=0.005;
disp('Start calibration')
it_counter = 1;
a2_current = a2_init;
theta_current = theta_init;
while err>tol

    it=1;

    tic;

    [a2_new, A_agg_wgt_mean, Y_inc_agg_wgt_mean, ...
        Y_inc_agg_wgt_median, mt_a2_solu_prime_store, mt_fl_err_rev] = ...
        snw_tax_steady_multitypes(...
        st_param_group_base, ...
        a2_current, theta_current, ...
        mp_controls, ...
        ls_fl_beta_val, ls_it_edu_simu_type, ...
        fl_wgt_college_frac, fl_wgt_low_beta_college, fl_wgt_low_beta_non_college, ...
        bl_print_calibrate_multitypes, bl_print_calibrate_multitypes_verbose);

    disp(['a2_init:' num2str(a2_current) ', a2_new:' num2str(a2_new)])
    toc;
    a2_current = a2_new;
    
    % Comparison
    name='Median household income (target=1.0)=';
    name2=[name,num2str(Y_inc_agg_wgt_median)];
    disp(name2);
    name='Aggregate wealth to aggregate income (target=3.0)=';
    name2=[name,num2str(A_agg_wgt_mean/Y_inc_agg_wgt_mean)];
    disp(name2);
    
    err1=abs(Y_inc_agg_wgt_median-1.0); % Target: Median household income (normalized to 1 in the model)
    err2=abs((A_agg_wgt_mean/Y_inc_agg_wgt_mean)-3.0); % Target: Annual capital/income ratio of 3
    
    % Display A/Y ratio, but only check on normalization
    if bl_check_theta_only
        err=max(err1);
    else
        err=max(err1,err2);
    end
    
    % Beta and Theta
    theta = theta_current;
    param_update=[theta];
    if ~bl_check_theta_only
        beta = mp_params('beta');
        param_update=[theta;beta];
    end
    
    if err>tol    
        theta=theta*((1.0/Y_inc_agg_wgt_median)^0.2); % Normalize theta such that median household income equals 1
        if ~bl_check_theta_only
            beta=beta*((3.0/(A_agg_wgt_mean/Y_inc_agg_wgt_mean))^0.025); % Calibrate beta such that annual capital/income ratio equals 3
        end
    end
    
%     mp_params('theta') = theta;
    theta_current = theta;
    param_update=[param_update(1,1),theta];
    if ~bl_check_theta_only
        mp_params('beta') = beta;
        param_update=[param_update(1,1),theta;param_update(2,1),beta];
    end

    name='Old and updated value of theta=';
    name2=[name,num2str(param_update(1,:))];
    disp(name2);
    if ~bl_check_theta_only
        name='Old and updated value of beta=';
        name2=[name,num2str(param_update(2,:))];
        disp(name2);
    end

    it_counter = it_counter + 1;
    disp(['it_counter=' num2str(it_counter)]);
    disp([err1])
    
end

% 
diary('off');
