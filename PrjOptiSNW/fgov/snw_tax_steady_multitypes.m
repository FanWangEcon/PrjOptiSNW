%% SNW_TAX_STEADY_MULTITYPES Solves Budget Clearning Tax Heterogeneous Beta Edu
%    There are four possible types of households, high and low discount
%    factor type, and high and low education. There are exogenously
%    determined proobabilities for high and low education types, and the
%    fraction of low beta type households in high and low education types.
%
%    The code can be used to solve for single beta type model as well, if
%    simply set the low and high beta to the same value.
%
%    Heterogeneous type parameters:
%
%    * BETA discount
%    * THETA total factor productivity normalizer
%    * R interest rate
%
%    Tax parameter:
%
%    * AGRID asset grid
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC(MP_PARAMS) invoke model with externally set
%    parameter map MP_PARAMS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC(MP_PARAMS, MP_CONTROLS) invoke model with
%    externally set parameter map MP_PARAMS as well as control mpa
%    MP_CONTROLS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC(MP_PARAMS, MP_CONTROLS, V_VFI_FIX) provides
%    existing value function. Suppose there is sudden shock, but future
%    value is preserved after one period. So now we have new value that is
%    specific to this period, that is the output V_VFI, the input V_VFI_FIX
%    is the value for all future periods. When this program is called with
%    V_VFI_FIX, the resource equation will use the unemployment shock
%    information.
%
%    See also SNWX_VFI_MAIN, SNW_MP_CONTROL, SNW_MP_PARAM
%

function varargout = snw_tax_steady_multitypes(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))
    
    a2_init = 1.2445;
    theta_current=0.565228521783443;
    mp_controls_ext = 0;
    if (length(varargin)==1)
        [st_param_group_base] = varargin{:};
    elseif (length(varargin)==3)
        [st_param_group_base, a2_init, theta_current] = varargin{:};
    elseif (length(varargin)==4)
        [st_param_group_base, a2_init, theta_current, mp_controls_ext] = varargin{:};
    elseif (length(varargin)==9)
        [st_param_group_base, a2_init, theta_current, mp_controls_ext, ...
            ls_fl_beta_val, ls_it_edu_simu_type, ...
            fl_wgt_college_frac, fl_wgt_low_beta_college, fl_wgt_low_beta_non_college] = varargin{:};
    elseif (length(varargin)==11)
        [st_param_group_base, a2_init, theta_current, mp_controls_ext, ...
            ls_fl_beta_val, ls_it_edu_simu_type, ...
            fl_wgt_college_frac, fl_wgt_low_beta_college, fl_wgt_low_beta_non_college, ...
            bl_print_calibrate_multitypes, bl_print_calibrate_multitypes_verbose] = varargin{:};
    else
        error('Need to provide 3 parameter inputs');
    end
    
else
    
    clc;
    clear all;
    
    % Defaults
    %     st_param_group_base = 'default_base';
    st_param_group_base = 'default_tiny';
    
    % beta and edu types
    %     ls_fl_beta_val = [0.60, 0.95];
    ls_fl_beta_val = [0.971162552785405, 0.971162552785405];
    
    % 1 is low education, 2 is high education
    ls_it_edu_simu_type = [1, 2];
    
    fl_wgt_college_frac = 0.303;
    fl_wgt_low_beta_college =  0.1;
    fl_wgt_low_beta_non_college = 0.4;
    
    %     a2_init = 1.528571486531964;
    a2_init = 1.2445;
    theta_current=0.565228521783443;
    
    bl_print_calibrate_multitypes = true;
    bl_print_calibrate_multitypes_verbose = true;
    mp_controls_ext = 0;
end

%% Parameters

fl_solu_tax_tol = 1e-4;
fl_solu_tax_prime_tol = 1e-4;
it_solu_tax_iter_max = 3;
it_solu_tax_prime_iter_max = 100;

% st_param_group_base = 'default_tiny';
%         st_param_group_base = 'default_docdense';
% st_param_groups_base = 'default_moredense_a65zh266zs5';

%% Parse Model Parameters Relevant Across Solutions
% st_param_group = 'default_moredense_a65zh266zs5';
mp_params = snw_mp_param(st_param_group_base);
params_group = values(mp_params, {'theta', 'r', 'g_n', 'g_cons', 'a2','jret'});
[theta, r, g_n, g_cons, a2_stored, jret] = params_group{:};
params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

% tax rate is equal to 1.528
% a2=1.528571486531964;
a2_solu_current = a2_init;

%% Store
fl_err_rev_iter_gap=Inf;

% store
mt_a2_solu_prime_store = NaN([it_solu_tax_prime_iter_max+1, it_solu_tax_iter_max+1]);
mt_fl_err_rev = NaN([it_solu_tax_prime_iter_max+1, it_solu_tax_iter_max+1]);

it_solu_tax_iter = 0;
while fl_err_rev_iter_gap > fl_solu_tax_tol && it_solu_tax_iter <= it_solu_tax_iter_max
    it_solu_tax_iter = it_solu_tax_iter + 1;
    
    %% Solve Model Given Tax
    % store results collections
    cl_mp_params = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    cl_mp_controls = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    cl_ap_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    cl_cons_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    cl_y_all_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    cl_Phi_true = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    
    % Aggregate outcome stores
    mt_Y_inc_median = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    mt_A_agg = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    mt_Y_inc_ori_agg = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
    
    % loop and solve
    it_beta_ctr = 0;
    for fl_beta_val = ls_fl_beta_val
        it_beta_ctr = it_beta_ctr + 1;
        
        it_edu_ctr = 0;
        for it_edu_simu_type = ls_it_edu_simu_type
            it_edu_ctr = it_edu_ctr + 1;
            
            % parameter method 1
            [mp_params, mp_controls] = get_param_maps(st_param_group_base, it_edu_simu_type, ...
            fl_beta_val, a2_solu_current, theta_current, mp_controls_ext);

%             % Set Up Parameters
%             mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
%             if (it_edu_simu_type == 1)
%                 mp_more_inputs('st_edu_simu_type') = 'low';
%                 st_param_group_suffix = '_e1lm2';
%             elseif (it_edu_simu_type == 2)
%                 mp_more_inputs('st_edu_simu_type') = 'high';
%                 st_param_group_suffix = '_e2hm2';
%             end
%             % param group name
%             st_param_group = [st_param_group_base st_param_group_suffix];
%             % get parametesr
%             mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
%             mp_controls = snw_mp_control('default_test');
%             % beta parameter update
%             mp_params('beta') = fl_beta_val;
%             % Tax in steady state world, don't need to define covid tax
%             mp_params('a2') = a2_solu_current;
%             mp_params('theta') = theta_current;
%             % override
%             if mp_controls_ext ~= 0
%                 mp_controls = [mp_controls; mp_controls_ext];
%             end
            
            % Solve for distributions and policies
            [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
            % solve for distributions
            [Phi_true, ~, A_agg, Y_inc_agg, ~, mp_dsvfi_results] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
            
            mp_cl_mt_xyz_of_s = mp_dsvfi_results('mp_cl_mt_xyz_of_s');
            tb_outcomes = mp_cl_mt_xyz_of_s('tb_outcomes');
            % Stats stored from snw_ds_main_vec
            Y_inc_median = tb_outcomes{'y_all', 'p50'};
            
            % Compute the same statistics here. (debug code)
            y_all_ss = mp_dsvfi_results('y_all_ss');
            if (bl_print_calibrate_multitypes_verbose)
                mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
                mp_cl_ar_xyz_of_s('y_all_ss') = {y_all_ss(:), zeros(1)};
                mp_cl_ar_xyz_of_s('ar_st_y_name') = ["y_all_ss"];
                mp_support = containers.Map('KeyType','char', 'ValueType','any');
                mp_support('ar_fl_percentiles') = [10 20 25 50 75 80 90];
                mp_support('bl_display_final') = false;
                mp_support('bl_display_detail') = false;
                mp_support('bl_display_drvm2outcomes') = false;
                mp_support('bl_display_drvstats') = false;
                mp_support('bl_display_drvm2covcor') = false;
                mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true(:)/sum(Phi_true,'all'), ...
                    mp_cl_ar_xyz_of_s, mp_support);
                tb_outcomes = mp_cl_mt_xyz_of_s('tb_outcomes');
                Y_inc_median_here = tb_outcomes{'y_all_ss', 'p50'};
                
                if (Y_inc_median_here ~= Y_inc_median)
                    error('error Y_inc_median_here ~= Y_inc_median');
                end
            end
            
            % Store results
            cl_mp_params{it_beta_ctr, it_edu_ctr} = mp_params;
            cl_mp_controls{it_beta_ctr, it_edu_ctr} = mp_controls;
            cl_ap_ss{it_beta_ctr, it_edu_ctr} = ap_ss;
            cl_cons_ss{it_beta_ctr, it_edu_ctr} = cons_ss;
            cl_y_all_ss{it_beta_ctr, it_edu_ctr} = y_all_ss;
            cl_Phi_true{it_beta_ctr, it_edu_ctr} = Phi_true;
            
            mt_Y_inc_median(it_beta_ctr, it_edu_ctr) = Y_inc_median;
            mt_A_agg(it_beta_ctr, it_edu_ctr) = A_agg;
            mt_Y_inc_ori_agg(it_beta_ctr, it_edu_ctr) = Y_inc_agg;
            
        end
    end
    
    %% Given P(tax), compute tax' that clears budget given P(tax)
    a2 = a2_solu_current;
    
    fl_err_rev_iter_gap=Inf;
    
    it_solu_tax_prime_iter=1;
    mt_a2_solu_prime_store(1, it_solu_tax_iter) = a2;
    mt_fl_err_rev(1, it_solu_tax_iter) = fl_err_rev_iter_gap;
    while fl_err_rev_iter_gap > fl_solu_tax_prime_tol && ...
            it_solu_tax_prime_iter <= it_solu_tax_prime_iter_max
        it_solu_tax_prime_iter = it_solu_tax_prime_iter + 1;
        
        % store results collections
        mt_SS_spend = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
        mt_Y_inc_agg = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
        mt_Tax_revenues_aux = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
        mt_Bequests_aux = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
        mt_grp_wgt = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
        
        % Solve and iterate over tax revenue
        it_beta_ctr = 0;
        for fl_beta_val = ls_fl_beta_val
            it_beta_ctr = it_beta_ctr + 1;
            
            it_edu_ctr = 0;
            for it_edu_simu_type = ls_it_edu_simu_type
                it_edu_ctr = it_edu_ctr + 1;
                               
                ap_ss = cl_ap_ss{it_beta_ctr, it_edu_ctr};
                cons_ss = cl_cons_ss{it_beta_ctr, it_edu_ctr};
                Phi_true = cl_Phi_true{it_beta_ctr, it_edu_ctr};
               
                % parameter method
                [mp_params, mp_controls] = get_param_maps(st_param_group_base, it_edu_simu_type, ...
                fl_beta_val, a2, theta_current, mp_controls_ext);
                
                % update tax
                mp_params('a2') = a2;
                
                % solve aggregation
                it_tax_iter_max_inner = 0;
                fl_rev_equal_tax_tol = fl_solu_tax_prime_tol^2;
                [SS_spend, Y_inc_agg, Tax_revenues_aux, Bequests_aux, it_hat_solu_tax_iter_inner] = ...
                    snw_ds_aggregation(mp_params, mp_controls, ap_ss, cons_ss, Phi_true, ...
                    it_tax_iter_max_inner, fl_rev_equal_tax_tol);
                
                % Construct weights
                if it_edu_ctr == 1 && it_beta_ctr == 1
                    fl_wgt = (1-fl_wgt_college_frac)*fl_wgt_low_beta_non_college;
                elseif it_edu_ctr == 1 && it_beta_ctr == 2
                    fl_wgt = (1-fl_wgt_college_frac)*(1-fl_wgt_low_beta_non_college);
                elseif it_edu_ctr == 2 && it_beta_ctr == 1
                    fl_wgt = fl_wgt_college_frac*fl_wgt_low_beta_college;
                elseif it_edu_ctr == 2 && it_beta_ctr == 2
                    fl_wgt = fl_wgt_college_frac*(1-fl_wgt_low_beta_college);
                end
                
                % Store results
                mt_SS_spend(it_beta_ctr, it_edu_ctr) = SS_spend;
                mt_Y_inc_agg(it_beta_ctr, it_edu_ctr) = Y_inc_agg;
                mt_Tax_revenues_aux(it_beta_ctr, it_edu_ctr) = Tax_revenues_aux;
                mt_Bequests_aux(it_beta_ctr, it_edu_ctr) = Bequests_aux;
                mt_grp_wgt(it_beta_ctr, it_edu_ctr) = fl_wgt;
                
            end
        end
        
        % Aggregation reweighting
        SS_spend_wgt = sum(mt_SS_spend.*mt_grp_wgt, 'all');
        SS_Y_inc_agg_wgt = sum(mt_Y_inc_agg.*mt_grp_wgt, 'all');
        SS_Tax_revenues_aux_wgt = sum(mt_Tax_revenues_aux.*mt_grp_wgt, 'all');
        SS_mt_Bequests_aux_wgt = sum(mt_Bequests_aux.*mt_grp_wgt, 'all');
        
        % Find new tax
        if bequests_option==1
            if throw_in_ocean==0
                SS_Tax_revenues_aux_wgt=SS_Tax_revenues_aux_wgt+SS_mt_Bequests_aux_wgt*(1+r);
            end
        end
        
        a2=a2*(((SS_spend_wgt+g_cons*SS_Y_inc_agg_wgt)/SS_Tax_revenues_aux_wgt)^0.75); % Find value of a2 that balances government budget
        
        fl_err_rev_iter_gap=abs((SS_Tax_revenues_aux_wgt/(SS_spend_wgt+g_cons*SS_Y_inc_agg_wgt))-1);
        
        if (bl_print_calibrate_multitypes_verbose)
            st_tax_iter = strjoin(...
                ["SNW_TAX_STEADY_MULTITYPES tax and spend", ...
                ['it=' num2str(it_solu_tax_prime_iter)], ...
                ['err=' num2str(fl_err_rev_iter_gap)] ...
                ], ";");
            disp(st_tax_iter);
        end
        
        % Track/store
        mt_a2_solu_prime_store(it_solu_tax_prime_iter, it_solu_tax_iter) = a2;
        mt_fl_err_rev(it_solu_tax_prime_iter, it_solu_tax_iter) = fl_err_rev_iter_gap;
    end
    
    fl_err_rev_iter_gap = Inf;
    a2_solu_current = a2;
    
end

%% Show Results
mt_fl_err_rev = mt_fl_err_rev(:,1:it_solu_tax_iter);
mt_fl_err_rev = mt_fl_err_rev(sum(~isnan(mt_fl_err_rev)==0,2)<it_solu_tax_iter,:);

mt_a2_solu_prime_store = mt_a2_solu_prime_store(:,1:it_solu_tax_iter);
mt_a2_solu_prime_store = mt_a2_solu_prime_store(sum(~isnan(mt_a2_solu_prime_store)==0,2)<it_solu_tax_iter,:);

if (bl_print_calibrate_multitypes)
    disp('');
    disp('======================================================================================');
    disp(['theta normalizer=' num2str(theta_current)]);
    disp('======================================================================================');
    disp('mt_fl_err_rev: abs((SS_Tax_revenues_aux_wgt/(SS_spend_wgt+g_cons*SS_Y_inc_agg_wgt))-1)');
    disp('triple nested GE calibrator: (gov clear) row=third-nest-holding-policy-iterate-tax');
    disp('triple nested GE calibrator: (policy converge) col=second-nest-holding-normalizer-iterate-policy');
    disp('triple nested GE calibrator: (match moments) matrix=vary-normalizer-and-possibly-other-params');
    disp(mt_fl_err_rev);
    disp('tax a2 updates: a2=a2*(((SS_spend_wgt+g_cons*SS_Y_inc_agg_wgt)/SS_Tax_revenues_aux_wgt)^0.75)');
    disp('triple nested GE calibrator: (gov clear) row=third-nest-holding-policy-iterate-tax');
    disp('triple nested GE calibrator: (policy converge) col=second-nest-holding-normalizer-iterate-policy');
    disp('triple nested GE calibrator: (match moments) matrix=vary-normalizer-and-possibly-other-params');
    disp(mt_a2_solu_prime_store);
    disp('======================================================================================');
    disp('');
end

%% Compute Distributional Statistics for All Groups given Converged policy and tax
% Average Output and Savings
A_agg_wgt_mean = sum(mt_A_agg.*mt_grp_wgt, 'all');
Y_inc_agg_wgt_mean = sum(mt_Y_inc_agg.*mt_grp_wgt, 'all');

% Median Income
y_all_ss_alltypes = [];
Phi_true_wgt_flat_alltypes = [];

it_beta_ctr = 0;
for fl_beta_val = ls_fl_beta_val
    it_beta_ctr = it_beta_ctr + 1;
    
    it_edu_ctr = 0;
    for it_edu_simu_type = ls_it_edu_simu_type
        it_edu_ctr = it_edu_ctr + 1;
        
        y_all_ss = cl_y_all_ss{it_beta_ctr, it_edu_ctr};
        Phi_true = cl_Phi_true{it_beta_ctr, it_edu_ctr};
        
        %         if it_edu_ctr == 1 && it_beta_ctr == 1
        %             fl_wgt = (1-fl_wgt_college_frac)*fl_wgt_low_beta_non_college;
        %         elseif it_edu_ctr == 1 && it_beta_ctr == 2
        %             fl_wgt = (1-fl_wgt_college_frac)*(1-fl_wgt_low_beta_non_college);
        %         elseif it_edu_ctr == 2 && it_beta_ctr == 1
        %             fl_wgt = fl_wgt_college_frac*fl_wgt_low_beta_college;
        %         elseif it_edu_ctr == 2 && it_beta_ctr == 2
        %             fl_wgt = fl_wgt_college_frac*(1-fl_wgt_low_beta_college);
        %         end
        
        fl_wgt_dup = mt_grp_wgt(it_beta_ctr, it_edu_ctr);
        
        Phi_true_wgt = (Phi_true(:)/sum(Phi_true,'all'))*fl_wgt;
        Phi_true_wgt_flat = Phi_true_wgt(:);
        
        y_all_ss_alltypes = [y_all_ss_alltypes y_all_ss(:)'];
        Phi_true_wgt_flat_alltypes = [Phi_true_wgt_flat_alltypes Phi_true_wgt_flat'];
        
    end
end

mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
mp_cl_ar_xyz_of_s('y_all_ss_alltypes') = {y_all_ss_alltypes(:), zeros(1)};
mp_cl_ar_xyz_of_s('ar_st_y_name') = ["y_all_ss_alltypes"];
mp_support = containers.Map('KeyType','char', 'ValueType','any');
mp_support('ar_fl_percentiles') = [10 20 25 50 75 80 90];
mp_support('bl_display_final') = false;
mp_support('bl_display_detail') = false;
mp_support('bl_display_drvm2outcomes') = false;
mp_support('bl_display_drvstats') = false;
mp_support('bl_display_drvm2covcor') = false;

% Call Function
mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true_wgt_flat_alltypes(:)/sum(Phi_true_wgt_flat_alltypes,'all'), ...
    mp_cl_ar_xyz_of_s, mp_support);
tb_outcomes = mp_cl_mt_xyz_of_s('tb_outcomes');
Y_inc_agg_wgt_median = tb_outcomes{'y_all_ss_alltypes', 'p50'};

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = a2;
    elseif (it_k==2)
        ob_out_cur = A_agg_wgt_mean;
    elseif (it_k==3)
        ob_out_cur = Y_inc_agg_wgt_mean;
    elseif (it_k==4)
        ob_out_cur = Y_inc_agg_wgt_median;
    elseif (it_k==5)
        ob_out_cur = mt_a2_solu_prime_store;
    elseif (it_k==6)
        ob_out_cur = mt_fl_err_rev;
    end
    varargout{it_k} = ob_out_cur;
end

end

% Parameter getter
function [mp_params, mp_controls] = get_param_maps(st_param_group_base, it_edu_simu_type, ...
    fl_beta_val, a2_solu_current, theta_current, mp_controls_ext)
% Separate function to avoid re-typing

% Set Up Parameters
mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
if (it_edu_simu_type == 1)
    mp_more_inputs('st_edu_simu_type') = 'low';
    st_param_group_suffix = '_e1lm2';
elseif (it_edu_simu_type == 2)
    mp_more_inputs('st_edu_simu_type') = 'high';
    st_param_group_suffix = '_e2hm2';
end
% param group name
st_param_group = [st_param_group_base st_param_group_suffix];
% get parametesr
mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);
mp_controls = snw_mp_control('default_test');
% beta parameter update
mp_params('beta') = fl_beta_val;
% Tax in steady state world, don't need to define covid tax
mp_params('a2') = a2_solu_current;
mp_params('theta') = theta_current;
% override
if mp_controls_ext ~= 0
    mp_controls = [mp_controls; mp_controls_ext];
end

end