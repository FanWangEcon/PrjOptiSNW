%% SNW_DS_MAIN_VEC_MULTITYPES Simulates Distributions Joint Beta Edu Types
%    For the version of the model with 2 beta types and 2 education types,
%    how do we analyze the distributional outcomes? We want to solve for
%    the beta types as well as the education types, combine the
%    distributions induced by induced by the four separate policy
%    functions, and generate a joint distributional that properly accounts
%    for different policy functions under beta/edu heterogeneities, and the
%    jointly induced distributional outcomes. 
%
%    The simple version for this is straight-forward to implement, solve the
%    model four times, and combine the resulting state-space, as well as
%    the resulting distributions with proper weights. This is already done
%    in the function SNW_TAX_STEADY_MULTITYPES, so can copy the code over. 
%
%    Several features needed: (1) allow for adjusting parameters, this is
%    straight-forward; (2) allow for considering TRUMP check; (3) same
%    output structure as current, but with expanded state-space. 
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] =
%    SNW_DS_MAIN_VEC_MULTITYPES() invoke model with externally set parameter map and
%    control map. Results outputed to a map containing various output
%    matrixes in mp_dsvfi_results, and also distributional matrixes.
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] =
%    SNW_DS_MAIN_VEC_MULTITYPES(MP_PARAMS, MP_CONTROLS, AP_SS, CONS_SS,
%    MP_VALPOL_MORE_SS, PHI_ADJ_BASE) solves for distribution induces by
%    AP_SS, CONS_SS policy functions and given existing distribution
%    PHI_ADJ_BASE
%
%    See also SNW_DS_MAIN, SNW_DS_GRID_SEARCH, SNWX_DS_BISEC_VEC,
%    SNWX_DS_BISEC_VEC_ONEPERIODPOLSHIFT
%

%%
function varargout=snw_ds_main_vec_multitypes(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))
    
    if (length(varargin)==8)
        [st_param_group_base, st_biden_or_trump, ...
            mp_params_ext, mp_controls_ext, ...
            ls_fl_beta_val, ls_it_edu_simu_type, ...
            fl_wgt_college_frac, fl_wgt_low_beta_college, fl_wgt_low_beta_non_college] = varargin{:};
    else
        error('Need to provide 8 parameter inputs');
    end
    
else
    
    clc;
    clear all;
    
    % 1. Defaults
    % st_param_group_base = 'default_base';
    st_param_group_base = 'default_tiny';
    % st_biden_or_trump = 'bidenchk';
    st_biden_or_trump = 'trumpchk';
    
    % 2a. beta and edu types
    ls_fl_beta_val = [0.60, 0.95];
    
    % 2b. 1 is low education, 2 is high education
    ls_it_edu_simu_type = [1, 2];
    
    % 3. edu and beta probabilities
    fl_wgt_college_frac = 0.303;
    fl_wgt_low_beta_college =  0.1;
    fl_wgt_low_beta_non_college = 0.4;
    
    % 4. mp_params and mp_controls
    mp_params_ext = containers.Map('KeyType','char', 'ValueType','any');
    mp_params_ext('a2_covidyr') = 1.6996;
    mp_controls_ext = containers.Map('KeyType','char', 'ValueType','any');
    
end

%% Solve Model Given Tax
% store results collections
cl_mp_params = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_mp_controls = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_a_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_ap_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_cons_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_y_all_ss = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
cl_Phi_true = cell(length(ls_fl_beta_val), length(ls_it_edu_simu_type));

% Aggregate outcome stores
mt_Y_inc_median = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
mt_a_mean = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
mt_ap_mean = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
mt_cons_mean = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
mt_Y_inc_mean = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));
mt_grp_wgt = zeros(length(ls_fl_beta_val), length(ls_it_edu_simu_type));

% loop and solve
it_beta_ctr = 0;
for fl_beta_val = ls_fl_beta_val
    it_beta_ctr = it_beta_ctr + 1;

    it_edu_ctr = 0;
    for it_edu_simu_type = ls_it_edu_simu_type
        it_edu_ctr = it_edu_ctr + 1;

        % parameter method 1
        [mp_params, mp_controls] = get_param_maps(...
            st_param_group_base, it_edu_simu_type, ...
            fl_beta_val, mp_params_ext, mp_controls_ext);

        % Solve for distributions and policies
        [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(...
            mp_params, mp_controls);
        % solve for distributions
        [Phi_true_ss, ~, A_agg_ss, Y_inc_agg_ss, ~, mp_dsvfi_results_ss] = snw_ds_main_vec(...
            mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
        
        % Trump of Biden Check Related Problems
        if strcmp(st_biden_or_trump, 'trumpchk')
            % no actions needed
            disp('Trump Check, do not need to resolve distribution')

            % Update policies
            Phi_true = Phi_true_ss;
            
            % Aggregate outcomes
            A_agg = A_agg_ss;
            Y_inc_agg = Y_inc_agg_ss;
            mp_dsvfi_results = mp_dsvfi_results_ss;
            
        elseif strcmp(st_biden_or_trump, 'bidenchk')
            % Use steady-state continuation value, to solve for optimal choices
            % given stimulus checks
            % Assume manna-from-heaven, (1) same tax in year-1 of covid as under
            % steady state. (2) assume income loss is fully covered by
            % "unemployment benefits". Because of this, do not have to solve
            % for unemployed and employed separately, their
            % resource-availability are identical. Approximating in effect
            % income in 2019, which is used to determine stimulus checks, with
            % income in 2020, and in effect ignoring one year kids transition
            % probability. Because of these assumptions, the effects of the
            % Trump-stimulus is to strictly increase resources availability for
            % households receiving stimulus, and no impact on those not
            % receiving stimulus. And those that do receive can save more,
            % leading to more savings than under steady-state. 

            disp('Biden Check, resolve for distributions given Trump check')
            xi = mp_params('xi');
            b = mp_params('b');
            mp_params('xi') = 0;
            mp_params('b') = 1;        
            [~,ap_xi0b1_trumpchecks, cons_xi0b1_trumpchecks, mp_valpol_more_trumpchecks] = ...
                snw_vfi_main_bisec_vec_stimulus(mp_params, mp_controls, v_ss);
            mp_params('xi') = xi;
            mp_params('b') = b;

            % Distribution upon arrival in 2nd-covid year, the biden year, update Phi_true_ss        
            [Phi_true_trumpcheck, ~, A_agg_trumpchecks, Y_inc_agg_trumpchecks, ~, mp_dsvfi_results_trumpchecks] = ...
                snw_ds_main_vec(mp_params, mp_controls, ...
                ap_xi0b1_trumpchecks, cons_xi0b1_trumpchecks, ...
                mp_valpol_more_trumpchecks, ...
                Phi_true_ss);
            
            % Update policies
            Phi_true = Phi_true_trumpcheck;

            % Aggregate outcomes
            A_agg = A_agg_trumpchecks;
            Y_inc_agg = Y_inc_agg_trumpchecks;
            mp_dsvfi_results = mp_dsvfi_results_trumpchecks;

        else 
            error(['st_biden_or_trump=' char(st_biden_or_trump) ' is not allowed, has to be trumpchk or bidenchk'])
        end
        
        
        mp_cl_mt_xyz_of_s = mp_dsvfi_results('mp_cl_mt_xyz_of_s');
        tb_outcomes = mp_cl_mt_xyz_of_s('tb_outcomes');
        % Stats stored from snw_ds_main_vec
        Y_inc_median = tb_outcomes{'y_all', 'p50'};

        % Compute the same statistics here. (debug code)
        y_all_ss = mp_dsvfi_results('y_all_ss');
        mn_a_ss = mp_dsvfi_results('a_ss');

        % Store results
        cl_mp_params{it_beta_ctr, it_edu_ctr} = mp_params;
        cl_mp_controls{it_beta_ctr, it_edu_ctr} = mp_controls;
        cl_a_ss{it_beta_ctr, it_edu_ctr} = mn_a_ss;
        cl_ap_ss{it_beta_ctr, it_edu_ctr} = ap_ss;
        cl_cons_ss{it_beta_ctr, it_edu_ctr} = cons_ss;
        cl_y_all_ss{it_beta_ctr, it_edu_ctr} = y_all_ss;
        cl_Phi_true{it_beta_ctr, it_edu_ctr} = Phi_true;

        % Compute the same statistics here. (debug code)
        bl_check_stats = true;
        if (bl_check_stats)
            mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
            mp_cl_ar_xyz_of_s('y_all') = {y_all_ss(:), zeros(1)};            
            mp_cl_ar_xyz_of_s('ap') = {ap_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('mn_a') = {mn_a_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('cons') = {cons_ss(:), zeros(1)};
            mp_cl_ar_xyz_of_s('ar_st_y_name') = ["y_all", "ap", "mn_a", "cons"];
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
            
            Y_inc_median_here = tb_outcomes{'y_all', 'p50'};
            
            % see also: snw_calibrate_tax_beta_normgdp
            Pop = mp_params('Pop');
            Y_inc_mean_here = tb_outcomes{'y_all', 'mean'};
            ap_mean_here = tb_outcomes{'ap', 'mean'};
            a_mean_here = tb_outcomes{'mn_a', 'mean'};
            cons_mean_here = tb_outcomes{'cons', 'mean'};

            if (Y_inc_median_here ~= Y_inc_median)
                error('error Y_inc_median_here ~= Y_inc_median');
            end
            if ((Y_inc_mean_here*sum(Pop)) == Y_inc_agg)
                error('error Y_inc_mean_here == Y_inc_agg, these are different should not be the same');
            end
            
            rel_tol=1e-09;
            abs_tol=0.0;
            if_is_close = @(a,b) (abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol));            
            if (~if_is_close((a_mean_here*sum(Pop)), A_agg))
                error('error Y_inc_mean_here ~= Y_inc_agg');
            end
        end
        
        % Store stats
        mt_Y_inc_median(it_beta_ctr, it_edu_ctr) = Y_inc_median;
        mt_a_mean(it_beta_ctr, it_edu_ctr) = a_mean_here;
        mt_ap_mean(it_beta_ctr, it_edu_ctr) = ap_mean_here;
        mt_cons_mean(it_beta_ctr, it_edu_ctr) = cons_mean_here;
        mt_Y_inc_mean(it_beta_ctr, it_edu_ctr) = Y_inc_mean_here;        
        
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
        mt_grp_wgt(it_beta_ctr, it_edu_ctr) = fl_wgt;
        
    end
end    

% Compute Distributional Statistics for All Groups given Converged policy and tax
% Average Output and Savings
fl_a_mean_alledubeta = sum(mt_a_mean.*mt_grp_wgt, 'all');
fl_ap_mean_alledubeta = sum(mt_ap_mean.*mt_grp_wgt, 'all');
fl_Y_inc_mean_alledubeta = sum(mt_Y_inc_mean.*mt_grp_wgt, 'all');
fl_cons_mean_alledubeta = sum(mt_cons_mean.*mt_grp_wgt, 'all');

% Store full states policy and distributions
[it_d1, it_d2, it_d3, it_d4, it_d5, it_d6]= size(ap_ss);
if (it_d4 ~= 1)
    error(['edu dimension should be 1'])
end

%% Expand Policy And States by One Dimension, weight Beta
% Note beta are rows, edu are columns. 
% first combine betas together for the same edu group
% sufficient for distributional statistics based on the state-space as well
% as the choice space, for weighted average type information, but actually
% not accurate for policy function related percentile, distributional
% non-mean information. Hiding within beta variations. 

% Initialize
a_wgtbeta = NaN([it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
ap_VFI_wgtbeta = NaN([it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
cons_VFI_wgtbeta = NaN([it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
y_all_wgtbeta = NaN([it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
Phi_true_wgtbeta = NaN([it_d1, it_d2, it_d3, 2, it_d5, it_d6]);

% Weighted average over beta, and proportional edu types
it_edu_ctr = 0;
for it_edu_simu_type = ls_it_edu_simu_type
    it_edu_ctr = it_edu_ctr + 1;
    
    % Weighted averaging over beta types
    it_beta_ctr = 0;
    a_beta_weighted = 0;
    ap_beta_weighted = 0;
    cons_beta_weighted = 0;
    y_all_beta_weighted = 0;
    Phi_true_beta_weighted = 0;    
    for fl_beta_val = ls_fl_beta_val
        it_beta_ctr = it_beta_ctr + 1;
            
        fl_wgt = mt_grp_wgt(it_beta_ctr, it_edu_ctr);
        fl_beta_condi_edu = fl_wgt/sum(mt_grp_wgt(:, it_edu_ctr), 'all');
        
        % Steady-state policy functions
        mn_y_all_ss = cl_y_all_ss{it_beta_ctr, it_edu_ctr};
        mn_a_ss = cl_a_ss{it_beta_ctr, it_edu_ctr};
        mn_ap_ss = cl_ap_ss{it_beta_ctr, it_edu_ctr};
        mn_cons_ss = cl_cons_ss{it_beta_ctr, it_edu_ctr};        
        mn_Phi_true = cl_Phi_true{it_beta_ctr, it_edu_ctr};
        % mn_Phi_true_1 = mn_Phi_true/sum(mn_Phi_true, 'all');
        
        y_all_beta_weighted = y_all_beta_weighted + mn_y_all_ss*fl_beta_condi_edu;
        a_beta_weighted = a_beta_weighted + mn_a_ss*fl_beta_condi_edu;
        ap_beta_weighted = ap_beta_weighted + mn_ap_ss*fl_beta_condi_edu;
        cons_beta_weighted = cons_beta_weighted + mn_cons_ss*fl_beta_condi_edu;
        Phi_true_beta_weighted = Phi_true_beta_weighted + mn_Phi_true*fl_beta_condi_edu;        
    end
    
    % Fill edu 1 and 2 as 4th dimension of outputs
    y_all_wgtbeta(:,:,:,it_edu_ctr,:,:) = y_all_beta_weighted;
    a_wgtbeta(:,:,:,it_edu_ctr,:,:) = a_beta_weighted;
    ap_VFI_wgtbeta(:,:,:,it_edu_ctr,:,:) = ap_beta_weighted;
    cons_VFI_wgtbeta(:,:,:,it_edu_ctr,:,:) = cons_beta_weighted;    
    Phi_true_wgtbeta(:,:,:,it_edu_ctr,:,:) = Phi_true_beta_weighted*sum(mt_grp_wgt(:, it_edu_ctr), 'all');
    
end

% Aggregate weighted average, note mn_a is state-space, common across types
% State-space
fl_Y_inc_mean_wgtbeta_innerwgt = sum(y_all_wgtbeta.*(Phi_true_wgtbeta./sum(Phi_true_wgtbeta, 'all')), 'all');
fl_a_mean_wgtbeta_innerwgt = sum(a_wgtbeta.*(Phi_true_wgtbeta./sum(Phi_true_wgtbeta, 'all')), 'all');
% Choice-space
fl_ap_mean_wgtbeta_innerwgt = sum(ap_VFI_wgtbeta.*(Phi_true_wgtbeta./sum(Phi_true_wgtbeta, 'all')), 'all');
fl_cons_mean_wgtbeta_innerwgt = sum(cons_VFI_wgtbeta.*(Phi_true_wgtbeta./sum(Phi_true_wgtbeta, 'all')), 'all');

% Collect results:
mp_dsvfi_results_wgtbeta = containers.Map('KeyType','char', 'ValueType','any');
mp_dsvfi_results_wgtbeta('y_all_wgtbeta') = y_all_wgtbeta;
mp_dsvfi_results_wgtbeta('a_wgtbeta') = a_wgtbeta;
mp_dsvfi_results_wgtbeta('ap_VFI_wgtbeta') = ap_VFI_wgtbeta;
mp_dsvfi_results_wgtbeta('cons_VFI_wgtbeta') = cons_VFI_wgtbeta;
mp_dsvfi_results_wgtbeta('Phi_true_wgtbeta') = Phi_true_wgtbeta;
mp_dsvfi_results_wgtbeta('fl_Y_inc_mean_wgtbeta_innerwgt') = fl_Y_inc_mean_wgtbeta_innerwgt;
mp_dsvfi_results_wgtbeta('fl_a_mean_wgtbeta_innerwgt') = fl_a_mean_wgtbeta_innerwgt;
mp_dsvfi_results_wgtbeta('fl_ap_mean_wgtbeta_innerwgt') = fl_ap_mean_wgtbeta_innerwgt;
mp_dsvfi_results_wgtbeta('fl_cons_mean_wgtbeta_innerwgt') = fl_cons_mean_wgtbeta_innerwgt;

%% Expand Policy and States with Beta Additional Dimension, Do not wgt beta first
% rather than 6, have 7 dimensions, the first dimension is low and high
% beta.

% Initialize
it_beta_n = length(ls_fl_beta_val);
a_betaedu = NaN([it_beta_n, it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
ap_VFI_betaedu = NaN([it_beta_n, it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
cons_VFI_betaedu = NaN([it_beta_n, it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
y_all_betaedu = NaN([it_beta_n, it_d1, it_d2, it_d3, 2, it_d5, it_d6]);
Phi_true_betaedu = NaN([it_beta_n, it_d1, it_d2, it_d3, 2, it_d5, it_d6]);

% Weighted average over beta, and proportional edu types
it_edu_ctr = 0;
for it_edu_simu_type = ls_it_edu_simu_type
    it_edu_ctr = it_edu_ctr + 1;
    
    % Weighted averaging over beta types
    it_beta_ctr = 0;
    for fl_beta_val = ls_fl_beta_val
        it_beta_ctr = it_beta_ctr + 1;           
        
        y_all_betaedu(it_beta_ctr, :,:,:,it_edu_ctr,:,:) = cl_y_all_ss{it_beta_ctr, it_edu_ctr};
        a_betaedu(it_beta_ctr, :,:,:,it_edu_ctr,:,:) = cl_a_ss{it_beta_ctr, it_edu_ctr};
        ap_VFI_betaedu(it_beta_ctr, :,:,:,it_edu_ctr,:,:) = cl_ap_ss{it_beta_ctr, it_edu_ctr};
        cons_VFI_betaedu(it_beta_ctr, :,:,:,it_edu_ctr,:,:) = cl_cons_ss{it_beta_ctr, it_edu_ctr};    
        Phi_true_betaedu(it_beta_ctr, :,:,:,it_edu_ctr,:,:) = cl_Phi_true{it_beta_ctr, it_edu_ctr}*mt_grp_wgt(it_beta_ctr, it_edu_ctr);
    end    
end

% Aggregate weighted average, note mn_a is state-space, common across types
% State-space
fl_Y_inc_mean_betaedu_innerwgt = sum(y_all_betaedu.*(Phi_true_betaedu./sum(Phi_true_betaedu, 'all')), 'all');
fl_a_mean_betaedu_innerwgt = sum(a_betaedu.*(Phi_true_betaedu./sum(Phi_true_betaedu, 'all')), 'all');
% Choice-space
fl_ap_mean_betaedu_innerwgt = sum(ap_VFI_betaedu.*(Phi_true_betaedu./sum(Phi_true_betaedu, 'all')), 'all');
fl_cons_mean_betaedu_innerwgt = sum(cons_VFI_betaedu.*(Phi_true_betaedu./sum(Phi_true_betaedu, 'all')), 'all');

% Collect results:
mp_dsvfi_results_betaedu = containers.Map('KeyType','char', 'ValueType','any');
mp_dsvfi_results_betaedu('y_all_betaedu') = y_all_betaedu;
mp_dsvfi_results_betaedu('a_betaedu') = a_betaedu;
mp_dsvfi_results_betaedu('ap_VFI_betaedu') = ap_VFI_betaedu;
mp_dsvfi_results_betaedu('cons_VFI_betaedu') = cons_VFI_betaedu;
mp_dsvfi_results_betaedu('Phi_true_betaedu') = Phi_true_betaedu;
mp_dsvfi_results_betaedu('fl_Y_inc_mean_betaedu_innerwgt') = fl_Y_inc_mean_betaedu_innerwgt;
mp_dsvfi_results_betaedu('fl_a_mean_betaedu_innerwgt') = fl_a_mean_betaedu_innerwgt;
mp_dsvfi_results_betaedu('fl_ap_mean_betaedu_innerwgt') = fl_ap_mean_betaedu_innerwgt;
mp_dsvfi_results_betaedu('fl_cons_mean_betaedu_innerwgt') = fl_cons_mean_betaedu_innerwgt;

%% Check Equality
bl_check_debug_equality = false;
if bl_check_debug_equality
    % state-space weighting works, but policy function do change, why is that?
    if (fl_Y_inc_mean_wgtbeta_innerwgt ~= fl_Y_inc_mean_alledubeta)
        error(['error fl_Y_inc_mean_wgtbeta_innerwgt ~= fl_Y_inc_mean_alledubeta']);
    end
    if (fl_a_mean_wgtbeta_innerwgt ~= fl_a_mean_alledubeta)
        error(['error fl_a_mean_betaedu_innerwgt ~= fl_a_mean_alledubeta']);
    end
    if (fl_ap_mean_wgtbeta_innerwgt ~= fl_ap_mean_alledubeta)
        warning(['error fl_ap_mean_betaedu_innerwgt ~= fl_ap_mean_alledubeta']);
    end
    if (fl_cons_mean_wgtbeta_innerwgt ~= fl_cons_mean_alledubeta)
        error(['error fl_cons_mean_betaedu_innerwgt ~= fl_cons_mean_alledubeta']);
    end
    
    % state-space weighting works, but policy function do change, why is that?
    if (fl_Y_inc_mean_betaedu_innerwgt ~= fl_Y_inc_mean_alledubeta)
        error(['error fl_Y_inc_mean_betaedu_innerwgt ~= fl_Y_inc_mean_alledubeta']);
    end
    if (fl_a_mean_betaedu_innerwgt ~= fl_a_mean_alledubeta)
        error(['error fl_a_mean_betaedu_innerwgt ~= fl_a_mean_alledubeta']);
    end
    if (fl_ap_mean_betaedu_innerwgt ~= fl_ap_mean_alledubeta)
        warning(['error fl_ap_mean_betaedu_innerwgt ~= fl_ap_mean_alledubeta']);
    end
    if (fl_cons_mean_betaedu_innerwgt ~= fl_cons_mean_alledubeta)
        warning(['error fl_cons_mean_betaedu_innerwgt ~= fl_cons_mean_alledubeta']);
    end
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = Phi_true_betaedu;
    elseif (it_k==2)
        ob_out_cur = ap_VFI_betaedu;
    elseif (it_k==3)
        ob_out_cur = cons_VFI_betaedu;
    elseif (it_k==4)
        ob_out_cur = y_all_betaedu;
    elseif (it_k==5)
        ob_out_cur = mp_dsvfi_results_wgtbeta;
    elseif (it_k==6)
        ob_out_cur = mp_dsvfi_results_betaedu;
    end
    varargout{it_k} = ob_out_cur;
end

end

% Parameter getter
function [mp_params, mp_controls] = get_param_maps(st_param_group_base, it_edu_simu_type, ...
    fl_beta_val, mp_params_ext, mp_controls_ext)
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

% override
mp_params = [mp_params; mp_params_ext];
mp_controls = [mp_controls; mp_controls_ext];

% beta parameter update
mp_params('beta') = fl_beta_val;

end
