%% SNW_VFI_MAIN Solves Policy/Value Function SNW (Bisection Vectorized)
%    Given parameters, iterate over life cycle, given age, marital status,
%    education level and child count, as well as persistent productivity
%    shock process, solve for optimal dynamic savings choices given
%    expectation of kid count transition and productivity shock transition.
%
%    Pref, Technology, and prices SCALARS:
%
%    * BETA discount
%    * THETA total factor productivity normalizer
%    * R interest rate
%
%    Vectorized State Space ARRAYS:
%
%    * AGRID asset grid
%    * ETA_GRID productivity shock grid
%
%    Transition Matrixes ARRAYS:
%
%    * PI_ETA shock productivity transition
%    * PI_KIDS shock kids count transition
%    * PSI shock survival probability
%
%    Permanent Education Type Heterogeneity ARRAYS:
%
%    * EPSILON perfect-foresight education type transition
%    * SS Social Security
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN(MP_PARAMS) invoke
%    model with externally set parameter map MP_PARAMS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN(MP_PARAMS,
%    MP_CONTROLS) invoke model with externally set parameter map MP_PARAMS
%    as well as control mpa MP_CONTROLS.
%
%    See also SNWX_VFI_MAIN, SNW_MP_CONTROL, SNW_MP_PARAM
%

%%
function [V_VFI,ap_VFI,cons_VFI,exitflag_VFI]=snw_vfi_main_bisec_vec(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==1)
        [mp_params] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
    end

else

    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');

end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global beta theta r agrid epsilon eta_grid SS pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in functions that are called by this code
global gamma g_n g_cons a2 cons_allocation_rule jret

%% Parse Model Parameters
params_group = values(mp_params, {'gamma', 'beta', 'theta', 'cons_allocation_rule', ...
    'r', 'g_n', 'g_cons', 'a2', 'jret'});
[gamma, beta, theta, cons_allocation_rule, ...
    r, g_n, g_cons, a2, jret] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_grid'});
[agrid, eta_grid] = params_group{:};

params_group = values(mp_params, {'pi_eta', 'pi_kids', 'psi'});
[pi_eta, pi_kids, psi] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

%% Parse Model Controls
% Minimizer Controls
params_group = values(mp_controls, ...
    {'A_aux', 'B_aux', ...
     'Aeq', 'Beq',...
     'nonlcon', 'options', 'options2'});
[A_aux, B_aux, ...
    Aeq, Beq, ...
    nonlcon, options, options2] = params_group{:};

% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_vfi', 'bl_print_vfi_verbose'});
[bl_print_vfi, bl_print_vfi_verbose] = params_group{:};

% Store Controls
params_group = values(mp_controls, {'bl_vfi_store_all'});
[bl_vfi_store_all] = params_group{:};

%% Define Functions

% Current Function and their Derivatives
if(gamma == 1) 
    f_util = @(c) log(c);
    f_du_da = @(c) -1./(c);
else
    f_util = @(c) (((c).^(1-gamma)-1)./(1-gamma));
    f_du_da = @(c) -1./(c.^gamma);
end

% Utility 
f_U = @(u, Ev) (u + beta.*Ev);
f_FOC = @(duda, devda) (duda + beta.*devda);

%% Timing and Profiling Start
if (bl_timer)
    tic
end

%% Solve optimization problem

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% if (bl_vfi_store_all)
%     y_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     tax_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     SS_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% end

exitflag_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Solve for value function and policy functions by means of backwards induction
for j=n_jgrid:(-1):1 % Age
    
    % A1. Generate the Resources Matrix
    mn_resources = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);    
    mn_z_ctr = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    for a=1:n_agrid % Assets
        it_z_ctr = 0;
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids    
                        % Resources
                        [inc,earn]=individual_income(j,a,eta,educ);
                        spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                        resources = (1+r)*agrid(a) ...
                                    + epsilon(j,educ)*theta*exp(eta_grid(eta)) ...
                                    + SS(j,educ) ...
                                    + (married-1)*spouse_inc ...
                                    - max(0,Tax(inc,(married-1)*spouse_inc));
                        mn_resources(a, eta, educ, married, kids) = resources;
                        % non-asset position
                        it_z_ctr = it_z_ctr + 1;
                        mn_z_ctr(a, eta, educ, married, kids) = it_z_ctr;
                    end
                end
            end
        end
    end
    % array states/shocks all
    ar_resources_amz = mn_resources(:);
    ar_z_ctr_amz = mn_z_ctr(:);
    
    % A2. Solve For EV(ap,z) = EV(ap,zp|z)f(zp|z) for all possible ap points
    % ev = 0 in final decision period.
    mn_ev_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    
    if (j~=n_jgrid)
        
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids
                        for a=1:n_agrid
                            % Add to each cell of mt_ev_ap_z, integrating over f(zp|z)
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    
                                    mn_ev_ap_z(a,eta,educ,married,kids) = mn_ev_ap_z(a,eta,educ,married,kids) ...
                                        + pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)...
                                          *V_VFI(j+1,a,etap,educ,married,kidsp);
                                end
                            end

                        end
                    end
                end
            end
        end

        % B1. z specific EV Slope: EV(ap,z)/d(ap)
        mn_deri_dev_dap = diff(mn_ev_ap_z, 1)./diff(agrid);    
        % B2. ND dimensional Array to 2D dimensional Array:
        mt_ev_ap_z = reshape(mn_ev_ap_z, n_agrid, []);
        mt_deri_dev_dap = reshape(mn_deri_dev_dap, n_agrid-1, []);

        % C1. Generate Resources and other matrixes and vectors
        % C. Generate Vectorized FOC Evaluator
        % x = fl_aprime_frac
        fc_ffi_vec_foc_u_v_ap = @(x) ffi_vec_foc_u_v_ap(...
            x, agrid', ...
            ar_resources_amz, ar_z_ctr_amz, mt_deri_dev_dap, ...
            f_du_da, f_FOC);
        
        % D. Solve via Bisection 
        [ar_opti_saveborr_frac_amz] = ff_optim_bisec_savezrone(fc_ffi_vec_foc_u_v_ap);

        % E. Evaluate at Bounds
        ar_nan_idx = isnan(ar_opti_saveborr_frac_amz);    
        if(sum(ar_nan_idx)>0)
            ar_min_max = [0, 1-1E-5];
            mt_val_min_max = zeros(sum(ar_nan_idx), length(ar_min_max));
            for it_minmax = [1,2]
                [~, mt_val_min_max(:,it_minmax), ~] = ffi_vec_u_v_ap(...
                    ar_min_max(it_minmax), agrid', ...
                    ar_resources_amz(ar_nan_idx), ar_z_ctr_amz(ar_nan_idx), ...
                    mt_ev_ap_z, mt_deri_dev_dap, ...
                    f_util, f_U);
            end
            [~, it_max] =  max(mt_val_min_max, [], 2);
            ar_opti_saveborr_frac_amz(ar_nan_idx) = ar_min_max(it_max);
        end

        % F. Evaluate 
        [ar_aprime_amz, ar_val_opti_amz, ar_c_opti_amz] = ffi_vec_u_v_ap(...
            ar_opti_saveborr_frac_amz, agrid', ...
            ar_resources_amz, ar_z_ctr_amz, ...
            mt_ev_ap_z, mt_deri_dev_dap, ...
            f_util, f_U);

        % G. Record Results    
        mn_val_cur = reshape(ar_val_opti_amz, size(mn_ev_ap_z));
        mn_aprime_cur = reshape(ar_aprime_amz, size(mn_ev_ap_z));
        mn_c_opti_amz = reshape(ar_aprime_amz, size(mn_ev_ap_z));
        
        % H. To main store:
        V_VFI(j,:,:,:,:,:) = mn_val_cur;
        ap_VFI(j,:,:,:,:,:) = mn_aprime_cur;
        cons_VFI(j,:,:,:,:,:) = mn_c_opti_amz;
        
    else

        for a=1:n_agrid % Assets
            for eta=1:n_etagrid % Productivity
                for educ=1:n_educgrid % Educational level
                    for married=1:n_marriedgrid % Marital status
                        for kids=1:n_kidsgrid % Number of kids
                            if j==n_jgrid
                                
                                ap_VFI(j,a,eta,educ,married,kids)=0;
                                cons_VFI(j,a,eta,educ,married,kids) = consumption(j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids));

                                if cons_VFI(j,a,eta,educ,married,kids)<=0
                                    disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                    error('Non-positive consumption')
                                end
                                V_VFI(j,a,eta,educ,married,kids)=utility(cons_VFI(j,a,eta,educ,married,kids),married,kids);
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    if (bl_print_vfi)
        disp(strcat(['SNW_VFI_MAIN: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end

end

%% Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_vfi = strjoin(...
        ["Completed SNW_VFI_MAIN", ...
         ['SNW_MP_PARAM=' mp_params('mp_params_name')], ...
         ['SNW_MP_CONTROL=' mp_controls('mp_params_name')] ...
        ], ";");
    disp(st_complete_vfi);
end

end

% Utility Maximization First Order Conditions
function [ar_dU_dap, ar_aprime] = ...
    ffi_vec_foc_u_v_ap(ar_aprime_frac_amz, ar_a, ...
                       ar_resources_amz, ar_z_ctr_amz, mt_deri_dev_dap,...
                       f_du_da, f_FOC)
% A. Percentage Asset Choice to Level Asset Choices
ar_aprime = ar_aprime_frac_amz.*(ar_resources_amz);

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_ap_near_lower_idx(ar_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% C. Current consumption
ar_c = ar_resources_amz - ar_aprime;

% D. Do not need to check fl_c > 0, because asset bound by 0 to 1 open set
ar_du_dap = f_du_da(ar_c);

% E. the marginal effects of additional asset is determined by the slope
% mt_z_ctr_amz = repmat(ar_z_ctr_amz, [1, size(ar_aprime_frac_amz,2)]);
ar_lin_idx = sub2ind(size(mt_deri_dev_dap), ar_ap_near_lower_idx, ar_z_ctr_amz);
ar_deri_dev_dap = mt_deri_dev_dap(ar_lin_idx);
% ar_deri_dev_dap = reshape(ar_deri_dev_dap, size(mt_z_ctr_amz));

% F. overall first order condition, this is the root search objective
ar_dU_dap = f_FOC(ar_du_dap, ar_deri_dev_dap);

end

% Utility given choices
function [ar_aprime, ar_val, ar_c] = ffi_vec_u_v_ap(...
    ar_aprime_frac, ar_a, ...
    ar_resources_amz, ar_z_ctr_amz, mt_ev_ap_z, mt_deri_dev_dap, ...
    f_util, f_U)
% A. Percentage Asset Choice to Level Asset Choices
ar_aprime = ar_aprime_frac.*(ar_resources_amz);

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_it_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_it_ap_near_lower_idx(ar_it_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% C. Current consumption
ar_c = ar_resources_amz - ar_aprime;

% D. Evaluate Value
ar_u_of_ap = f_util(ar_c);

% the marginal effects of additional asset is determined by the slope
ar_deri_lin_idx = sub2ind(size(mt_deri_dev_dap), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_ev_lin_idx = sub2ind(size(mt_ev_ap_z), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_deri_dev_dap = mt_deri_dev_dap(ar_deri_lin_idx);
ar_ev_ap_lower_idx = mt_ev_ap_z(ar_ev_lin_idx);

% Ev(a_lower_idx,z) + slope*(fl_aprime - fl_a_lower)
ar_ev_aprime_z = ar_ev_ap_lower_idx + (ar_aprime - ar_a(ar_it_ap_near_lower_idx)').*ar_deri_dev_dap;

% overall utility at choice
ar_val = f_U(ar_u_of_ap, ar_ev_aprime_z);
end
