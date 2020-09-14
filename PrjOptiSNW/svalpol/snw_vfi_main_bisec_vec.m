%% SNW_VFI_MAIN_BISEC_VEC Solves Policy/Value Function SNW (Bisection Vectorized)
%    Given parameters, iterate over life cycle, given age, marital status,
%    education level and child count, as well as persistent productivity
%    shock process, solve for optimal dynamic savings choices given
%    expectation of kid count transition and productivity shock transition.
%
%    Note that the solution algorithm is exact given the discrete asset
%    grid, but not exact with respect to the underlying continuous
%    state/choice problem we are trying to approximate.
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

%%
function [varargout]=snw_vfi_main_bisec_vec(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    bl_covid_year = false;
    if (length(varargin)==1)
        [mp_params] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==3)
        [mp_params, mp_controls, V_VFI_POSTSHOCK] = varargin{:};
    end
    
    if (length(varargin)==3)
        bl_covid_year = true;
    end
    
else
    
    clc;
    bl_covid_year = false;
    mp_params = snw_mp_param('default_tiny', false, 'tauchen', false, 8, 8);
    mp_controls = snw_mp_control('default_test');
    false;

end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
% global beta theta r agrid epsilon SS pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in functions that are called by this code
% global gamma g_n g_cons a2 cons_allocation_rule jret
% July 1st new parameters
% global eta_H_grid eta_S_grid Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {...
    'gamma', 'beta', 'theta', 'cons_allocation_rule', ...
    'r', 'g_n', 'g_cons', 'a2', 'a2_covidyr', 'jret'});
[gamma, beta, theta, cons_allocation_rule, ...
    r, g_n, g_cons, a2, a2_covidyr, jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'bl_store_shock_trans', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, bl_store_shock_trans, cl_mt_pi_jem_kidseta, psi] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

% unemployment parameters under covid
if (length(varargin)==3)
    params_group = values(mp_params, {'xi','b'});
    [xi, b] = params_group{:};
end

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_vfi', 'bl_print_vfi_verbose'});
[bl_print_vfi, bl_print_vfi_verbose] = params_group{:};

% Store Controls
bl_vfi_store_all = false;
if (nargout >= 4)
    bl_vfi_store_all = true;
end

%% Define Functions

% Current Function and their Derivatives
hh_power=1/cons_allocation_rule;
if(gamma == 1)
    f_util = @(c,hh_size) log(c./(hh_size.^hh_power));
    f_du_da = @(c,hh_size) (-1)./(c);
else
    f_util = @(c,hh_size) ((c./(hh_size.^hh_power)).^(1-gamma)-1)./(1-gamma);
    f_du_da = @(c,hh_size) (-(hh_size.^(hh_power.*(gamma-1))))./(c.^gamma);
end

% Utility
f_U = @(u, Ev) (u + beta.*Ev);
f_FOC = @(duda, devda) (duda + beta.*devda);

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Solve optimization problem
% 1. resources (savings choices are shares of these)
% 2. Ev(ap,z)
% 3. derivatives: du/dap, dev/dap (linear spline)
% 4. Exact solution given derivatives
% 5. If resource max ap > agrid bound, extrapolate based on final spline

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
if (bl_vfi_store_all)
    inc_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    earn_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    spouse_inc_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    SS_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    tax_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
end

if (length(varargin)==3)
    ar_j_seq = 1:n_jgrid;
else
    ar_j_seq = n_jgrid:(-1):1; % Age
end

if (bl_covid_year)
    a2 = a2_covidyr;
    if (isnan(a2))
        error('a2 can not be NaN');
    end
end

% Solve for value function and policy functions by means of backwards induction
for j=ar_j_seq % Age
    
    % Age Timer
    if (bl_print_vfi) tm_vfi_age = tic; end
    
    % A1. Generate the Resources Matrix
    mn_resources = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_z_ctr = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_hh_size = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

    for a=1:n_agrid % Assets
        it_z_ctr = 0;
        % order here is reverse of order of loops at L181, so that it_z_ctr
        % match up with column order after reshape(mn_ev_ap_z, n_agrid, [])
        for kids=1:n_kidsgrid % Number of kids
            for married=1:n_marriedgrid % Marital status
                for educ=1:n_educgrid % Educational level
                    for eta=1:n_etagrid % Productivity

                        % Resources
                        if (length(varargin)>=3)
                            % one period unemployed shock
                            [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                                theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
                                xi,b);
                            % do not earn one hundred percent
                            fl_earn_ratio = (xi+b*(1-xi));
                        else
                            [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                                theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                            fl_earn_ratio = 1;
                        end

                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        resources = (1+r)*(agrid(a)+Bequests*(bequests_option-1)) ...
                                    + epsilon(j,educ)*theta*exp(eta_H_grid(eta))*fl_earn_ratio ...
                                    + SS(j,educ) ...
                                    + (married-1)*spouse_inc*exp(eta_S_grid(eta))*fl_earn_ratio ...
                                    - max(0,snw_tax_hh(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*fl_earn_ratio,a2));
                        mn_resources(a, eta, educ, married, kids) = resources;

                        % Non-asset position Counter
                        it_z_ctr = it_z_ctr + 1;
                        mn_z_ctr(a, eta, educ, married, kids) = it_z_ctr;

                        % Household Size: 1 kid married, hh_size = 2 + 2 -
                        % 1 = 3. 0 kid and unmarried hh_size = 1 + 1 - 1 =
                        % 1
                        hh_size = married + kids - 1; % m=1 if single; m=2 if married; k=1 if 0 children
                        mn_hh_size(a, eta, educ, married, kids) = hh_size;

                        % Store More Information for Analysis
                        if (bl_vfi_store_all)
                            inc_VFI(j,a,eta,educ,married,kids) = inc;
                            earn_VFI(j,a,eta,educ,married,kids) = earn*fl_earn_ratio;
                            spouse_inc_VFI(j,a,eta,educ,married,kids) = (married-1)*spouse_inc*exp(eta_S_grid(eta));
                            SS_VFI(j,a,eta,educ,married,kids) = SS(j,educ);
                            tax_VFI(j,a,eta,educ,married,kids) = max(0,snw_tax_hh(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)),a2));
                        end

                    end
                end
            end
        end
    end

    % array states/shocks all
    ar_resources_amz = mn_resources(:);
    ar_z_ctr_amz = mn_z_ctr(:);
    ar_hh_size_amz = mn_hh_size(:);

    % A2. Solve For EV(ap,z) = EV(ap,zp|z)f(zp|z) for all possible ap points
    % ev = 0 in final decision period.
    mn_ev_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

    if (j~=n_jgrid)
           
        if (length(varargin)==3)
            V_VFI_FUTURE = V_VFI_POSTSHOCK;
        else
            V_VFI_FUTURE = V_VFI;
        end
 
        for educ=1:n_educgrid % Educational level
            for married=1:n_marriedgrid % Marital status

                % A2a. Get P(S'|S), S = [eta x kids] by [eta x kids] transition matrix
                % do not pre-compute save full transition of shocks, too
                % much memory, requires 100s GB for full model
                if (bl_store_shock_trans)
                    mt_pi_jem_kidseta = cl_mt_pi_jem_kidseta{j,educ,married};
                else
                    mt_pi_jem_kidseta = kron(pi_kids(:,:,j,educ,married), pi_eta);
                end

                % A2b. Get age/edu/marry submatrix
                mn_ev20_jem = permute(V_VFI_FUTURE(j+1,:,:,educ,married,:), [2,3,6,1,4,5]);

                % A2c. 2D rows are savings states, columns are [eta x kids]
                mt_ev20_jem = reshape(mn_ev20_jem, n_agrid, []);

                % A2d. EV = V([a] by [eta x kids]) x TRANS([eta x kids] by [eta x kids])
                mt_ev19_jem = psi(j)*mt_ev20_jem*mt_pi_jem_kidseta';

                % A2e. Reshape Back and Store
                mn_ev19_jem = reshape(mt_ev19_jem, [n_agrid, n_etagrid, 1, 1, n_kidsgrid]);
                mn_ev_ap_z(:, :, educ, married, :) = mn_ev19_jem;
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
            ar_resources_amz, ar_z_ctr_amz, ar_hh_size_amz, ...
            mt_deri_dev_dap, ...
            f_du_da, f_FOC);

        % D. Solve via Bisection
        mp_bisec_ctrlinfo = containers.Map('KeyType','char', 'ValueType','any');
        mp_bisec_ctrlinfo('it_bisect_max_iter') = 30;
        mp_bisec_ctrlinfo('fl_x_left_start') = 10e-6;
        mp_bisec_ctrlinfo('fl_x_right_start') = 1-10e-6;

        [ar_opti_saveborr_frac_amz] = ...
            ff_optim_bisec_savezrone(fc_ffi_vec_foc_u_v_ap,...
            false, false, mp_bisec_ctrlinfo);

        % E. Evaluate at Bounds
        ar_nan_idx = isnan(ar_opti_saveborr_frac_amz);
        if(sum(ar_nan_idx)>0)
            ar_min_max = [0, 1-1E-5];
            mt_val_min_max = zeros(sum(ar_nan_idx), length(ar_min_max));
            for it_minmax = [1,2]
                [~, mt_val_min_max(:,it_minmax), ~] = ffi_vec_u_v_ap(...
                    ar_min_max(it_minmax), agrid', ...
                    ar_resources_amz(ar_nan_idx), ar_z_ctr_amz(ar_nan_idx), ar_hh_size_amz(ar_nan_idx), ...
                    mt_ev_ap_z, mt_deri_dev_dap, ...
                    f_util, f_U);
            end
            [~, it_max] =  max(mt_val_min_max, [], 2);
            ar_opti_saveborr_frac_amz(ar_nan_idx) = ar_min_max(it_max);
        end

        % F. Evaluate
        [ar_aprime_amz, ar_val_opti_amz, ar_c_opti_amz] = ffi_vec_u_v_ap(...
            ar_opti_saveborr_frac_amz, agrid', ...
            ar_resources_amz, ar_z_ctr_amz, ar_hh_size_amz, ...
            mt_ev_ap_z, mt_deri_dev_dap, ...
            f_util, f_U);

        % G. Record Results
        mn_val_cur = reshape(ar_val_opti_amz, size(mn_ev_ap_z));
        mn_aprime_cur = reshape(ar_aprime_amz, size(mn_ev_ap_z));
        mn_c_opti_amz = reshape(ar_c_opti_amz, size(mn_ev_ap_z));

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
                                                                
                                if (length(varargin)==3)
                                    cons_VFI(j,a,eta,educ,married,kids) = snw_hh_consumption(...
                                        j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids),...
                                        theta, r, agrid, epsilon, eta_H_grid, eta_S_grid, SS, Bequests, bequests_option,...
                                        jret,a2,...
                                        xi,b);
                                else
                                    cons_VFI(j,a,eta,educ,married,kids) = snw_hh_consumption(...
                                        j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids),...
                                        theta, r, agrid, epsilon, eta_H_grid, eta_S_grid, SS, Bequests, bequests_option,...
                                        jret,a2);                                    
                                end


                                if cons_VFI(j,a,eta,educ,married,kids)<=0
                                    disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                    error('Non-positive consumption')
                                end

                                V_VFI(j,a,eta,educ,married,kids)=snw_hh_utility(...
                                    cons_VFI(j,a,eta,educ,married,kids),married,kids,...
                                    gamma, cons_allocation_rule);

                            end
                        end
                    end
                end
            end
        end

    end

    if (bl_print_vfi)
        tm_vfi_age_end = toc(tm_vfi_age);
        if (length(varargin)==3)
            disp(strcat(['SNW_VFI_MAIN_BISEC_VEC 1 Period Unemp Shock: Age ' ...
                num2str(j) ' of ' num2str(n_jgrid-1) ...
                ', time-this-age:' num2str(tm_vfi_age_end)]));        
        else
            disp(strcat(['SNW_VFI_MAIN_BISEC_VEC: Finished Age Group:' ...
                num2str(j) ' of ' num2str(n_jgrid-1) ...
                ', time-this-age:' num2str(tm_vfi_age_end)]));        
        end
    end

end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    if (length(varargin)==3)
        st_complete_vfi = strjoin(...
            ["Completed SNW_VFI_MAIN_BISEC_VEC 1 Period Unemp Shock", ...
             ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
             ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
             ['time=' num2str(tm_end)] ...
            ], ";");
    else
        st_complete_vfi = strjoin(...
            ["Completed SNW_VFI_MAIN_BISEC_VEC", ...
             ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
             ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
             ['time=' num2str(tm_end)] ...
            ], ";");
    end
    disp(st_complete_vfi);
end

%% Print
if (bl_print_vfi_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('V_VFI') = V_VFI;
    mp_outcomes('ap_VFI') = ap_VFI;
    mp_outcomes('cons_VFI') = cons_VFI;
    ff_container_map_display(mp_outcomes, 9, 9);
end

%% Return 
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_VFI;
    elseif (it_k==2)
        ob_out_cur = ap_VFI;
    elseif (it_k==3)
        ob_out_cur = cons_VFI;
    elseif (it_k==4)
        mp_valpol_more = containers.Map('KeyType','char', 'ValueType','any');        
        mp_valpol_more('inc_VFI') = inc_VFI;
        mp_valpol_more('earn_VFI') = earn_VFI;
        mp_valpol_more('spouse_inc_VFI') = spouse_inc_VFI;
        mp_valpol_more('SS_VFI') = SS_VFI;
        mp_valpol_more('tax_VFI') = tax_VFI;
        ob_out_cur = mp_valpol_more;
    end
    varargout{it_k} = ob_out_cur;
end

end

% Utility Maximization First Order Conditions
function [ar_dU_dap, ar_aprime] = ...
    ffi_vec_foc_u_v_ap(ar_aprime_frac_amz, ar_a, ...
                       ar_resources_amz, ar_z_ctr_amz, ar_hh_size_amz, ...
                       mt_deri_dev_dap,...
                       f_du_da, f_FOC)
% A. Percentage Asset Choice to Level Asset Choices
ar_aprime = ar_aprime_frac_amz.*(ar_resources_amz);

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_ap_near_lower_idx(ar_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% C. Current consumption
ar_c = ar_resources_amz - ar_aprime;

% D. Do not need to check fl_c > 0, because asset bound by 0 to 1 open set
ar_du_dap = f_du_da(ar_c, ar_hh_size_amz);

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
    ar_resources_amz, ar_z_ctr_amz, ar_hh_size_amz, ...
    mt_ev_ap_z, mt_deri_dev_dap, ...
    f_util, f_U)
% A. Percentage Asset Choice to Level Asset Choices
ar_aprime = ar_aprime_frac.*(ar_resources_amz);

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_it_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_it_ap_near_lower_idx(ar_it_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% C. Current consumption
ar_c = ar_resources_amz - ar_aprime;

% D. Evaluate Value
ar_u_of_ap = f_util(ar_c, ar_hh_size_amz);

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
