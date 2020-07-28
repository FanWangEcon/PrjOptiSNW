%% SNW_EVUVW19_JAEEMKF_FOC Solves for EVUVW from 2019 Conditional JAEEMK
%    Given 2020 JAEEMK, what is the expectated value for the planner given
%    2020 JAEEMK and transition between 2019 to 2020 JAEEMK. For one Check,
%    the check number set by WELF_CHECKS.
%
%    [EV19_JAEEMK, EC19_JAEEMK, EV20_JAEEMK, EC20_JAEEMK] =
%    SNW_EVUVW19_JAEEMK_FOC(WELF_CHECKS, ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS,
%    V_SS, CONS_SS, V_UNEMP, CONS_UNEMP, MP_PRECOMPUTE_RES) provide V_SS
%    and V_UNEMP solved out elsewhere, and only get two outputs out.
%
%    See also SNW_EVUVW19_JYMK, SNW_EVUVW19_JAEEMK, SNW_EVUVW20_JAEEMK,
%    SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jaeemk_foc(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==10)
        [welf_checks, st_solu_type, mp_params, mp_controls, ...
            V_ss, ap_ss, cons_ss, ...
            V_unemp, cons_unemp, ...
            mp_precompute_res] = varargin{:};
    else
        error('Need to provide 10 parameter inputs');
    end

else
    clc;
    close all;

%     st_solu_type = 'matlab_minimizer';
    st_solu_type = 'bisec_vec';
%     st_solu_type = 'grid_search';

    % Solve the VFI Problem and get Value Function
%     mp_params = snw_mp_param('default_tiny');
%     mp_params = snw_mp_param('default_dense');
    mp_params = snw_mp_param('default_moredense');
    mp_controls = snw_mp_control('default_test');

    % The Number of Checks to Provide
    welf_checks = 2;

    % set Unemployment Related Variables
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values

    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_params('TR') = TR;

    % Solve for Unemployment Values
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_precompute') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_a4chk_verbose') = false;

    % Solve the Model to get V working and unemployed
    [V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [V_unemp,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);

    % Get Matrixes
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid',...
        'ar_z_ctr_amz'};
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);

end

%% Parse Pre-Computes
params_group = values(mp_precompute_res, ...
    {'ar_z_ctr_amz'});
[ar_z_ctr_amz] = params_group{:};

%% Parse Model Parameters
params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, cl_mt_pi_jem_kidseta, psi] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_evuvw19_jaeemk', 'bl_print_evuvw19_jaeemk_verbose'});
[bl_print_evuvw19_jaeemk, bl_print_evuvw19_jaeemk_verbose] = params_group{:};

%% Solve evuvw19 Given Current Check
[ev20_jaeemk, ec20_jaeemk] = snw_evuvw20_jaeemk(...
    welf_checks, ...
    st_solu_type, mp_params, mp_controls, ...
    V_ss, cons_ss, ...
    V_unemp, cons_unemp, ...
    mp_precompute_res);

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Planner val2019(J,A,E,E,M,K) = E(val2020(J+1,A',E',E,M,K')) given trans.
ev19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ec19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid-1 % Age
    
    % A1. mn_z_ctr and mn_aprime
    % array states/shocks all
    ar_aprime_amz = reshape(ap_ss(j,:,:,:,:,:), [], 1);
    
    % B1. Solve For EV(ap,z) = EV(ap,zp|z)f(zp|z) for all possible ap points
    mn_ev_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_ec_ap_z = zeros(n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    for educ=1:n_educgrid % Educational level
        for married=1:n_marriedgrid % Marital status
            
            % B1a. Get P(S'|S), S = [eta x kids] by [eta x kids] transition matrix
            mt_pi_jem_kidseta = cl_mt_pi_jem_kidseta{j,educ,married};
            
            % B1b. Get age/edu/marry submatrix
            mn_ev20_jem = permute(ev20_jaeemk(j+1,:,:,educ,married,:), [2,3,6,1,4,5]);
            mn_ec20_jem = permute(ec20_jaeemk(j+1,:,:,educ,married,:), [2,3,6,1,4,5]);
            
            % B1c. 2D rows are savings states, columns are [eta x kids]
            mt_ev20_jem = reshape(mn_ev20_jem, n_agrid, []);
            mt_ec20_jem = reshape(mn_ec20_jem, n_agrid, []);
            
            % B1d. EV = V([a] by [eta x kids]) x TRANS([eta x kids] by [eta x kids])
            mt_ev19_jem = mt_ev20_jem*mt_pi_jem_kidseta';
            mt_ec19_jem = mt_ec20_jem*mt_pi_jem_kidseta';
            
            % B1e. Reshape Back and Store
            mn_ev19_jem = reshape(mt_ev19_jem, [n_agrid, n_etagrid, 1, 1, n_kidsgrid]);
            mn_ec19_jem = reshape(mt_ec19_jem, [n_agrid, n_etagrid, 1, 1, n_kidsgrid]);
            mn_ev_ap_z(:, :, educ, married, :) = mn_ev19_jem;
            mn_ec_ap_z(:, :, educ, married, :) = mn_ec19_jem;
            
        end
    end
    
    % C1. z specific EV Slope: EV(ap,z)/d(ap)
    mn_deri_dev_dap = diff(mn_ev_ap_z, 1)./diff(agrid);
    mn_deri_dec_dap = diff(mn_ec_ap_z, 1)./diff(agrid);
    % C2. ND dimensional Array to 2D dimensional Array:
    mt_ev_ap_z = reshape(mn_ev_ap_z, n_agrid, []);
    mt_ec_ap_z = reshape(mn_ec_ap_z, n_agrid, []);
    mt_deri_dev_dap = reshape(mn_deri_dev_dap, n_agrid-1, []);
    mt_deri_dec_dap = reshape(mn_deri_dec_dap, n_agrid-1, []);
    
    % D, Evaluate
    [ar_ev_aprime_z] = ffi_vec_v_ap(...
        ar_aprime_amz, agrid', ...
        ar_z_ctr_amz, ...
        mt_ev_ap_z, mt_deri_dev_dap);
    [ar_ec_aprime_z] = ffi_vec_v_ap(...
        ar_aprime_amz, agrid', ...
        ar_z_ctr_amz, ...
        mt_ec_ap_z, mt_deri_dec_dap);

    % G. Record Results
    mn_ev_aprime_z = reshape(ar_ev_aprime_z, size(mn_ev_ap_z));
    mn_ec_aprime_z = reshape(ar_ec_aprime_z, size(mn_ev_ap_z));

    % H. To main store:
    ev19_jaeemk(j,:,:,:,:,:) = mn_ev_aprime_z;
    ec19_jaeemk(j,:,:,:,:,:) = mn_ec_aprime_z;

end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_EVUVW19_JAEEMK", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_evuvw19_jaeemk_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('ev19_jaeemk') = ev19_jaeemk;
    mp_outcomes('ec19_jaeemk') = ec19_jaeemk;
    mp_outcomes('ev20_jaeemk') = ev20_jaeemk;
    mp_outcomes('ec20_jaeemk') = ec20_jaeemk;
    ff_container_map_display(mp_outcomes, 9, 9);    
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev19_jaeemk;
    elseif (it_k==2)
        ob_out_cur = ec19_jaeemk;
    elseif (it_k==3)
        ob_out_cur = ev20_jaeemk;
    elseif (it_k==4)
        ob_out_cur = ec20_jaeemk;
    end
    varargout{it_k} = ob_out_cur;
end

end

% Utility given choices
function [ar_ev_aprime_z] = ffi_vec_v_ap(...
    ar_aprime, ar_a, ar_z_ctr_amz, ...
    mt_ev_ap_z, mt_deri_dev_dap)

% B. Identify the Closest ar_a point to fl_aprime, this is spline knot point
ar_it_ap_near_lower_idx = sum(ar_a <= ar_aprime, 2);
ar_it_ap_near_lower_idx(ar_it_ap_near_lower_idx == length(ar_a)) = length(ar_a) - 1;

% the marginal effects of additional asset is determined by the slope
ar_deri_lin_idx = sub2ind(size(mt_deri_dev_dap), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_ev_lin_idx = sub2ind(size(mt_ev_ap_z), ar_it_ap_near_lower_idx, ar_z_ctr_amz);
ar_deri_dev_dap = mt_deri_dev_dap(ar_deri_lin_idx);
ar_ev_ap_lower_idx = mt_ev_ap_z(ar_ev_lin_idx);

% Ev(a_lower_idx,z) + slope*(fl_aprime - fl_a_lower)
ar_ev_aprime_z = ar_ev_ap_lower_idx + (ar_aprime - ar_a(ar_it_ap_near_lower_idx)').*ar_deri_dev_dap;
end
