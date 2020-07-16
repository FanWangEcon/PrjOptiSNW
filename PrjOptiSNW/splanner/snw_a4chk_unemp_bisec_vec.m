%% SNW_A4CHK_UNEMP_BISEC_VEC (Exact Vectorized) Asset Position for Check
%    What is the value of a check? From the perspective of the value
%    function? We have Asset as a state variable, in a cash-on-hand sense,
%    how much must the asset (or think cash-on-hand) increase by, so that
%    it is equivalent to providing the household with a check? This is not
%    the same as the check amount because of tax as well as interest rates.
%    Interest rates means that you might need to offer a smaller a than the
%    check amount. The tax rate means that we might need to shift a by
%    larger than the check amount.
%
%    This function is slightly different from SNW_A4CHK_WRK. There we solve
%    the problem for working individuals who do not have an unemployment
%    shock. Here, we have the unemployment shock. This is the slower looped
%    code fzero version.
%
%    This is the faster vectorized solution. It takes as given pre-computed
%    household head and spousal income that is state-space specific, which
%    does not need to be recomputed.
%
%    * WELF_CHECKS integer the number of checks
%    * TR float the value of each check
%    * V_SS ndarray the value matrix along standard state-space dimensions:
%    (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid)
%    * MP_PARAMS map with model parameters
%    * MP_CONTROLS map with control parameters
%
%    [V_W, EXITFLAG_FSOLVE] = SNW_A4CHK_UNEMP(WELF_CHECKS, TR, V_SS,
%    MP_PARAMS, MP_CONTROLS) solves for working value given V_SS value
%    function results, for number of check WELF_CHECKS, and given the value
%    of each check equal to TR.
%
%    See also SNWX_A4CHK_WRK_SMALL, SNWX_A4CHK_WRK_DENSE,
%    SNW_A4CHK_UNEMP_BISEC, SNW_A4CHK_UNEMP_BISEC_VEC, FIND_A_WORKING
%

%%
function [V_U, C_U]=snw_a4chk_unemp_bisec_vec(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==3)
        [welf_checks, V_unemp, cons_unemp] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==5)
        [welf_checks, V_unemp, cons_unemp, mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==8)
        [welf_checks, V_unemp, cons_unemp, mp_params, mp_controls, ...
            ar_a_amz, ar_inc_unemp_amz, ar_spouse_inc_unemp_amz] = varargin{:};        
    else
        error('Need to provide 3/5/8 parameter inputs');
    end
    
else
    close all;
                
    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    
    % Solve for Value Function, Without One Period Unemployment Shock
    [V_ss,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    
    % The number of checks
    welf_checks = 2;
    
    % Solve for Value of One Period Unemployment Shock
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    TR = 100/58056;

    mp_params('TR') = TR;
    mp_params('xi') = xi;
    mp_params('b') = b;
    
    [V_unemp,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);    
    
end

%% Reset All globals
% Parameters used in this code directly
% global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in find_a_working function
% global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r', 'a2', 'jret'});
[theta,  r, a2, jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'TR', 'xi','b'});
[TR, xi, b] = params_group{:};    

%% Parse Model Controls
% Minimizer Controls
params_group = values(mp_controls, {'fl_max_trchk_perc_increase'});
[fl_max_trchk_perc_increase] = params_group{:};
params_group = values(mp_controls, {'options2'});
[options2] = params_group{:};

% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_a4chk','bl_print_a4chk_verbose'});
[bl_print_a4chk, bl_print_a4chk_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tic
end

%% A. Compute Household-Head and Spousal Income

% this is only called when the function is called without mn_inc_plus_spouse_inc
if ~exist('ar_inc_unemp_amz','var')
    
    % initialize
    mn_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_spouse_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_a = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    
    % Txable Income at all state-space points
    for j=1:n_jgrid % Age
        for a=1:n_agrid % Assets
            for eta=1:n_etagrid % Productivity
                for educ=1:n_educgrid % Educational level
                    for married=1:n_marriedgrid % Marital status
                        for kids=1:n_kidsgrid % Number of kids

                            % [inc,earn]=individual_income(j,a,eta,educ,xi,b);
                            % spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                            [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                                theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
                                xi,b);
                            spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                            
                            mn_inc_unemp(j,a,eta,educ,married,kids) = inc;
                            mn_spouse_inc_unemp(j,a,eta,educ,married,kids) = (married-1)*spouse_inc*exp(eta_S_grid(eta));
                            mn_a(j,a,eta,educ,married,kids) = agrid(a);
                            
                        end
                    end
                end
            end
        end
    end
    
    % flatten the nd dimensional array
    ar_inc_unemp_amz = mn_inc_unemp(:);
    ar_spouse_inc_unemp_amz = mn_spouse_inc_unemp(:);
    ar_a_amz = mn_a(:);
    
end

%% B. Vectorized Solution for Optimal Check Adjustments

% B1. Anonymous Function where X is fraction of addition given bounds
fc_ffi_frac0t1_find_a_working = @(x) ffi_frac0t1_find_a_unemp_vec(...
    x, ...
    ar_a_amz, ar_inc_unemp_amz, ar_spouse_inc_unemp_amz, ...
    welf_checks, TR, r, a2, fl_max_trchk_perc_increase);

% B2. Solve with Bisection
[~, ar_a_aux_unemp_bisec_amz] = ...
    ff_optim_bisec_savezrone(fc_ffi_frac0t1_find_a_working);

% B3. Reshape
mn_a_aux_unemp_bisec = reshape(ar_a_aux_unemp_bisec_amz, [n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid]);

%% C. Interpolate and Extrapolate All Age Check Values
% mn = nd dimensional array, mt = 2d
% ad1 = a is dimension 1

% C1. Change matrix order so asset becomes the first dimension
mn_v_ad1 = permute(V_unemp, [2,1,3,4,5,6]);
mn_c_ad1 = permute(cons_unemp, [2,1,3,4,5,6]);
mn_a_aux_bisec_ad1 = permute(mn_a_aux_unemp_bisec, [2,1,3,4,5,6]);

% C2. Reshape so that asset is dim 1, all other dim 2
mt_v_ad1 = reshape(mn_v_ad1, n_agrid, []);
mt_c_ad1 = reshape(mn_c_ad1, n_agrid, []);
mt_a_aux_bisec_ad1 = reshape(mn_a_aux_bisec_ad1, n_agrid, []);

% C3. Derivative dv/da
mt_dv_da_ad1 = diff(mt_v_ad1, 1)./diff(agrid);    
mt_dc_da_ad1 = diff(mt_c_ad1, 1)./diff(agrid);

% C4. Optimal aux and closest a index
ar_a_aux_bisec = mt_a_aux_bisec_ad1(:);
ar_it_a_near_lower_idx = sum(agrid' <= ar_a_aux_bisec, 2);
ar_it_a_near_lower_idx(ar_it_a_near_lower_idx == length(agrid)) = length(agrid) - 1;

% C5. Z index
mt_z_ctr = repmat((1:size(mt_v_ad1, 2)), size(mt_v_ad1, 1), 1);
ar_z_ctr = mt_z_ctr(:);

% C6. the marginal effects of additional asset is determined by the slope
ar_deri_lin_idx = sub2ind(size(mt_dv_da_ad1), ar_it_a_near_lower_idx, ar_z_ctr);
ar_v_lin_idx = sub2ind(size(mt_v_ad1), ar_it_a_near_lower_idx, ar_z_ctr);
ar_deri_dev_dap = mt_dv_da_ad1(ar_deri_lin_idx);
ar_v_a_lower_idx = mt_v_ad1(ar_v_lin_idx);
ar_deri_dc_da = mt_dc_da_ad1(ar_deri_lin_idx);
ar_c_a_lower_idx = mt_c_ad1(ar_v_lin_idx);

% C7. v(a_lower_idx,z) + slope*(fl_a_aux - fl_a)
ar_v_u_a_aux_j = ar_v_a_lower_idx + (ar_a_aux_bisec - agrid(ar_it_a_near_lower_idx)).*ar_deri_dev_dap;
ar_c_u_a_aux_j = ar_c_a_lower_idx + (ar_a_aux_bisec - agrid(ar_it_a_near_lower_idx)).*ar_deri_dc_da;
mn_v_u_a_aux_ad1 = reshape(ar_v_u_a_aux_j, n_agrid, n_jgrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid);
mn_c_u_a_aux_ad1 = reshape(ar_c_u_a_aux_j, n_agrid, n_jgrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid);

% C8. Permute age as dim 1
V_U = permute(mn_v_u_a_aux_ad1, [2,1,3,4,5,6]);
C_U = permute(mn_c_u_a_aux_ad1, [2,1,3,4,5,6]);

%% D. Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_a4chk = strjoin(...
        ["Completed SNW_A4CHK_UNEMP_BISEC_VEC", ...
         ['welf_checks=' num2str(welf_checks)], ...
         ['TR=' num2str(TR)], ...
         ['xi=' num2str(xi)], ...
         ['b=' num2str(b)], ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_a4chk);
end


%% Compare Difference between V_ss and V_W

if (bl_print_a4chk_verbose)
    
    mn_V_gain_check = V_U - V_unemp;
    mn_C_gain_check = C_U - cons_unemp;
    mn_MPC = (C_U - cons_unemp)./(welf_checks*TR);    
    mp_container_map = containers.Map('KeyType','char', 'ValueType','any');
    mp_container_map('V_U') = V_U;
    mp_container_map('C_U') = C_U;
    mp_container_map('V_U_minus_V_unemp') = mn_V_gain_check;
    mp_container_map('C_U_minus_C_unemp') = mn_C_gain_check;    
    mp_container_map('mn_MPC_unemp') = mn_MPC;    
    ff_container_map_display(mp_container_map);
    
end

end

% this function is in fact identical to ffi_frac0t1_find_a_working_vec from
% snw_a4chk_wrk_bisec_vec.m
function [ar_root_zero, ar_a_aux_amz] = ...
    ffi_frac0t1_find_a_unemp_vec(...
    ar_aux_change_frac_amz, ...
    ar_a_amz, ar_inc_unemp_amz, ar_spouse_unemp_amz, ...
    welf_checks, TR, r, a2, fl_max_trchk_perc_increase)
    
    % Max A change to account for check
    fl_a_aux_max = TR*welf_checks*fl_max_trchk_perc_increase;    
    
    % Level of A change
    ar_a_aux_amz = ar_a_amz + ar_aux_change_frac_amz.*fl_a_aux_max;
    
    % Account for Interest Rates
    ar_r_gap = (1+r).*(ar_a_amz - ar_a_aux_amz);
    
    % Account for tax, inc changes by r
    ar_tax_gap = ...
          max(0, snw_tax_hh(ar_inc_unemp_amz, ar_spouse_unemp_amz, a2)) ...
        - max(0, snw_tax_hh(ar_inc_unemp_amz - ar_a_amz*r + ar_a_aux_amz*r, ar_spouse_unemp_amz, a2));
    
    % difference equation f(a4chkchange)=0
    ar_root_zero = TR*welf_checks + ar_r_gap - ar_tax_gap;
end