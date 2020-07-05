%% SNW_A4CHK_WRK_BISEC_VEC (Exact Vectorized) Asset Position Corresponding to Check Level
%    What is the value of a check? From the perspective of the value
%    function? We have Asset as a state variable, in a cash-on-hand sense,
%    how much must the asset (or think cash-on-hand) increase by, so that
%    it is equivalent to providing the household with a check? This is not
%    the same as the check amount because of tax as well as interest rates.
%    Interest rates means that you might need to offer a smaller a than the
%    check amount. The tax rate means that we might need to shift a by
%    larger than the check amount.
%
%    This is the faster vectorized solution. It takes as given pre-computed
%    household head and spousal income that is state-space specific, which
%    does not need to be recomputed.
%
%    * WELF_CHECKS integer the number of checks
%    * TR float the value of each check
%    * V_SS ndarray the value matrix along standard state-space dimensions:
%    (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid)
%    * MN_IT_MESH_CTR ND array of meshed counters, this can be produced
%    once outside so that it does not need to be regenerated each time. It
%    is the ND mesh of all the looping counters, 1:ctr as array for eahc.
%    This has six dimensions, meshed together following:
%    j,a,eta,educ,married,kids dimensions
%    * MP_PARAMS map with model parameters
%    * MP_CONTROLS map with control parameters
%
%    % bl_fzero is true if solve via fzero
%    mp_controls('bl_fzero') = false;
%    $ bl_ff_bisec is true is solve via ff_optim_bisec_savezrone, if both
%    true solve both and return error message if answers do not match up.
%    mp_controls('bl_ff_bisec') = true;
%
%    [V_W, EXITFLAG_FSOLVE] = SNW_A4CHK_WRK_BISEC(WELF_CHECKS, TR, V_SS,
%    MP_PARAMS, MP_CONTROLS) solves for working value given V_SS value
%    function results, for number of check WELF_CHECKS, and given the value
%    of each check equal to TR.
%
%    [V_W, EXITFLAG_FSOLVE] = SNW_A4CHK_WRK_BISEC(WELF_CHECKS, TR, V_SS,
%    MP_PARAMS, MP_CONTROLS, AR_A_AMZ, AR_INC_AMZ, AR_SPOUSE_INC_AMZ)
%    AR_A_AMZ is the flattened nd dimensional array that stores the asset
%    state value for each element of the state space AR_INC_AMZ is the
%    flattend nd dimensional array of all household income head together.
%    AR_SPOUSE_INC_AMZ is the flattend nd dimensional array of spousal
%    income. Kept separately in case they are not linearly additive for the
%    tax function.
%
%    See also SNWX_A4CHK_WRK_BISEC_VEC_DENSE,
%    SNWX_A4CHK_WRK_BISEC_VEC_SMALL, SNW_A4CHK_WRK, SNW_A4CHK_WRK_BISEC
%

%%
function [V_W]=snw_a4chk_wrk_bisec_vec(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==2)
        [welf_checks, V_ss] = varargin{:};
    elseif (length(varargin)==4)
        [welf_checks, V_ss, mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==7)
        [welf_checks, V_ss, mp_params, mp_controls, ...
            ar_a_amz, ar_inc_amz, ar_spouse_inc_amz] = varargin{:};
    else
        error('Need to provide 2/4/7 parameter inputs');
    end
    
else
    close all;
    
    % A1. Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    [V_ss,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    
    % Solve for Value of One Period Unemployment Shock
    welf_checks = 2;    
    TR = 100/58056;
    mp_params('TR') = TR;
        
end

%% Reset All globals
% Parameters used in this code directly
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in find_a_working function
global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r'});
[theta,  r] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'TR'});
[TR] = params_group{:};    

%% Parse Model Controls
% Minimizer Controls
params_group = values(mp_controls, {'fl_max_trchk_perc_increase'});
[fl_max_trchk_perc_increase] = params_group{:};

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
if ~exist('ar_inc_amz','var')
    
    % initialize
    mn_inc = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_spouse_inc = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_a = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    
    % Txable Income at all state-space points
    for j=1:n_jgrid % Age
        for a=1:n_agrid % Assets
            for eta=1:n_etagrid % Productivity
                for educ=1:n_educgrid % Educational level
                    for married=1:n_marriedgrid % Marital status
                        for kids=1:n_kidsgrid % Number of kids

                            [inc,earn]=individual_income(j,a,eta,educ);
                            spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

                            mn_inc(j,a,eta,educ,married,kids) = inc;
                            mn_spouse_inc(j,a,eta,educ,married,kids) = (married-1)*spouse_inc*exp(eta_S_grid(eta));
                            mn_a(j,a,eta,educ,married,kids) = agrid(a);
                            
                        end
                    end
                end
            end
        end
    end
    
    % flatten the nd dimensional array
    ar_inc_amz = mn_inc(:);
    ar_spouse_inc_amz = mn_spouse_inc(:);
    ar_a_amz = mn_a(:);
    
end

%% B. Vectorized Solution for Optimal Check Adjustments

% B1. Anonymous Function where X is fraction of addition given bounds
fc_ffi_frac0t1_find_a_working = @(x) ffi_frac0t1_find_a_working_vec(...
    x, ...
    ar_a_amz, ar_inc_amz, ar_spouse_inc_amz, ...
    welf_checks, TR, r, fl_max_trchk_perc_increase);

% B2. Solve with Bisection
[~, ar_a_aux_bisec_amz] = ...
    ff_optim_bisec_savezrone(fc_ffi_frac0t1_find_a_working);

% B3. Reshape
mn_a_aux_bisec = reshape(ar_a_aux_bisec_amz, [n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid]);

%% C. Interpolate and Extrapolate All Age Check Values
% mn = nd dimensional array, mt = 2d
% ad1 = a is dimension 1

% C1. Change matrix order so asset becomes the first dimension
mn_v_ad1 = permute(V_ss, [2,1,3,4,5,6]);
mn_a_aux_bisec_ad1 = permute(mn_a_aux_bisec, [2,1,3,4,5,6]);

% C2. Reshape so that asset is dim 1, all other dim 2
mt_v_ad1 = reshape(mn_v_ad1, n_agrid, []);
mt_a_aux_bisec_ad1 = reshape(mn_a_aux_bisec_ad1, n_agrid, []);

% C3. Derivative dv/da
mt_dv_da_ad1 = diff(mt_v_ad1, 1)./diff(agrid);    

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

% C7. v(a_lower_idx,z) + slope*(fl_a_aux - fl_a)
ar_v_w_a_aux_j = ar_v_a_lower_idx + (ar_a_aux_bisec - agrid(ar_it_a_near_lower_idx)).*ar_deri_dev_dap;
mn_v_w_a_aux_ad1 = reshape(ar_v_w_a_aux_j, n_agrid, n_jgrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid);

% C8. Permute age as dim 1
V_W = permute(mn_v_w_a_aux_ad1, [2,1,3,4,5,6]);

%% D. Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_a4chk = strjoin(...
        ["Completed SNW_A4CHK_WRK_BISEC_VEC", ...
         ['welf_checks=' num2str(welf_checks)], ...
         ['TR=' num2str(TR)], ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_a4chk);
end

%% E. Compare Difference between V_ss and V_W

if (bl_print_a4chk_verbose)
    mn_V_gain_check = V_W - V_ss;
    mn_V_gain_frac_check = (V_W - V_ss)./V_ss;
    mp_container_map = containers.Map('KeyType','char', 'ValueType','any');
    mp_container_map('V_W') = V_W;
    mp_container_map('V_ss') = V_ss;
    mp_container_map('V_W_minus_V_ss') = mn_V_gain_check;
    mp_container_map('V_W_minus_V_ss_divide_V_ss') = mn_V_gain_frac_check;
    ff_container_map_display(mp_container_map);
end

end

function [ar_root_zero, ar_a_aux_amz] = ...
    ffi_frac0t1_find_a_working_vec(...
    ar_aux_change_frac_amz, ...
    ar_a_amz, ar_inc_amz, ar_spouse_amz, ...
    welf_checks, TR, r, fl_max_trchk_perc_increase)
    
    % Max A change to account for check
    fl_a_aux_max = TR*welf_checks*fl_max_trchk_perc_increase;    
    
    % Level of A change
    ar_a_aux_amz = ar_a_amz + ar_aux_change_frac_amz.*fl_a_aux_max;
    
    % Account for Interest Rates
    ar_r_gap = (1+r).*(ar_a_amz - ar_a_aux_amz);
    
    % Account for tax, inc changes by r
    ar_tax_gap = ...
          max(0, Tax(ar_inc_amz, ar_spouse_amz)) ...
        - max(0, Tax(ar_inc_amz - ar_a_amz*r + ar_a_aux_amz*r, ar_spouse_amz));
    
    % difference equation f(a4chkchange)=0
    ar_root_zero = TR*welf_checks + ar_r_gap - ar_tax_gap;
end