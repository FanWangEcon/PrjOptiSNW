%% SNW_EVUVW20_JAEEMK Solves for EVUVW given Unemployment Risk
%    Solve for Value working and unemployed and consumption working and
%    unemployed, and get their average based on unemployment shock. This
%    is the average value in 2020. For one Check, the check number set by
%    WELF_CHECKS.
%
%    [ev20_jaeemk, ec20_jaeemk] = SNW_EVUVW20_JAEEMK(WELF_CHECKS,
%    ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS, V_SS, V_UNEMP) provide V_SS and
%    V_UNEMP solved out elsewhere, and only get two outputs out.
%
%    See also SNW_EVUVW19_JYMK, SNW_EVUVW19_JAEEMK, SNW_EVUVW20_JAEEMK,
%    SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw20_jaeemk(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==9)
        [welf_checks, st_solu_type, mp_params, mp_controls, ...
            V_ss, cons_ss, ...
            V_unemp, cons_unemp, ...
            mp_precompute_res] = varargin{:};
    else
        error('Need to provide 9 parameter inputs');
    end
    
else
    close all;
        
%     st_solu_type = 'matlab_minimizer';
    st_solu_type = 'bisec_vec';
%     st_solu_type = 'grid_search';
    
    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
%     mp_params = snw_mp_param('default_moredense');
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
    mp_controls('bl_print_ds_verbose') = false;
    
    % Solve the Model to get V working and unemployed
    [V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [V_unemp,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);        
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
    
    % Get Matrixes
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid'};        
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);
    
end

%% Parse Pre-Computes
params_group = values(mp_precompute_res, {'ar_a_amz'});
[ar_a_amz] = params_group{:};
params_group = values(mp_precompute_res, {'ar_inc_amz', 'ar_inc_unemp_amz'});
[ar_inc_amz, ar_inc_unemp_amz] = params_group{:};
params_group = values(mp_precompute_res, {'ar_spouse_inc_amz', 'ar_spouse_inc_unemp_amz'});
[ar_spouse_inc_amz, ar_spouse_inc_unemp_amz] = params_group{:};
params_group = values(mp_precompute_res, {'ref_earn_wageind_grid'});
[ref_earn_wageind_grid] = params_group{:};

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp'});
[pi_unemp] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_evuvw20_jaeemk', 'bl_print_evuvw20_jaeemk_verbose'});
[bl_print_evuvw20_jaeemk, bl_print_evuvw20_jaeemk_verbose] = params_group{:};

%% Compute Check V_W and check V_U Values
    
if (bl_print_evuvw20_jaeemk)
    disp('Solve for V_W and V_U for different number of welfare checks')
end

if (strcmp(st_solu_type, 'matlab_minimizer'))
    [V_U,C_U]=snw_a4chk_unemp(...
        welf_checks, V_unemp, cons_unemp, mp_params, mp_controls);            
    [V_W,C_W]=snw_a4chk_wrk(...
        welf_checks, V_ss, cons_ss, mp_params, mp_controls);
else
    % always use bisec_vec unless specified to use matlab_minimizer, no
    % grid_search option for this
    [V_U,C_U]=snw_a4chk_unemp_bisec_vec( ...
        welf_checks, V_unemp, cons_unemp, mp_params, mp_controls, ...
        ar_a_amz, ar_inc_unemp_amz, ar_spouse_inc_unemp_amz);            
    [V_W,C_W]=snw_a4chk_wrk_bisec_vec( ...
        welf_checks, V_ss, cons_ss, mp_params, mp_controls, ...
        ar_a_amz, ar_inc_amz, ar_spouse_inc_amz);
end

%% Timing and Profiling Start

if (bl_timer)
    tm_start = tic;
end

%% Take Expectation

% EV Store, expected check values
ev20_jaeemk=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ec20_jaeemk=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Loop Over States and Pre-Store
for j=1:n_jgrid
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid
                        
                        % 1. Get Wage Index
                        wage_ind = ref_earn_wageind_grid(j,a,eta,educ,married,kids);
                        
                        % 2. Expected Value of Checks
                        ev20_jaeemk(j,a,eta,educ,married,kids)=...
                            pi_unemp(j,wage_ind)*V_U(j,a,eta,educ,married,kids)...
                            +(1-pi_unemp(j,wage_ind))*V_W(j,a,eta,educ,married,kids);
                        
                        ec20_jaeemk(j,a,eta,educ,married,kids)=...
                            pi_unemp(j,wage_ind)*C_U(j,a,eta,educ,married,kids)...
                            +(1-pi_unemp(j,wage_ind))*C_W(j,a,eta,educ,married,kids);  
                        
                    end
                end
            end
        end
    end
end   

%% Timing and Profiling End

if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_EVUVW20_JAEEMK", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['timeEUEC=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);    
end

%% Print 
if (bl_print_evuvw20_jaeemk_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('ev20_jaeemk') = ev20_jaeemk;
    mp_outcomes('ec20_jaeemk') = ec20_jaeemk;
    mp_outcomes('V_W') = V_W;
    mp_outcomes('C_W') = C_W;    
    mp_outcomes('V_U') = V_U;
    mp_outcomes('C_U') = C_U;    
    ff_container_map_display(mp_outcomes, 9, 9);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev20_jaeemk;
    elseif (it_k==2)
        ob_out_cur = ec20_jaeemk;        
    elseif (it_k==3)
        ob_out_cur = V_W;
    elseif (it_k==4)
        ob_out_cur = V_U;
    elseif (it_k==5)
        ob_out_cur = C_W;
    elseif (it_k==6)
        ob_out_cur = C_U;
    end
    varargout{it_k} = ob_out_cur;
end

end
