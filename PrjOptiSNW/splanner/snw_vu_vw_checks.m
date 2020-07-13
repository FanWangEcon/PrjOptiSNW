%% SNW_VU_VW_CHECKS Solve for V, Vumpe, and VU(check) VW(check)
%    What are the V and V_unemploy at all state space points?
%
%    For unemployed and employed individuals in 2020, what is the value of
%    of the check? Solve for V(s,check|working), V(s,check|unemployed). 
%
%    * ST_SOLU_TYPE: 'matlab_minimizer', 'grid_search', 'bisec_vec', three
%    alternative solution methods
%
%    [V_VFI, ap_VFI, cons_VFI] = SNW_VU_VW_CHECKS(ST_SOLU_TYPE) Solves
%    given solution method for value function, policy functions. Using
%    default parameters.
%
%    [V_VFI, ap_VFI, cons_VFI, V_VFI_unemp, ap_VFI_unemp, cons_VFI_unemp] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS) Solves given
%    solution method for value function, policy functions under both
%    employed and unemployed situations, using default parameters. Change
%    model parameters through MP_PARAMS.
%
%    [V_VFI, ap_VFI, cons_VFI, V_VFI_unemp, ap_VFI_unemp, cons_VFI_unemp,
%    V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS) Solves the
%    model for some number of checks for employed and unemployed. The
%    number of checks is set through the 'n_welfchecksgrid' parameters in
%    MP_PARAMS.
%
%    [V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] =
%    SNW_VU_VW_CHECKS(ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS, V_SS, V_UNEMP)
%    provide V_SS and V_UNEMP solved out elsewhere, and only get two
%    outputs out.
%
%    See also SNWX_A4CHK_WRK_BISEC_VEC_DENSE,
%    SNWX_A4CHK_WRK_BISEC_VEC_SMALL, SNW_A4CHK_WRK, SNW_A4CHK_WRK_BISEC
%

%%
function [varargout]=snw_vu_vw_checks(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==1)
        [st_solu_type] = varargin{:};
        mp_params = snw_mp_param('default_tiny');
        mp_controls = snw_mp_control('default_test');        
    elseif (length(varargin)==3)
        [st_solu_type, mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==5)
        [st_solu_type, mp_params, mp_controls, V_ss, cons_ss] = varargin{:};
    elseif (length(varargin)==7)
        [st_solu_type, mp_params, mp_controls, V_ss, cons_ss, V_unemp, cons_unemp] = varargin{:};
    else
        error('Need to provide 1/3/5/7 parameter inputs');
    end
    
else
    close all;
    
%     st_solu_type = 'matlab_minimizer';
    st_solu_type = 'bisec_vec';
%     st_solu_type = 'grid_search';
    
    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    
    % set Unemployment Related Variables
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
    n_welfchecksgrid=3; % Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars
    
    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_params('TR') = TR;
    mp_params('n_welfchecksgrid') = n_welfchecksgrid;
    
    % Solve for Unemployment Values
    mp_controls('bl_print_a4chk') = false;
        
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

params_group = values(mp_params, {'n_welfchecksgrid', 'xi','b'});
[n_welfchecksgrid, xi, b] = params_group{:};  

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_vu_vw', 'bl_print_vu_vw_verbose'});
[bl_print_vu_vw, bl_print_vu_vw_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tic
end

%% Compute V_ss If not Provided as a Parameter

% Compute policy functions in the event of unemployment. Required to compute V_U in the Planner's problem
if ~exist('V_ss','var')
    if (bl_print_vu_vw)
        disp('Compute value function and policy functions: V_ss')
    end
    if (strcmp(st_solu_type, 'matlab_minimizer')) 
        [V_ss,ap_ss,cons_ss,~] = snw_vfi_main(mp_params, mp_controls);
    elseif (strcmp(st_solu_type, 'grid_search'))
        [V_ss,ap_idx_VFI,cons_ss,~] = snw_vfi_main_grid_search(mp_params, mp_controls);
        ap_ss = ap_idx_VFI;
    elseif (strcmp(st_solu_type, 'bisec_vec'))
        [V_ss,ap_ss,cons_ss,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    else
        error(['SNW_VU_VW_CHECKS V, st_solu_type=' st_solu_type ' is not known'])
    end    
end

%% Compute V_unemp If not Provided as a Parameter

if (nargout ~= 3)
    
    % Compute policy functions in the event of unemployment. Required to compute V_U in the Planner's problem
    if ~exist('V_unemp','var')
        if (bl_print_vu_vw)
            disp('Compute value function and policy functions in the event of unemployment: V_unemp (given V_ss)')
        end
        
        if (strcmp(st_solu_type, 'matlab_minimizer'))
            [V_unemp,ap_unemp,cons_unemp,~] = snw_vfi_unemp(mp_params, mp_controls, V_ss);
        elseif (strcmp(st_solu_type, 'grid_search'))
            [V_unemp,ap_idx_VFI_unemp,cons_unemp,~]= snw_vfi_main_grid_search(mp_params, mp_controls, V_ss);
            ap_unemp = ap_idx_VFI_unemp;
        elseif (strcmp(st_solu_type, 'bisec_vec'))
            [V_unemp,ap_unemp,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
        else
            error(['SNW_VU_VW_CHECKS V_unemp, st_solu_type=' st_solu_type ' is not known'])
        end
        
    end
end

%% A1. Compute EMPLOYED Household-Head and Spousal Income

if (nargout ~= 3 && nargout ~= 6)
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
    
    % Print
    if (bl_print_vu_vw)
        disp(['SNW_VU_VW_CHECKS: Finished Employed income AMZ']);
    end
end

%% A2. Compute UNEMPLOYED Household-Head and Spousal Income 

if (nargout ~= 3 && nargout ~= 6)
    % initialize
    mn_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    mn_spouse_inc_unemp = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    
    % Txable Income at all state-space points
    for j=1:n_jgrid % Age
        for a=1:n_agrid % Assets
            for eta=1:n_etagrid % Productivity
                for educ=1:n_educgrid % Educational level
                    for married=1:n_marriedgrid % Marital status
                        for kids=1:n_kidsgrid % Number of kids
                            
                            [inc,earn]=individual_income(j,a,eta,educ,xi,b);
                            spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                            
                            mn_inc_unemp(j,a,eta,educ,married,kids) = inc;
                            mn_spouse_inc_unemp(j,a,eta,educ,married,kids) = (married-1)*spouse_inc*exp(eta_S_grid(eta));
                            
                        end
                    end
                end
            end
        end
    end
    
    % flatten the nd dimensional array
    ar_inc_unemp_amz = mn_inc_unemp(:);
    ar_spouse_inc_unemp_amz = mn_spouse_inc_unemp(:);
    
    % Print
    if (bl_print_vu_vw)
        disp(['SNW_VU_VW_CHECKS: Finished Unemployed income AMZ']);
    end
    
end

%% Compute Check V_W and check V_U Values

if (nargout ~= 3 && nargout ~= 6)
    
    V_W_allchecks=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
    V_U_allchecks=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
    
    C_W_allchecks=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
    C_U_allchecks=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
    
    if (bl_print_vu_vw)
        disp('Solve for V_W and V_U for different number of welfare checks')
    end
    
    for welf_checks=0:(n_welfchecksgrid-1)
        
        if (strcmp(st_solu_type, 'matlab_minimizer'))
            [V_W,C_W]=snw_a4chk_wrk(...
                welf_checks, V_ss, cons_ss, mp_params, mp_controls);
            [V_U,C_U]=snw_a4chk_unemp(...
                welf_checks, V_unemp, cons_unemp, mp_params, mp_controls);
        else
            % always use bisec_vec unless specified to use matlab_minimizer, no
            % grid_search option for this
            [V_W,C_W]=snw_a4chk_wrk_bisec_vec( ...
                welf_checks, V_ss, cons_ss, mp_params, mp_controls, ...
                ar_a_amz, ar_inc_amz, ar_spouse_inc_amz);
            [V_U,C_U]=snw_a4chk_unemp_bisec_vec( ...
                welf_checks, V_unemp, cons_unemp, mp_params, mp_controls, ...
                ar_a_amz, ar_inc_unemp_amz, ar_spouse_inc_unemp_amz);
        end
        
        % Update Collection:
        V_W_allchecks(:,:,:,:,:,:,welf_checks+1) = V_W;
        V_U_allchecks(:,:,:,:,:,:,welf_checks+1) = V_U;
        C_W_allchecks(:,:,:,:,:,:,welf_checks+1) = C_W;
        C_U_allchecks(:,:,:,:,:,:,welf_checks+1) = C_U;
        
        % Print
        if (bl_print_vu_vw)
            disp(strcat(['SNW_VU_VW_CHECKS: Finished checks:' ...
                num2str(welf_checks) ' of ' num2str(n_welfchecksgrid-1)]));
        end
        
    end
    
end

%% Timing and Profiling End

if (bl_timer)
    toc
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_VU_VW_CHECKS", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_vu_vw_checks);    
end

%% Return

if (length(varargin)==5)
    varargout = cell(nargout,0);
    for it_k = 1:nargout
        if (it_k==1)
            ob_out_cur = V_W_allchecks;
        elseif (it_k==2)
            ob_out_cur = C_W_allchecks;
        elseif (it_k==3)
            ob_out_cur = V_U_allchecks;
        elseif (it_k==4)
            ob_out_cur = C_U_allchecks;
        end
        varargout{it_k} = ob_out_cur;
    end    
else
    varargout = cell(nargout,0);
    for it_k = 1:nargout
        if (it_k==1)
            ob_out_cur = V_ss;
        elseif (it_k==2)
            ob_out_cur = ap_ss;
        elseif (it_k==3)
            ob_out_cur = cons_ss;
        elseif (it_k==4)
            ob_out_cur = V_unemp;
        elseif (it_k==5)
            ob_out_cur = ap_unemp;
        elseif (it_k==6)
            ob_out_cur = cons_unemp;
        elseif (it_k==7)
            ob_out_cur = V_W_allchecks;
        elseif (it_k==8)
            ob_out_cur = C_W_allchecks;
        elseif (it_k==9)
            ob_out_cur = V_U_allchecks;
        elseif (it_k==10)
            ob_out_cur = C_U_allchecks;
        end
        varargout{it_k} = ob_out_cur;
    end
end

end
