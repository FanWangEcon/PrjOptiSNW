%% SNW_VFI_MAIN_BISEC_VEC_STIMULUS Solves Policy/Value Function SNW (Bisection Vectorized)
%    This is a wrapper for the policy function solution function
%    SNW_VFI_MAIN_BISEC_VEC. The wrapper provides a fourth input to that
%    function that provides income, marital and kids specific stimulus
%    check provisions. 
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC_STIMULUS(MP_PARAMS) invoke model with
%    externally set parameter map MP_PARAMS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC_STIMULUS(MP_PARAMS, MP_CONTROLS) invoke model
%    with externally set parameter map MP_PARAMS as well as control mpa
%    MP_CONTROLS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC_STIMULUS(MP_PARAMS, MP_CONTROLS, V_VFI_FIX)
%    provides existing value function. Suppose there is sudden shock, but
%    future value is preserved after one period. So now we have new value
%    that is specific to this period, that is the output V_VFI, the input
%    V_VFI_FIX is the value for all future periods. When this program is
%    called with V_VFI_FIX, the resource equation will use the unemployment
%    shock information.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] =
%    SNW_VFI_MAIN_BISEC_VEC_STIMULUS(MP_PARAMS, MP_CONTROLS, V_VFI_FIX,
%    CL_FC_WELFARE_CHECKS) CL_FC_WELFARE_CHECKS provides a cell array of
%    anonymous functions where each cell corresponds to a particular
%    child-count and marital status group, and each anonymous function is a
%    function of a single variable, aggregate household realized income,
%    and the output of the function is the amount of stimulus check, in
%    dollar units (which has to be converted to model units).
%
%    See also SNWX_VFI_MAIN_BISEC_VEC_STIMULUS, SNW_VFI_MAIN_BISEC_VEC
%

%%
function [varargout]=snw_vfi_main_bisec_vec_stimulus(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==3)
        [mp_params, mp_controls, v_vfi_postshock] = varargin{:};
    end
    
else
    
    clc;
    mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
    mp_more_inputs('st_edu_simu_type') = 'both';
    mp_params = snw_mp_param('default_tiny', false, 'tauchen', false, 8, 8, mp_more_inputs);    
%     mp_more_inputs('st_edu_simu_type') = 'low';
%     mp_params = snw_mp_param('default_tiny_e1l', false, 'tauchen', false, 8, 8, mp_more_inputs);
%     mp_more_inputs('st_edu_simu_type') = 'high';
%     mp_params = snw_mp_param('default_tiny_e2h', false, 'tauchen', false, 8, 8, mp_more_inputs);
   
    mp_controls = snw_mp_control('default_test');
    false;
    
    % Solve steady-state problem
    [v_vfi_postshock,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    
    % Covid-Period 
    mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');    
    mp_params('xi') = 0.5;
    mp_params('b') = 0;

end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, ...
    {'fl_stimulus_adult_first', 'fl_stimulus_child_first', ...
    'fl_stimulus_adult_second', 'fl_stimulus_child_second'});
[fl_stimulus_adult_first, fl_stimulus_child_first, ...
    fl_stimulus_adult_second, fl_stimulus_child_second] = params_group{:};

%% Collect stimulus check functions
% Collection cell
cl_fc_welfare_checks = cell([n_kidsgrid, n_marriedgrid]);

% Store functions
for married=1:n_marriedgrid
    
    % Marital status
    if (married == 1)
        bl_marital = false;
    else
        bl_marital = true;
    end
    
    for kids=1:n_kidsgrid
        
        % Kids count
        it_kids=kids-1;
        
        % Run functions
        fc_stimulus_check_firstsecond = snw_stimulus_checks(it_kids, bl_marital, ...
            fl_stimulus_adult_first, fl_stimulus_child_first, ...
            fl_stimulus_adult_second, fl_stimulus_child_second);
        
        % Store results
        cl_fc_welfare_checks{kids, married} = fc_stimulus_check_firstsecond;
        
    end
end

%% Call to solve the VFI problem
% calling unemployment problem, however, by adjust xi and b, unemployment
% result can be the same as employed result. 
[V_VFI_unemp, ap_VFI_unemp, cons_VFI_unemp, mp_valpol_more_unemp] = ...
    snw_vfi_main_bisec_vec(mp_params, mp_controls, v_vfi_postshock, cl_fc_welfare_checks);

%% Return 
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_VFI_unemp;
    elseif (it_k==2)
        ob_out_cur = ap_VFI_unemp;
    elseif (it_k==3)
        ob_out_cur = cons_VFI_unemp;
    elseif (it_k==4)
        ob_out_cur = mp_valpol_more_unemp;
    end
    varargout{it_k} = ob_out_cur;
end

end
