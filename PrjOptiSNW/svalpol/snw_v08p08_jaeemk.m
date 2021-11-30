%% SNW_V08P08_JAEEMK Value and Policy in 2008 Given 09 Unemployment Uncertainty
%    Solve for 2008 policy and value functions. 2009 arrives with some
%    probability of unemployment shock, that is specific to education and
%    age groups. The 2009 shock is not a MIT shock, and households in 2008
%    make optimal consumption and savings decisions taking the one-period
%    expected 2009 unemployment shock into consideration.
%
%    Below, to solve for the 2008 problem, we proceed in three steps:
%
%    1. Provide as inputs 2009 V for employed and 2009 V for unemployed.
%    Since the 2009 shock lasts for only 1 period, the 2009 employed V and
%    policy function follows as in steady state.
%
%    2. Household EV' in 2009, integrating over two employment states,
%    where X includes age and education
%       - EV^{09;X} = P(U|X)*V^{U,09} + (1-P(U|X))*V^{E,09}
%
%    3. Solve for Value/Policy in 2008 given EV'
%       - V^{08} = U(C^{08}) + \beta*EV^{09;X}(A^{08,\prime})
%
%    Note that in our structure, the EV09 is solved once, checks are given
%    ex-ante of the realization of shocks, not ex-post. This is a key
%    difference with respect to the COVID19 checks, where checks are
%    received ex-post of shocks, but determined based on ex-ante
%    information. Here, checks are also based on ex-ante information.
%
%    Note that this function is not named well, would be better named as
%    SNW_VFI_MAIN_BISEC_VEC_2008 to be consistent with over VFI functions.
%
%    [V_2008, AP_2008, CONS_2008, EV_EMPSHK_2009] =
%    SNW_v08p08_jaeemk(MP_PARAMS, MP_CONTROLS, V_EMP_2009, V_UNEMP_2009)
%    provide V_EMP_2009 and V_UNEMP_2009 pre-solved.
%
%    See also SNWX_V0808_JAEEMK, SNW_V0808_JAEEMK, SNWX_UNEMP_2008, SNW_UNEMP_2008, SNW_UNEMP_2008
%

%%
function [varargout]=snw_v08p08_jaeemk(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==4)
        [mp_params, mp_controls, V_emp_2009, V_unemp_2009] = varargin{:};
    else
        error('Need to provide 4 parameter inputs');
    end

else
    clc;
    close all;

    % Solve the VFI Problem and get Value Function
    mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
    mp_more_inputs('fl_ss_non_college') = 0.225;
    mp_more_inputs('fl_ss_college') = 0.271;
    mp_more_inputs('fl_scaleconvertor') = 54831;
    st_param_group = 'default_small';
%     st_param_group = 'default_dense';
    mp_params = snw_mp_param(st_param_group, false, 'tauchen', false, 8, 8, mp_more_inputs);

    mp_controls = snw_mp_control('default_test');

    % set Unemployment Related Variables
    xi=0.532; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0.37992; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)

    mp_params('xi') = xi;
    mp_params('b') = b;

    % Solve for Unemployment Values
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_v08p08_jaeemk') = true;
    mp_controls('bl_print_v08p08_jaeemk_verbose') = true;

    % Solve the Model to get V working and unemployed
    [V_ss,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);

    % 2020 V and C same as V_SS and cons_ss if tax the same
    V_emp_2009 = V_ss;

    % 2020 Unemployed results
    % a2_covid_year  means a2 crisis year, this tax rate set to normal for 08 crisis.
    mp_params('a2_covidyr') = mp_params('a2_greatrecession_2009');
    [V_unemp_2009,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);

end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp_2009_edu_age'});
[pi_unemp_2009_edu_age] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_v08p08_jaeemk', 'bl_print_v08p08_jaeemk_verbose'});
[bl_print_v08p08_jaeemk, bl_print_v08p08_jaeemk_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Solve 08 Policy and Value given 09 Unemployment shock
% 09 has two states, and two associated values, employed or NOT
% 09 shock is NOT MIT shock, expected from 08 perspective
% resolve the 08 problem given shock expectation, unemployment probability.

% EV Store, expected check values
% EV_empshk_2009 integrates by states the unemployment shock uncertainty.
EV_empshk_2009 = NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Loop Over States and Pre-Store
for j=1:n_jgrid
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid

                        % 1. Get unemployment probability conditional on age and edu
                        % columns are edu: 1st column high school; 2nd col college
                        % rows are age groups: 18--24; 25-54; 55 to 65
                        if (j<=7)
                            it_age_grp = 1;
                        elseif (j<=37)
                            it_age_grp = 2;
                        else
                            it_age_grp = 3;
                        end
                        fl_unemp_edu_age = pi_unemp_2009_edu_age(it_age_grp, educ);

                        % 2. Expected Value at each state, integrating over employment shock
                        % employment shock is not observed.
                        EV_empshk_2009(j,a,eta,educ,married,kids) = ...
                            fl_unemp_edu_age*V_unemp_2009(j,a,eta,educ,married,kids)...
                            +(1 - fl_unemp_edu_age)*V_emp_2009(j,a,eta,educ,married,kids);

                    end
                end
            end
        end
    end
end

%% 2008 Optimal choices

% a2_covidyr, here meant to be tax in 2008, during disbursment of checks
% when there is a 3rd V parameter in SNW_VFI_MAIN_BISEC_VEC call, needs to
% specify A2_COVIDYR as the tax rate for this particular year.
mp_params('a2_covidyr') = mp_params('a2_bushchkyr_2008');

xi = mp_params('xi');
b = mp_params('b');
mp_params('xi') = 1;
mp_params('b') = 0;
[V_2008, ap_2008, cons_2008, ~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, EV_empshk_2009);
mp_params('xi') = xi;
mp_params('b') = b;

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_V08P08_JAEEMK", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_v08p08_jaeemk_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('V_2008') = V_2008;
    mp_outcomes('ap_2008') = ap_2008;
    mp_outcomes('cons_2008') = cons_2008;
    ff_container_map_display(mp_outcomes, 9, 9);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_2008;
    elseif (it_k==2)
        ob_out_cur = ap_2008;
    elseif (it_k==3)
        ob_out_cur = cons_2008;
    elseif (it_k==4)
        ob_out_cur = EV_empshk_2009;
    end
    varargout{it_k} = ob_out_cur;
end

end
