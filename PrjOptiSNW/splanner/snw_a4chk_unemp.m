%% SNW_A4CHK_UNEMP (loop fzero) Asset Position Corresponding to Check Level
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
%    * WELF_CHECKS integer the number of checks
%    * V_UNEMP ndarray the value matrix along standard state-space
%    dimensions:
%    (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid) for
%    value with associated one period unemployment shock
%    * MP_PARAMS map with model parameters
%    * MP_CONTROLS map with control parameters
%
%    [V_U, EXITFLAG_FSOLVE] = SNW_A4CHK_UNEMP(WELF_CHECKS, V_UNEMP)
%    solves for working value given V_UNEMP value function results, value
%    function given one period unemployment shock. WELF_CHECKS is the
%    number of checks. TR, the value of each check, XI, intensity of
%    of check WELF_CHECKS, and given the value of each check equal to TR.
%    Unemployment XI and B variable will be automatically grabbed out from the , XI is the intensity of
%    unemployment proportional earning shock, B is the replacement rate.
%
%    [V_U, EXITFLAG_FSOLVE] = SNW_A4CHK_UNEMP(WELF_CHECKS, TR, V_SS, XI, B)
%    solves for working value given V_SS value function results, for number
%    of check WELF_CHECKS, and given the value of each check equal to TR,
%    with unemployment XI and B variable, XI is the intensity of
%    unemployment proportional earning shock, B is the replacement rate.
%
%    [V_W, EXITFLAG_FSOLVE] = SNW_A4CHK_UNEMP(WELF_CHECKS, TR, V_SS, XI, B,
%    MP_PARAMS, MP_CONTROLS) control parameters and controls.
%
%    See also SNWX_A4CHK_WRK_SMALL, SNWX_A4CHK_WRK_DENSE,
%    SNW_A4CHK_UNEMP_BISEC, SNW_A4CHK_UNEMP_BISEC_VEC, FIND_A_WORKING
%

%%
function [V_U, C_U]=snw_a4chk_unemp(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==3)
        [welf_checks, V_unemp, cons_unemp] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==5)
        [welf_checks, V_unemp, cons_unemp, mp_params, mp_controls] = varargin{:};
    else
        error('Need to provide 3/5 parameter inputs');
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
% globals = who('global');
% clear(globals{:});
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

params_group = values(mp_params, {'TR', 'xi','b'});
[TR, xi, b] = params_group{:};

%% Parse Model Controls
% Minimizer Controls
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

%% Loop over states and find A state for a Particular Check Level

V_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
exitflag_fsolve=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids

                        % Find value of assets that approximates the value of the welfare checks
                        x0=agrid(a)+TR*welf_checks; % Initial guess for a

                        [a_aux,~,exitflag_fsolve(j,a,eta,educ,married,kids)]=fsolve(@(x)find_a_unemp(x,j,a,eta,educ,married,kids,TR,xi,b,welf_checks),x0,options2);

                        if a_aux<0
                            disp(a_aux)
                            error('Check code! Should not allow for negative welfare checks')
                        elseif a_aux>agrid(n_agrid)
                            a_aux=agrid(n_agrid);
                        end

                        % Linear interpolation
                        ind_aux=find(agrid<=a_aux,1,'last');

                        if a_aux==0
                            inds(1)=1;
                            inds(2)=1;
                            vals(1)=1;
                            vals(2)=0;

                        elseif a_aux==agrid(n_agrid)
                            inds(1)=n_agrid;
                            inds(2)=n_agrid;
                            vals(1)=1;
                            vals(2)=0;

                        else
                            inds(1)=ind_aux;
                            inds(2)=ind_aux+1;
                            vals(1)=1-((a_aux-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                            vals(2)=1-vals(1);

                        end

                        V_U(j,a,eta,educ,married,kids)=vals(1)*V_unemp(j,inds(1),eta,educ,married,kids)+vals(2)*V_unemp(j,inds(2),eta,educ,married,kids);
                        C_U(j,a,eta,educ,married,kids)=vals(1)*(cons_unemp(j,inds(1),eta,educ,married,kids)/(married+kids-1))+vals(2)*(cons_unemp(j,inds(2),eta,educ,married,kids)/(married+kids-1));
                    end
                end
            end
        end
    end

    if (bl_print_a4chk)
        disp(strcat(['SNW_A4CHK_UNEMP: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end

end

%% Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_a4chk = strjoin(...
        ["Completed SNW_A4CHK_UNEMP", ...
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
