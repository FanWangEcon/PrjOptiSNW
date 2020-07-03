%% SNW_A4CHK_WRK solves for Asset Position Corresponding to Check Level
%    What is the value of a check? From the perspective of the value
%    function? We have Asset as a state variable, in a cash-on-hand sense,
%    how much must the asset (or think cash-on-hand) increase by, so that
%    it is equivalent to providing the household with a check?
%
%    * WELF_CHECKS integer the number of checks
%    * TR float the value of each check
%    * V_SS ndarray the value matrix along standard state-space dimensions:
%    (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid)
%    * MP_PARAMS map with model parameters
%    * MP_CONTROLS map with control parameters
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_A4CHK_WRK(MP_PARAMS) invoke
%    model with externally set parameter map MP_PARAMS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_A4CHK_WRK(MP_PARAMS,
%    MP_CONTROLS) invoke model with externally set parameter map MP_PARAMS
%    as well as control mpa MP_CONTROLS.
%
%    See also SNWX_A4CHK_WRK, FIND_A_WORKING
%

%%
function [V_W, exitflag_fsolve]=snw_a4chk_wrk(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==3)
        [welf_checks, TR, V_ss] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==5)
        [welf_checks, TR, V_ss, mp_params, mp_controls] = varargin{:};
    end
    
else
    close all;
    
    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    [V_ss,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    welf_checks = 2;
    TR = 100/58056;
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in find_a_working function
global theta r agrid epsilon eta_grid SS

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r'});
[theta,  r] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_grid'});
[agrid, eta_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

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

%% Loop over states and find A state for a Particular Check Level

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
exitflag_fsolve=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids
                        
                        % Find value of assets that approximates the value of the welfare checks
                        x0=agrid(a)+TR*welf_checks; % Initial guess for a
                        
                        [a_aux,~,exitflag_fsolve(j,a,eta,educ,married,kids)]=fsolve(@(x)find_a_working(x,j,a,eta,educ,married,kids,TR,welf_checks),x0,options2);
                        
                        if a_aux<0
                            a_aux=0;
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
                        
                        V_W(j,a,eta,educ,married,kids)=vals(1)*V_ss(j,inds(1),eta,educ,married,kids)+vals(2)*V_ss(j,inds(2),eta,educ,married,kids);
                        
                    end
                end
            end
        end
    end
    
    if (bl_print_a4chk)
        disp(strcat(['SNW_A4CHK_WRK: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end
    
end


end