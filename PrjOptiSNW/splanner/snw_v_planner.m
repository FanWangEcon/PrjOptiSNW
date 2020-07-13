%% SNW_V_PLANNER Planner Last Year Value Given Income, Marriage, Age, Kids
%    Obtains the mass at each age, marriage status, child count, and income
%    level. And calculates planner value integrating over V_U and V_W,
%    where the state-space fall within the planner's
%    age/marriage/child-count/income-level group, and V_U and V_W are
%    specific to each check. Calculate the planner Value associated with
%    each check level. 
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
%    See also SNW_VU_VW_CHECKS
%

%%
function [varargout]=snw_v_planner(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==8)
        [mp_params, mp_controls, Phi_true, ap_ss, ...
            V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] = varargin{:};
    else
        error('Need to provide 8 parameter inputs');
    end

else
    close all;

    st_solu_type = 'bisec_vec';    
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');

    mp_params('TR') = 100/58056;
    mp_params('xi') = 0.5;
    mp_params('b') = 0;    
    mp_params('n_welfchecksgrid') = 3;
    n_incgrid = 5;
    mp_params('n_incgrid') = n_incgrid;
    mp_params('inc_grid') = linspace(0, 7, n_incgrid)'; % 7 refers to 7*58056=406392 dollars in 2012USD    
        
    
    mp_controls('bl_timer') = false;
    
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_vu_vw') = false;
    
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_vu_vw_verbose') = false;
    
    mp_controls('bl_print_v_planner') = true;
    mp_controls('bl_print_v_planner_verbose') = true;      
    
    % Calculate Values for All Checks    
    [~, ap_ss, ~, ...
        ~, ~, ~, ...
        V_W_allchecks, C_W_allchecks, V_U_allchecks, C_U_allchecks] = ...
        snw_vu_vw_checks(st_solu_type, mp_params, mp_controls);
    
    % Distribution
    [Phi_true] = snw_ds_main(mp_params, mp_controls);
    
end

%% Reset All globals
% Parameters used in this code directly
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in find_a_working function
global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option throw_in_ocean
global pi_unemp

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

params_group = values(mp_params, {'n_welfchecksgrid', 'scaleconvertor'});
[n_welfchecksgrid, scaleconvertor] = params_group{:};

params_group = values(mp_params, {'pi_unemp', 'n_incgrid', 'inc_grid'});
[pi_unemp, n_incgrid, inc_grid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_v_planner', 'bl_print_v_planner_verbose'});
[bl_print_v_planner, bl_print_v_planner_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tic
end

%% Wage Cut-Offs
cutoffs=wage_cutoffs(Phi_true);

%% Solve planner's problem

C_planner=NaN(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid);
V_planner=NaN(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid);
Phi_mass=NaN(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

for c_or_v=[1,2]
    
    if (bl_print_v_planner)
        disp(strcat(['SNW_V_PLANNER: c(1) or v(2):' num2str(c_or_v)]));
    end
    
    if (c_or_v == 1)
        W_allchecks = V_W_allchecks;
        U_allchecks = V_U_allchecks;
    elseif(c_or_v == 2)
        W_allchecks = C_W_allchecks;
        U_allchecks = C_U_allchecks;
    end
    
    for j=1:(n_jgrid-1) % Age
        for married=1:n_marriedgrid % Marital status
            for kids=1:n_kidsgrid % Number of kids
                for welf_checks=0:(n_welfchecksgrid-1)
                    for inc_group=1:n_incgrid
                        
                        if inc_group<n_incgrid
                            [V_or_C,Phi_mass(j,married,kids,inc_group)]=...
                                Planner(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),inc_grid(inc_group+1),...
                                U_allchecks,...
                                W_allchecks,...
                                ap_ss,...
                                cutoffs);
                            
                        elseif inc_group==n_incgrid
                            [V_or_C,Phi_mass(j,married,kids,inc_group)]=...
                                Planner(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),10E30,...
                                U_allchecks,...
                                W_allchecks,...
                                ap_ss,...
                                cutoffs);
                            
                        end
                        
                        if (c_or_v == 1)
                            V_planner(j,married,kids,welf_checks+1,inc_group) = V_or_C;
                        elseif(c_or_v == 2)
                            C_planner(j,married,kids,welf_checks+1,inc_group) = V_or_C;
                        end
                        
                    end
                end
            end
        end
        
        if (bl_print_v_planner)
            disp(strcat(['SNW_V_PLANNER: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid-1)]));
        end
        
    end
end

%% Output for computing optimal allocation
Output=NaN((n_jgrid-1)*n_marriedgrid*n_kidsgrid*n_welfchecksgrid*inc_group,9);
counter=0;
for j=1:(n_jgrid-1) % Age
    for married=1:n_marriedgrid % Marital status
        for kids=1:n_kidsgrid % Number of kids
            for welf_checks=0:(n_welfchecksgrid-1) % Number of welfare checks
                for inc_group=1:n_incgrid
                    
                    counter=counter+1;
                    
                    Output(counter,1)=17+j;
                    Output(counter,2)=married-1;
                    Output(counter,3)=kids-1;
                    Output(counter,4)=welf_checks;
                    
                    Output(counter,5)=inc_grid(inc_group)*scaleconvertor;
                    Output(counter,6)=10E30;
                    Output(counter,7)=Phi_mass(j,married,kids,inc_group);
                    Output(counter,8)=psi(j);
                    Output(counter,9)=V_planner(j,married,kids,welf_checks+1,inc_group);
                    Output(counter,10)=C_planner(j,married,kids,welf_checks+1,inc_group);
                    
                end
            end
        end
    end
    
end


%% Timing and Profiling End

if (bl_timer)
    toc
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_V_PLANNER", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_vu_vw_checks);    
end

%% Print
if (bl_print_v_planner_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('V_planner') = V_planner;
    mp_outcomes('C_planner') = C_planner;
    mp_outcomes('Phi_mass') = Phi_mass;
    mp_outcomes('Output') = Output;    
    ff_container_map_display(mp_outcomes);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = V_planner;
    elseif (it_k==2)
        ob_out_cur = C_planner;
    elseif (it_k==3)
        ob_out_cur = Phi_mass;
    elseif (it_k==4)
        ob_out_cur = Output;
    end
    varargout{it_k} = ob_out_cur;
end

end
