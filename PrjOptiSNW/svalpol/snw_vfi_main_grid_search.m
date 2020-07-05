%% SNW_VFI_MAIN_GRID_SEARCH Solves Policy/Value Function SNW (Grid Search)
%    Given parameters, iterate over life cycle, given age, marital status,
%    education level and child count, as well as persistent productivity
%    shock process, solve for optimal dynamic savings choices given
%    expectation of kid count transition and productivity shock transition.
%
%    This is the grid search solution. Loop over the state space,
%    vectorized solve over the choice space. Choice space is a level grid
%    shared by all state-space points.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN_GRID_SEARCH(MP_PARAMS) invoke
%    model with externally set parameter map MP_PARAMS.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN_GRID_SEARCH(MP_PARAMS,
%    MP_CONTROLS) invoke model with externally set parameter map MP_PARAMS
%    as well as control mpa MP_CONTROLS.
%
%    See also SNW_VFI_MAIN, SNW_VFI_BISEC_VEC, SNW_VFI_UNEMP,
%    SNW_MP_CONTROL, SNW_MP_PARAM
%

%%
function [V_VFI,ap_idx_VFI,cons_VFI,exitflag_VFI]=snw_vfi_main_grid_search(varargin)

%% Default and Parse 
if (~isempty(varargin))
    
    if (length(varargin)==1)
        [mp_params] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==3)
        [mp_params, mp_controls, V_VFI_POSTSHOCK] = varargin{:};        
    end
    
else
    
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global beta theta r agrid epsilon SS pi_eta pi_kids psi n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in functions that are called by this code
global gamma g_n g_cons a2 cons_allocation_rule jret
% July 1st new parameters
global eta_H_grid eta_S_grid Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {...
    'gamma', 'beta', 'theta', 'cons_allocation_rule', ...
    'r', 'g_n', 'g_cons', 'a2', 'jret'});
[gamma, beta, theta, cons_allocation_rule, ...
    r, g_n, g_cons, a2, jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'psi'});
[pi_eta, pi_kids, psi] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

%% Parse Model Controls
% Minimizer Controls
params_group = values(mp_controls, ...
    {'A_aux', 'B_aux', ...
     'Aeq', 'Beq',...
     'nonlcon', 'options', 'options2'});
[A_aux, B_aux, ...
    Aeq, Beq, ...
    nonlcon, options, options2] = params_group{:};

% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_vfi', 'bl_print_vfi_verbose'});
[bl_print_vfi, bl_print_vfi_verbose] = params_group{:};

% Store Controls
params_group = values(mp_controls, {'bl_vfi_store_all'});
[bl_vfi_store_all] = params_group{:};

% unemployment parameters under covid
if (length(varargin)==3)    
    params_group = values(mp_params, {'xi','b'});
    [xi, b] = params_group{:};    
end

%% Timing and Profiling Start
if (bl_timer)
    tic
end

%% Solve optimization problem

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_idx_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% if (bl_vfi_store_all)
%     y_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     tax_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     SS_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% end

exitflag_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

if (length(varargin)==3)
    ar_j_seq = 1:n_jgrid;
else
    ar_j_seq = n_jgrid:(-1):1; % Age
end

% Solve for value function and policy functions by means of backwards induction
for j=ar_j_seq % Age

    % Solving for VFI with or Without Shock
    if (length(varargin)==3)
        V_VFI_FUTURE = V_VFI_POSTSHOCK;
    else
        V_VFI_FUTURE = V_VFI;
    end
    
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids
                        
                        if j==n_jgrid
                            
                            ap_idx_VFI(j,a,eta,educ,married,kids)=1;
                            if (length(varargin)==3)
                                cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_idx_VFI(j,a,eta,educ,married,kids),xi,b);
                            else
                                cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_idx_VFI(j,a,eta,educ,married,kids));
                            end
                            
                            if cons_VFI(j,a,eta,educ,married,kids)<=0
                                disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                error('Non-positive consumption')
                            end
                            
                            V_VFI(j,a,eta,educ,married,kids)=utility_grid_search(cons_VFI(j,a,eta,educ,married,kids),married,kids);
                            
                        else
                            
                            if (length(varargin)==3)
                                c_aux=consumption_grid_search(j,a,eta,educ,married,kids,1:n_agrid,xi,b);
                            else
                                c_aux=consumption_grid_search(j,a,eta,educ,married,kids,1:n_agrid);
                            end                            
                                                        
                            cont=zeros(n_agrid,n_etagrid);
                            for kidsp=1:n_kidsgrid
                                cont(1:n_agrid,1:n_etagrid)=cont(1:n_agrid,1:n_etagrid)+pi_kids(kids,kidsp,j,educ,married)*reshape(V_VFI_FUTURE(j+1,1:n_agrid,1:n_etagrid,educ,married,kidsp),n_agrid,n_etagrid);
                            end
                            
                            cont=pi_eta(eta,:)*cont(1:n_agrid,:)';
                            
                            V_VFI_aux=utility_grid_search(c_aux,married,kids)+beta*psi(j)*cont';
                            
                            [max_val,max_ind]=max(V_VFI_aux);
                            
                            ap_idx_VFI(j,a,eta,educ,married,kids)=max_ind;
                            if (length(varargin)==3)
                                cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_idx_VFI(j,a,eta,educ,married,kids),xi,b);
                            else
                                cons_VFI(j,a,eta,educ,married,kids)=consumption_grid_search(j,a,eta,educ,married,kids,ap_idx_VFI(j,a,eta,educ,married,kids));
                            end                            
                            
                            V_VFI(j,a,eta,educ,married,kids)=max_val;
                            
                            if cons_VFI(j,a,eta,educ,married,kids)<=0
                                disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                error('Non-positive consumption')
                            end
                            
                        end
                        
                    end
                end
            end
        end
    end
    
    if (bl_print_vfi)
        disp(strcat(['SNW_VFI_MAIN_GRID_SEARCH: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end
    
end

%% Timing and Profiling End
if (bl_timer)
    toc;
    if (length(varargin)==3)
        st_complete_vfi = strjoin(...
            ["Completed SNW_VFI_MAIN_GRID_SEARCH 1 PERIOD UNEMP SHK", ...
             ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
             ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
            ], ";");
    else
        st_complete_vfi = strjoin(...
            ["Completed SNW_VFI_MAIN_GRID_SEARCH", ...
             ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
             ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
            ], ";");
    end
    disp(st_complete_vfi);
end

end
