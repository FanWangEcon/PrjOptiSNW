%% SNW_VFI_UNEMP Solves Policy/Value Function SNW (Loop)
%    What is the value given an one period unemployment shock? We have the
%    same value in the future after this period, but for the current
%    period, there is a sudden bad shock. Given this sudden bad shock,
%    current resources are reduced, current consumption as well as savigs
%    choices are both likely to go down. This is the fmincon results
%
%    % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
%    mp_params('xi')=0.5; 
%    % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
%    mp_params('b')=0; 
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN(V_SS, MP_PARAMS)
%    invoke model with externally set parameter map MP_PARAMS, given V_SS
%    future value solved before the shock.
%
%    [V_VFI,AP_VFI,CONS_VFI,EXITFLAG_VFI] = SNW_VFI_MAIN(V_SS. MP_PARAMS,
%    MP_CONTROLS) invoke model with externally set parameter map MP_PARAMS
%    as well as control mpa MP_CONTROLS.
%
%    See also SNW_VFI_MAIN, SNW_MP_CONTROL, SNW_MP_PARAM
%

%%
function [V_VFI,ap_VFI,cons_VFI,exitflag_VFI]=snw_vfi_unemp(varargin)

%% Default and Parse 
if (~isempty(varargin))
    
    if (length(varargin)==2)
        [V_ss, mp_params] = varargin{:};
        mp_controls = snw_mp_control('default_base');
    elseif (length(varargin)==3)
        [V_ss, mp_params, mp_controls] = varargin{:};
    end
    
else
    
    mp_params = snw_mp_param('default_small');
    mp_controls = snw_mp_control('default_test');
    
    [V_ss,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    
    % Solve for Value of One Period Unemployment Shock
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)

    mp_params('xi') = xi;
    mp_params('b') = b;
    
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

params_group = values(mp_params, {'xi','b'});
[xi, b] = params_group{:};    

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

%% Timing and Profiling Start
if (bl_timer)
    tic
end

%% Solve optimization problem

V_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ap_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
cons_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% if (bl_vfi_store_all)
%     y_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     tax_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
%     SS_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% end

exitflag_VFI=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Solve for value function and policy functions by means of backwards induction
for j=1:n_jgrid % Age
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids
                        
                        if j<n_jgrid
                            
                            % Solve for next period's assets
                            x0=agrid(a); % Initial guess for ap
                            
                            amin=0;
                            
                            [inc,earn]=individual_income(j,a,eta,educ,xi,b);
                            spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                            
                            amax = min(agrid(end), ...
                                (1+r)*(agrid(a) ...
                                + Bequests*(bequests_option-1)) ...                                
                                + (epsilon(j,educ)*theta*exp(eta_H_grid(eta)))*(xi+b*(1-xi)) ...
                                + SS(j,educ) ...
                                + (married-1)*spouse_inc*exp(eta_S_grid(eta)) ...
                                - max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)))));
                            
                            [ap_aux,~,exitflag_VFI(j,a,eta,educ,married,kids)]=fmincon(@(x)value_func_aux_unemp(x,j,a,eta,educ,married,kids,V_ss,xi,b),x0,A_aux,B_aux,Aeq,Beq,amin,amax,nonlcon,options);
                            
                            ind_aux=find(agrid<=ap_aux,1,'last');
                            
                            % Linear interpolation
                            if ap_aux==0
                                inds(1)=1;
                                inds(2)=1;
                                vals(1)=1;
                                vals(2)=0;
                                
                            elseif ap_aux==agrid(n_agrid)
                                inds(1)=n_agrid;
                                inds(2)=n_agrid;
                                vals(1)=1;
                                vals(2)=0;
                                
                            else
                                inds(1)=ind_aux;
                                inds(2)=ind_aux+1;
                                vals(1)=1-((ap_aux-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                                vals(2)=1-vals(1);
                                
                            end
                            
                            cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    cont=cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)*(vals(1)*V_ss(j+1,inds(1),etap,educ,married,kidsp)+vals(2)*V_ss(j+1,inds(2),etap,educ,married,kidsp));
                                end
                            end
                            
                            c_aux=consumption(j,a,eta,educ,married,kids,ap_aux,xi,b);
                            %                          c_aux=(1+r)*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc))-ap_aux;
                            
                            ap_VFI(j,a,eta,educ,married,kids)=ap_aux;
                            cons_VFI(j,a,eta,educ,married,kids)=c_aux;
                            
                            V_VFI(j,a,eta,educ,married,kids)=utility(c_aux,married,kids)+beta*psi(j)*cont;
                            
                            % Check end point of asset grid (ap=0)
                            c_aux3=consumption(j,a,eta,educ,married,kids,0,xi,b);
                            %                          c_aux3=(1+r)*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc));
                            
                            cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    cont=cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)*V_ss(j+1,1,etap,educ,married,kidsp);
                                end
                            end
                            V_aux3=utility(c_aux3,married,kids)+beta*psi(j)*cont;
                            
                            if V_aux3>V_VFI(j,a,eta,educ,married,kids)
                                ap_VFI(j,a,eta,educ,married,kids)=0;
                                cons_VFI(j,a,eta,educ,married,kids)=c_aux3;
                                
                                V_VFI(j,a,eta,educ,married,kids)=V_aux3;
                            end
                            
                            if cons_VFI(j,a,eta,educ,married,kids)<=0
                                disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                error('Non-positive consumption')
                            end
                            
                        elseif j==n_jgrid
                            
                            %                          inc=r*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ);
                            %                          spouse_inc=spousal_income(j,educ,kids,( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi)),SS(j,educ));
                            
                            ap_VFI(j,a,eta,educ,married,kids)=0;
                            cons_VFI(j,a,eta,educ,married,kids)=consumption(j,a,eta,educ,married,kids,ap_VFI(j,a,eta,educ,married,kids),xi,b);
                            %                          cons_VFI(j,a,eta,educ,married,kids)=(1+r)*agrid(a)+( epsilon(j,educ)*theta*exp(eta_grid(eta)) )*(xi+b*(1-xi))+SS(j,educ)+(married-1)*spouse_inc-max(0,Tax(inc,(married-1)*spouse_inc));
                            
                            if cons_VFI(j,a,eta,educ,married,kids)<=0
                                disp([j,a,eta,educ,married,kids,cons_VFI(j,a,eta,educ,married,kids)])
                                error('Non-positive consumption')
                            end
                            
                            V_VFI(j,a,eta,educ,married,kids)=utility(cons_VFI(j,a,eta,educ,married,kids),married,kids);
                            
                        end
                        
                    end
                end
            end
        end
    end
    
    if (bl_print_vfi)
        disp(strcat(['SNW_VFI_MAIN: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end
    
end

%% Timing and Profiling End
if (bl_timer)
    toc;    
    st_complete_vfi = strjoin(...
        ["Completed SNW_VFI_MAIN", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_vfi);
end

end
