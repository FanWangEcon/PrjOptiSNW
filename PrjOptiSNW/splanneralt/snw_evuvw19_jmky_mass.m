%% SNW_EVUVW19_JMKY_MASS Mass At each Age, Income, Marry, Kids Group
%    Mass at each Age, Income, Marriage Kids group.
%
%    [Phi_true_jmky] = SNW_EVUVW19_JMKY_MASS(MP_PARAMS, MP_CONTROLS,
%    PHI_TRUE, INC_TOT_YGROUP_GRID) calling the function retures the
%    population mass Phi_true array, and the pre-computed
%    INC_TOT_YGROUP_GRID.
%
%    See also SNW_EVUVW19_JYMK, SNW_EVUVW19_JAEEMK, SNW_EVUVW20_JAEEMK,
%    SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jmky_mass(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==4)
        [mp_params, mp_controls, Phi_true, inc_tot_ygroup_grid] = varargin{:};
    else
        error('Need to provide 4 parameter inputs');
    end
    
else
    close all;
    
% %     % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
%     mp_params = snw_mp_param('default_small');
    mp_params = snw_mp_param('default_dense');
    mp_controls = snw_mp_control('default_test');
    
    % Solve for Unemployment Values
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;    
    mp_controls('bl_print_precompute') = false;
    
    % Solve the Model to get V working and unemployed
    [~,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [Phi_true, Phi_adj] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
    
    % Get Matrixes
    cl_st_precompute_list = {'inc_tot_ygroup_grid'};
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);    
    params_group = values(mp_precompute_res, {'inc_tot_ygroup_grid'});
    [inc_tot_ygroup_grid] = params_group{:};
    
end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_jgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'n_incgrid'});
[n_incgrid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_evuvw19_jmky_mass', 'bl_print_evuvw19_jmky_mass_verbose'});
[bl_print_evuvw19_jmky_mass, bl_print_evuvw19_jmky_mass_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Solve planner's problem
Phi_true_jmky=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

if (bl_print_evuvw19_jmky_mass)
    disp(strcat(['SNW_EVUVW19_JMKY_MASS Start']));
end

for j=1:(n_jgrid-1) % Age    
    for married=1:n_marriedgrid % Marital status
        for kids=1:n_kidsgrid % Number of kids
            
            % income groups that are reached by this age/marry/kids group
            Phi_true_jmk = Phi_true(j,:,:,:,married,kids);
            inc_tot_ygroup_jmk = inc_tot_ygroup_grid(j,:,:,:,married,kids);
            inc_tot_ygroup_jmk_unique = unique(inc_tot_ygroup_jmk(:))';
            
            % loop over income groups that are specific to this amk
            for inc_group=inc_tot_ygroup_jmk_unique                
                % zero wage on in the income group, but no mass
                if(inc_group~=0)
                    % Mass at this amk point
                    Phi_true_yjmk_ingroup = (inc_tot_ygroup_jmk == inc_group);
                    fl_Phi_mass_yjmk = sum(Phi_true_jmk(Phi_true_yjmk_ingroup), 'all');
                    Phi_true_jmky(j,married,kids,inc_group) = fl_Phi_mass_yjmk;                    
                end
            end
            
        end
    end    
end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_EVUVW19_JMKY_MASS", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete_vu_vw_checks);
end

%% Print
if (bl_print_evuvw19_jmky_mass_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');    
    mp_outcomes('Phi_true') = Phi_true;
%     mp_outcomes('Phi_adj') = Phi_adj;
    mp_outcomes('Phi_true_jmky') = Phi_true_jmky;
    ff_container_map_display(mp_outcomes);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = Phi_true_jmky;
    end
    varargout{it_k} = ob_out_cur;
end

end
