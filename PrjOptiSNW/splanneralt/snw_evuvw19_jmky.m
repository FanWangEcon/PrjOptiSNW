%% SNW_EVUVW19_JMKY Expected 2019 Planner Value given Age, Y, Marry, Kids
%    Given 2019 JAEEMK, average over income bins. Income level is a
%    function of savings and shock in 2020. The planner looks at age, y,
%    marry and kids count in 2019, given expectation of shock transition
%    and optimal savings choices, and unemployment probability, solve for
%    the planner C and V expected values.
%
%    This function is also used by the Bush Check problem, no separate
%    tester is needed because the inputs are not 2019 specific. The
%    function works identically for our Bush stimulus and the Trump/Biden
%    stimulus problems.
%
%    This is for one check, given the check we are solving for EV19_JAEEMK.
%
%    [EV19_JMKY, EC19_JMKY] = SNW_EVUVW19_JMKY(WELF_CHECKS, ST_SOLU_TYPE,
%    MP_PARAMS, MP_CONTROLS, V_SS, V_UNEMP) provide V_SS and V_UNEMP solved
%    out elsewhere, and only get two outputs out.
%
%    See also SNW_EVUVW19_JMKY_MASS, SNW_EVUVW19_JMKY_ALLCHECKS,
%    SNW_EVUVW19_JAEEMK, SNW_EVUVW20_JAEEMK, SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jmky(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==7)
        [mp_params, mp_controls, ...
            ev19_jaeemk, ec19_jaeemk, ...
            Phi_true, Phi_true_jmky, ...
            inc_tot_ygroup_grid] = varargin{:};
    else
        error('Need to provide 7 parameter inputs');
    end

else
    close all;

%     st_solu_type = 'matlab_minimizer';
    st_solu_type = 'bisec_vec';
%     st_solu_type = 'grid_search';

    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
%     mp_params = snw_mp_param('default_dense');
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

    % Solve planner's problem
    n_incgrid=201; % Number of income groups
    n_incgrid_aux=round(0.75*n_incgrid);
    inc_grid1=linspace(0,4,n_incgrid_aux)'; % 4 refers to 4*58056=232224 dollars in 2012USD
    inc_grid=[inc_grid1;linspace(4+((7-4)/(n_incgrid-n_incgrid_aux)),7,n_incgrid-n_incgrid_aux)']; % 7 refers to 7*58056=406392 dollars in 2012USD

    mp_params('n_incgrid') = n_incgrid;
    mp_params('inc_grid') = inc_grid;

    % Solve for Unemployment Values
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_precompute') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jaeemk') = false;
    mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jmky_mass') = false;
    mp_controls('bl_print_evuvw19_jmky_mass_verbose') = false;
    mp_controls('bl_print_a4chk_verbose') = false;

    % Solve the Model to get V working and unemployed
    [V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [V_unemp,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);

    % Get Matrixes
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid',...
        'inc_tot_ygroup_grid', ...
        'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss'};
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);

    % Solve evuvw19 Given Current Check
    [ev19_jaeemk, ec19_jaeemk] = snw_evuvw19_jaeemk(...
        welf_checks, ...
        st_solu_type, mp_params, mp_controls, ...
        V_ss, cons_ss, ...
        V_unemp, cons_unemp, ...
        mp_precompute_res);

    % Get Mass
    params_group = values(mp_precompute_res, {'inc_tot_ygroup_grid'});
    [inc_tot_ygroup_grid] = params_group{:};
    [Phi_true_jmky] = snw_evuvw19_jmky_mass(...
        mp_params, mp_controls, Phi_true, inc_tot_ygroup_grid);

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
params_group = values(mp_controls, {'bl_print_evuvw19_jmky', 'bl_print_evuvw19_jmky_verbose'});
[bl_print_evuvw19_jmky, bl_print_evuvw19_jmky_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tm_start = tic;
end

%% Solve planner's problem

ev19_jmky=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);
ec19_jmky=zeros(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

if (bl_print_evuvw19_jmky)
    disp(strcat(['SNW_EVUVW19_JMKY Start']));
end

for j=1:(n_jgrid-1) % Age

    if (bl_print_evuvw19_jmky) tm_start_age = tic; end

    for married=1:n_marriedgrid % Marital status
        for kids=1:n_kidsgrid % Number of kids

            % income groups that are reached by this age/marry/kids group
            inc_tot_ygroup_jmk = inc_tot_ygroup_grid(j,:,:,:,married,kids);
            inc_tot_ygroup_jmk_unique = unique(inc_tot_ygroup_jmk(:))';

            % loop over income groups that are specific to this amk
            for inc_group=inc_tot_ygroup_jmk_unique
                % zero wage on in the income group, but no mass
                if(inc_group~=0)

                    % Get a, eta, educ nonzero index:
                    fl_Phi_mass_yjmk = Phi_true_jmky(j,married,kids,inc_group);
                    inc_tot_ytropu_jmky_ingroup = (inc_tot_ygroup_jmk == inc_group);
                    [~, a_idx, eta_idx, educ_idx] = ind2sub(size(inc_tot_ygroup_jmk), find(inc_tot_ytropu_jmky_ingroup));

                    % Within Y group weighted average
                    fl_ev19_jmky = 0;
                    fl_ec19_jmky = 0;

                    % Loop over valid points
                    for it_ctr=1:length(a_idx)
                        % Get Index
                        a = a_idx(it_ctr);
                        eta = eta_idx(it_ctr);
                        educ = educ_idx(it_ctr);

                        % Conditional Probablity Weights
                        condi_mass = Phi_true(j,a,eta,educ,married,kids)/fl_Phi_mass_yjmk;

                        % Weighted Average
                        fl_ev19_jmky = fl_ev19_jmky + condi_mass*ev19_jaeemk(j,a,eta,educ,married,kids);
                        fl_ec19_jmky = fl_ec19_jmky + condi_mass*ec19_jaeemk(j,a,eta,educ,married,kids);
                    end

                    % Store
                    ev19_jmky(j,married,kids,inc_group) = fl_ev19_jmky;
                    ec19_jmky(j,married,kids,inc_group) = fl_ec19_jmky;

                end
            end
        end
    end

    if (bl_print_evuvw19_jmky)
        tm_end_age = toc(tm_start_age);
        disp(strcat(['SNW_EVUVW19_JMKY: Finished Age Group:' ...
            num2str(j) ' of ' num2str(n_jgrid-1) ...
            ', time-this-age:' num2str(tm_end_age)]));
    end

end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_EVUVW19_JMKY", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete_vu_vw_checks);
end

%% Print
if (bl_print_evuvw19_jmky_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('ev19_jmky') = ev19_jmky;
    mp_outcomes('ec19_jmky') = ec19_jmky;
    mp_outcomes('ev19_jaeemk') = ev19_jaeemk;
    mp_outcomes('ec19_jaeemk') = ec19_jaeemk;
    mp_outcomes('Phi_true') = Phi_true;
    mp_outcomes('Phi_true_jmky') = Phi_true_jmky;
    ff_container_map_display(mp_outcomes);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev19_jmky;
    elseif (it_k==2)
        ob_out_cur = ec19_jmky;
    end
    varargout{it_k} = ob_out_cur;
end

end
