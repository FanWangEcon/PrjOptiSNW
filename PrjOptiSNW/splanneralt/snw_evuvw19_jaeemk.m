%% SNW_EVUVW19_JAEEMK Solves for EVUVW from 2019 Conditional JAEEMK
%    Given 2020 JAEEMK, what is the expectated value for the planner given
%    2020 JAEEMK and transition between 2019 to 2020 JAEEMK. For one Check,
%    the check number set by WELF_CHECKS.
%
%    [EV19_JAEEMK, EC19_JAEEMK, EV20_JAEEMK, EC20_JAEEMK] =
%    SNW_EVUVW19_JAEEMK(WELF_CHECKS, ST_SOLU_TYPE, MP_PARAMS, MP_CONTROLS,
%    V_SS, CONS_SS, V_UNEMP, CONS_UNEMP, MP_PRECOMPUTE_RES) provide V_SS
%    and V_UNEMP solved out elsewhere, and only get two outputs out.
%
%    See also SNW_EVUVW19_JYMK, SNW_EVUVW19_JAEEMK, SNW_EVUVW20_JAEEMK,
%    SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jaeemk(varargin)

%% Default and Parse
if (~isempty(varargin))

    if (length(varargin)==9)
        [welf_checks, st_solu_type, mp_params, mp_controls, ...
            V_ss, cons_ss, ...
            V_unemp, cons_unemp, ...
            mp_precompute_res] = varargin{:};
    else
        error('Need to provide 9 parameter inputs');
    end

else
    clc;
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

    % Solve for Unemployment Values
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_precompute') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_a4chk_verbose') = false;

    % Solve the Model to get V working and unemployed
    [V_ss,ap_ss,cons_ss,mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    [V_unemp,~,cons_unemp,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
    [Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);

    % Get Matrixes
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid',...
        'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss'};
    mp_controls('bl_print_precompute_verbose') = false;
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);

end

%% Parse Pre-Computes
params_group = values(mp_precompute_res, ...
    {'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss'});
[ap_idx_lower_ss, ap_idx_higher_ss, ap_idx_lower_weight_ss] = params_group{:};

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, cl_mt_pi_jem_kidseta, psi] = params_group{:};
params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_evuvw19_jaeemk', 'bl_print_evuvw19_jaeemk_verbose'});
[bl_print_evuvw19_jaeemk, bl_print_evuvw19_jaeemk_verbose] = params_group{:};

%% Solve evuvw19 Given Current Check
[ev20_jaeemk, ec20_jaeemk] = snw_evuvw20_jaeemk(...
    welf_checks, ...
    st_solu_type, mp_params, mp_controls, ...
    V_ss, cons_ss, ...
    V_unemp, cons_unemp, ...
    mp_precompute_res);

%% Timing and Profiling Start

if (bl_timer)
    tm_start = tic;
end

%% Planner val2019(J,A,E,E,M,K) = E(val2020(J+1,A',E',E,M,K')) given trans.
ev19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
ec19_jaeemk=NaN(n_jgrid-1,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

bl_vectorize = false;
bl_vectorize_a_only = false;

% Solve for value function and policy functions by means of backwards induction
for j=1:n_jgrid-1 % Age
    for educ=1:n_educgrid % Educational level
        for married=1:n_marriedgrid % Marital status
            for kids=1:n_kidsgrid % Number of kids

                % Age, married and kids specific transition
                mt_pi_jem_kidseta = kron(pi_kids(:,:,j,educ,married), pi_eta);

                % loop of current uncertainty
                for eta=1:n_etagrid % Productivity

                    if (bl_vectorize)

                        if (bl_vectorize_a_only)

                            % Solve for value function and policy functions by means of backwards induction

                            % Choice Index and Weights
                            inds_l=ap_idx_lower_ss(j,:,eta,educ,married,kids);
                            inds_h=ap_idx_higher_ss(j,:,eta,educ,married,kids);
                            % Linear interpolation
                            vals_l=ap_idx_lower_weight_ss(j,:,eta,educ,married,kids);
                            vals_h=1-vals_l;

                            % Continuation Value V = Planner Value 2019
                            val_cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    val_cont=val_cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)...
                                        *(vals_l.*ev20_jaeemk(j+1,inds_l,etap,educ,married,kidsp)...
                                        +vals_h.*ev20_jaeemk(j+1,inds_h,etap,educ,married,kidsp));
                                end
                            end

                            % Continuation Value C = Planner Value 2019
                            cons_cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    cons_cont=cons_cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)...
                                        *(vals_l.*ec20_jaeemk(j+1,inds_l,etap,educ,married,kidsp)...
                                        +vals_h.*ec20_jaeemk(j+1,inds_h,etap,educ,married,kidsp));
                                end
                            end

                            ev19_jaeemk(j,:,eta,educ,married,kids) = val_cont;
                            ec19_jaeemk(j,:,eta,educ,married,kids) = cons_cont;

                        else
                            % Get P(s'|s) conditional on current s
                            ar_pi_eta_kids = mt_pi_jem_kidseta((n_etagrid*(kids-1)+eta),:);

                            % Choice Index and Weights
                            mn_idx_lower=ap_idx_lower_ss(j,:,eta,educ,married,kids);
                            mn_idx_higher=ap_idx_higher_ss(j,:,eta,educ,married,kids);

                            % Linear interpolation
                            mn_consprime_lower=ap_idx_lower_weight_ss(j,:,eta,educ,married,kids);
                            mn_consprime_higher=1-mn_consprime_lower;

                            % flatten
                            ar_idx_lower = mn_idx_lower(:);
                            ar_idx_higher = mn_idx_higher(:);
                            ar_val_lower=mn_consprime_lower(:);
                            ar_val_higher=mn_consprime_higher(:);

                            % Continuation Value V = Planner Value 2019
                            mn_consprime_lower = ev20_jaeemk(j+1,ar_idx_lower,:,educ,married,:);
                            mn_consprime_higher = ev20_jaeemk(j+1,ar_idx_higher,:,educ,married,:);
                            mn_consprime_weighted_average = (ar_val_lower'.*mn_consprime_lower + ar_val_higher'.*mn_consprime_higher);
                            % Reshape so Shocks are Columns and a are rows
                            mt_cons_weighted_average = reshape(mn_consprime_weighted_average, n_agrid, []);
                            % EV along asset states
                            ar_ev_cont=mt_cons_weighted_average*ar_pi_eta_kids';

                            % Continuation Value C = Planner C 2019
                            mn_consprime_lower = ec20_jaeemk(j+1,ar_idx_lower,:,educ,married,:);
                            mn_consprime_higher = ec20_jaeemk(j+1,ar_idx_higher,:,educ,married,:);
                            mn_consprime_weighted_average = (ar_val_lower'.*mn_consprime_lower + ar_val_higher'.*mn_consprime_higher);
                            % Reshape so Shocks are Columns and a are rows
                            mt_cons_weighted_average = reshape(mn_consprime_weighted_average, n_agrid, []);
                            % EV along asset states
                            ar_cons_cont=mt_cons_weighted_average*ar_pi_eta_kids';

                            % Store
                            ev19_jaeemk(j,:,eta,educ,married,kids) = ar_ev_cont;
                            ec19_jaeemk(j,:,eta,educ,married,kids) = ar_cons_cont;
                        end

                    else

                        % Solve for value function and policy functions by means of backwards induction
                        ar_ev_cont2 = zeros([n_agrid,1]);
                        ar_ec_cont2 = zeros([n_agrid,1]);
                        for a=1:n_agrid

                            % Choice Index and Weights
                            inds(1)=ap_idx_lower_ss(j,a,eta,educ,married,kids);
                            inds(2)=ap_idx_higher_ss(j,a,eta,educ,married,kids);
                            % Linear interpolation
                            vals(1)=ap_idx_lower_weight_ss(j,a,eta,educ,married,kids);
                            vals(2)=1-vals(1);

                            % Continuation Value V = Planner Value 2019
                            val_cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    val_cont=val_cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)...
                                        *(vals(1)*ev20_jaeemk(j+1,inds(1),etap,educ,married,kidsp)...
                                        +vals(2)*ev20_jaeemk(j+1,inds(2),etap,educ,married,kidsp));
                                end
                            end

                            % Continuation Value C = Planner Value 2019
                            cons_cont=0;
                            for etap=1:n_etagrid
                                for kidsp=1:n_kidsgrid
                                    cons_cont=cons_cont+pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married)...
                                        *(vals(1)*ec20_jaeemk(j+1,inds(1),etap,educ,married,kidsp)...
                                        +vals(2)*ec20_jaeemk(j+1,inds(2),etap,educ,married,kidsp));
                                end
                            end

                            % Store
                            ar_ev_cont2(a) = val_cont;
                            ar_ec_cont2(a) = cons_cont;
                        end

                        ev19_jaeemk(j,:,eta,educ,married,kids) = ar_ev_cont2;
                        ec19_jaeemk(j,:,eta,educ,married,kids) = ar_ec_cont2;
                    end

                end
            end
        end
    end
end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_EVUVW19_JAEEMK", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Print
if (bl_print_evuvw19_jaeemk_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('ev19_jaeemk') = ev19_jaeemk;
    mp_outcomes('ec19_jaeemk') = ec19_jaeemk;
    mp_outcomes('ev20_jaeemk') = ev20_jaeemk;
    mp_outcomes('ec20_jaeemk') = ec20_jaeemk;
    ff_container_map_display(mp_outcomes, 9, 9);
%     disp(ec19_jaeemk);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev19_jaeemk;
    elseif (it_k==2)
        ob_out_cur = ec19_jaeemk;
    elseif (it_k==3)
        ob_out_cur = ev20_jaeemk;
    elseif (it_k==4)
        ob_out_cur = ec20_jaeemk;
    end
    varargout{it_k} = ob_out_cur;
end

end
