%% SNW_DS_AGGREGATION Given Distribution and Policy obtain Aggregates
%    This was previously a part of the SNW_DS_AGGREGATION and SNW_DS_MAIN
%    functions. This is separated out to allow for aggregation over
%    different education and discount factor types. 
%
%    The version inside SNW_DS_MAIN iterates over tax rates, this version
%    allows for controlling whether iteration should take place. Should not
%    if there is an outter loop for iteration over different types. 
%
%    * IT_TAX_ITER_MAX int maximum number of tax iteration. Should set to
%    0 if do not want to iterate to find tax rate.
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] =
%    SNW_DS_AGGREGATION() invoke model with externally set parameter map and
%    control map. Results outputed to a map containing various output
%    matrixes in mp_dsvfi_results, and also distributional matrixes.
%
%    See also SNW_DS_MAIN, SNW_DS_GRID_SEARCH
%

%%
function varargout=snw_ds_aggregation(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))

    it_tax_iter_max = 1;
    fl_rev_equal_tax_tol = 1e-4;
    
    if (length(varargin)==5)
        [mp_params, mp_controls, ap_ss, cons_ss, Phi_true] = varargin{:};
    elseif (length(varargin)==6)
        [mp_params, mp_controls, ap_ss, cons_ss, Phi_true, it_tax_iter_max] = varargin{:};
    elseif (length(varargin)==7)
        [mp_params, mp_controls, ap_ss, cons_ss, Phi_true, ...
            it_tax_iter_max, fl_rev_equal_tax_tol] = varargin{:};
    else
        error('Need to provide 2/4 parameter inputs');
    end

else

    clc;
    clear all;
%     mp_params = snw_mp_param('default_docdense');
%     mp_params = snw_mp_param('default_small');
%     mp_params = snw_mp_param('default_tiny', false, 'tauchen', true, 8, 8);    
%     mp_params = snw_mp_param('default_tiny');
    
    mp_more_inputs = containers.Map('KeyType','char', 'ValueType','any');
%     mp_more_inputs('st_edu_simu_type') = 'both';
%     mp_params = snw_mp_param('default_tiny', false, 'tauchen', false, 8, 8, mp_more_inputs);
    mp_more_inputs('st_edu_simu_type') = 'low';
    mp_params = snw_mp_param('default_tiny_e1l', false, 'tauchen', false, 8, 8, mp_more_inputs);
%     mp_more_inputs('st_edu_simu_type') = 'high';
%     mp_params = snw_mp_param('default_tiny_e2h', false, 'tauchen', false, 8, 8, mp_more_inputs);

    mp_params('a2') = 1.2445;
    mp_params('theta') = 0.565228521783443;
    mp_params('beta') = 0.971162552785405;

    mp_controls = snw_mp_control('default_test');
    [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);

    % distribution
    [Phi_true] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss);
    
    % max tax iteration 
    it_tax_iter_max = 0;
    fl_rev_equal_tax_tol = 1e-4;
    
end

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r', 'g_n', 'g_cons', 'a2','jret'});
[theta, r, g_n, g_cons, a2,jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'bl_store_shock_trans', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, bl_store_shock_trans, cl_mt_pi_jem_kidseta, psi] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, ...
    {'stat_distr_eta', 'stat_distr_educ', 'stat_distr_married', 'stat_distr_kids', 'Pop'});
[stat_distr_eta, stat_distr_educ, stat_distr_married, stat_distr_kids, Pop] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_ds_aggregation', 'bl_print_ds_aggregation_verbose'});
[bl_print_ds_aggregation, bl_print_ds_aggregation_verbose] = params_group{:};

%% Timer
if (bl_timer)
    tm_start = tic;
end

%% Explicit Looped Aggregate variables Calculation
A_agg=0;
Aprime_agg=0;
C_agg=0;
Y_inc_agg=0;
Tax_revenues=0;
SS_spend=0;
Bequests_aux=0;

for j=1:n_jgrid
    for eta=1:n_etagrid
        for educ =1:n_educgrid
            for married=1:n_marriedgrid
                for kids=1:n_kidsgrid

                    A_agg=A_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*agrid(1:n_agrid); % Aggregate wealth
                    Aprime_agg=Aprime_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*ap_ss(j,1:n_agrid,eta,educ,married,kids)'; % Aggregate saving
                    C_agg=C_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*cons_ss(j,1:n_agrid,eta,educ,married,kids)'; % Aggregate consumption

                    % [inc,earn]=individual_income(j,1:n_agrid,eta,educ);
                    [inc,earn]=snw_hh_individual_income(j,1:n_agrid,eta,educ,...
                        theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                    % spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                    spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);

                    inc_aux=r*(agrid(1:n_agrid)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta)); % Income (excluding Social Security benefits)

                    Y_inc_agg=Y_inc_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*( inc_aux+(married-1)*spouse_inc*exp(eta_S_grid(eta)) ); % Aggregate income (labor earnings, spousal income, interest earnings)

                    Tax_revenues=Tax_revenues+Phi_true(j,1:n_agrid,eta,educ,married,kids)*max(0,snw_tax_hh(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)),a2)); % Tax revenues

                    SS_spend=SS_spend+sum(Phi_true(j,1:n_agrid,eta,educ,married,kids)*SS(j,educ)); % Total spending on Social Security

                    Bequests_aux=Bequests_aux+Phi_true(j,1:n_agrid,eta,educ,married,kids)*ap_ss(j,1:n_agrid,eta,educ,married,kids)'*(1-psi(j)); % Accidental Bequests*(bequests_option-1)

                end
            end
        end
    end
end

% Alternative 1: Suppose accidental bequests go to the government
% Alternative 2: Allocate accidental bequests uniformly across the population
if bequests_option==1
    if throw_in_ocean==0
        if (bl_print_ds_aggregation)
            disp('SNW_DS_MAIN: Accidental bequests are part of government revenues');
        end
        Tax_revenues=Tax_revenues+Bequests_aux*(1+r);
    elseif throw_in_ocean==1
        if (bl_print_ds_aggregation)
            disp('SNW_DS_MAIN: Accidental bequests are thrown in the ocean');
        end
    end
elseif bequests_option==2
    if (bl_print_ds_aggregation)
        disp('SNW_DS_MAIN: Accidental bequests are uniformly distributed across the population');
    end
    Bequests=Bequests_aux/(sum(Pop)*(1+g_n));
end

% Update guess for a2 (determines average level of income taxation)
% Assuming government balances its budget period-by-period

% tol=10^-4; %10^-3; %5*10^-4;
fl_err_rev_iter_gap=abs((Tax_revenues/(SS_spend+g_cons*Y_inc_agg))-1);

a2_init = a2;
a2_update=a2;
ar_a2_store = NaN([it_tax_iter_max+1,1]);
ar_a2_store(1,1) = a2_init;

it_tax=0;

% turn off tolerance, which is controlled outside
while fl_err_rev_iter_gap>fl_rev_equal_tax_tol && it_tax <= it_tax_iter_max
% while it_tax <= it_tax_iter_max

    it_tax=it_tax+1;

    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid

                        % [inc,earn]=individual_income(j,1:n_agrid,eta,educ);
                        [inc,earn]=snw_hh_individual_income(j,1:n_agrid,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                        % spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,1:n_agrid,eta,educ,married,kids)*max(0,snw_tax_hh(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)),a2)); % Tax revenues

                    end
                end
            end
        end
    end

    if bequests_option==1
        if throw_in_ocean==0
            Tax_revenues_aux=Tax_revenues_aux+Bequests_aux*(1+r);
        end
    end

    a2=a2*(((SS_spend+g_cons*Y_inc_agg)/Tax_revenues_aux)^0.75); % Find value of a2 that balances government budget

    fl_err_rev_iter_gap=abs((Tax_revenues_aux/(SS_spend+g_cons*Y_inc_agg))-1);

    if (bl_print_ds_aggregation)
        st_tax_iter = strjoin(...
            ["SNW_DS_AGGREGATION tax and spend", ...
            ['it=' num2str(it_tax)], ...
            ['err=' num2str(fl_err_rev_iter_gap)] ...
            ], ";");
        disp(st_tax_iter);
    end
    
    % Track/store
    ar_a2_store(it_tax+1,1) = a2;    
end

% keep non NaN
ar_a2_store = ar_a2_store(~isnan(ar_a2_store));

% display
if (bl_print_ds_aggregation)
    name='SNW_DS_AGGREGATION: Number of a2-adjustments (for taxation) used to balance the government budget= ';
    name2=[name,num2str(it_tax)];
    disp(name2);

    a2_update=[a2_update,a2];

    name='SNW_DS_AGGREGATION: Old and updated value of a2=';
    name2=[name,num2str(a2_update)];
    disp(name2);

    name='SNW_DS_AGGREGATION: Aggregates: Cons., Gov. cons., Save, Assets, Income, Bequests ';
    aggregates=[C_agg,g_cons*Y_inc_agg,A_agg,Aprime_agg,Y_inc_agg,Bequests_aux];
    name2=[name,num2str(aggregates)];
    disp(name2);

    name='SNW_DS_AGGREGATION: Resource constraint: C_t+A_{t+1}+G_t=A_t+Y_t ';
    name2=[name,num2str([C_agg+g_cons*Y_inc_agg+Aprime_agg,A_agg+Y_inc_agg])];
    disp(name2);
end

a2_final = a2;

%% Timing and Profiling End
if (bl_timer && bl_print_ds_aggregation)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_DS_AGGREGATION", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
end

%% Collect Results if Require Extra Outputs
mp_agg_results = containers.Map('KeyType', 'char', 'ValueType', 'any');    
mp_agg_results('A_agg') = A_agg;
mp_agg_results('Aprime_agg') = Aprime_agg;
mp_agg_results('Y_inc_agg') = Y_inc_agg;
mp_agg_results('C_agg') = C_agg;
mp_agg_results('Tax_revenues') = Tax_revenues;
mp_agg_results('SS_spend') = SS_spend;
mp_agg_results('Bequests_aux') = Bequests_aux;

mp_agg_results('A_agg_perhh') = A_agg/sum(Pop);
mp_agg_results('Aprime_agg_perhh') = Aprime_agg/sum(Pop);
mp_agg_results('Y_inc_agg_perhh') = Y_inc_agg/sum(Pop);

mp_agg_results('C_agg_perhh') = C_agg/sum(Pop);
mp_agg_results('Tax_revenues_perhh') = Tax_revenues/sum(Pop);
mp_agg_results('SS_spend_perhh') = SS_spend/sum(Pop);
mp_agg_results('Bequests_aux_perhh') = Bequests_aux/sum(Pop);

mp_agg_results('it_tax') = it_tax;
mp_agg_results('a2') = a2;
mp_agg_results('a2_init') = a2_init;
mp_agg_results('a2_final') = a2_final;
mp_agg_results('ar_a2_store') = ar_a2_store;

mp_agg_results('fl_err_rev_iter_gap') = fl_err_rev_iter_gap;

if (bl_print_ds_aggregation_verbose)
    ff_container_map_display(mp_agg_results, it_tax+1, 1);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = SS_spend;
    elseif (it_k==2)
        ob_out_cur = Y_inc_agg;
    elseif (it_k==3)
        ob_out_cur = Tax_revenues_aux;
    elseif (it_k==4)
        ob_out_cur = Bequests_aux;        
    elseif (it_k==5)
        ob_out_cur = it_tax;
    elseif (it_k==6)
        % update the tax parameter
        ob_out_cur = a2;
    elseif (it_k==7)
        ob_out_cur = mp_agg_results;
    end
    varargout{it_k} = ob_out_cur;
end

end
