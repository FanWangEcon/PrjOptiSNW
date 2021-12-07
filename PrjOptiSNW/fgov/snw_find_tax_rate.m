%% SNW_FIND_TAX_RATE Solves for government budget clearning tax rate
%    Shift the tax policy variable to pay for various government outlays.
%    Solves for this problem in 2020 context given total COVID checks
%    outlays.
%
%    This function takes the as given distributional results as well as
%    policy functions. Note that we are, for the COVID welfare checks,
%    imposing a one-shot/one-period tax to pay for it as a boundary point
%    test-case. Given the distribution of the state-space that is
%    determined pre-covid, and the realization of the COVID shock, each
%    individual in the economy has some realized levels of resources. These
%    resources are invariant with respect to the one period government tax
%    policy to pay for COVID checks. The government taxes these resources
%    sufficiently to pay for COVID checks.
%
%    Given the checks that the government hands out and the taxes imposed,
%    individual resources post-tax are different in 2020. Households'
%    savings decisions in 2020 vary with taxes and checks. However, the
%    policy function post 2020 shifts back to thte previous non-COVID
%    world's policy function because the COVID shock is an one period
%    surprise shock.
%
%    * Covid_checks_per_capita: based on the actual allocation of covid
%    checks using actual policy criteria
%    * bl_adjust_a0: adjust a0, max tax bound, rather than a2.
%    * bl_load_existing: if to load existing saved vfi results distribution
%    and policy functions. this will save some time. the file is saved to
%    the output folder, and has only one name *snw_find_tax_rate*
%
%    The planner takes into consideration the responses of the households
%    to the COVID checks/taxes/shocks in building expected
%    values/consumptions for each household.
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] =
%    SNW_FIND_TAX_RATE()
%
%    See also SNW_DS_MAIN, SNW_DS_GRID_SEARCH
%

%%
function a2=snw_find_tax_rate(varargin)

%% Default and Parse Inputs

if (~isempty(varargin))
    
    bl_load_existing = false;
    
    if (length(varargin)==3)
        [mp_params, mp_controls, Covid_checks_per_capita] = varargin{:};
    elseif (length(varargin)==4)
        [mp_params, mp_controls, Covid_checks_per_capita, bl_adjust_a0] = varargin{:};
    elseif (length(varargin)==5)
        [mp_params, mp_controls, Covid_checks_per_capita, bl_adjust_a0, bl_load_existing] = varargin{:};
    else
        error('Need to provide 3 parameter inputs');
    end
    
else
    
    clc;
    clear all;
    %   mp_params = snw_mp_param('default_tiny');
    %   mp_params = snw_mp_param('default_small');
    mp_params = snw_mp_param('default_dense');
    %   mp_params = snw_mp_param('default_docdense');
    %   mp_params = snw_mp_param('default_base');
    %   mp_params = snw_mp_param('default_moredense', false, 'tauchen', true, 8, 8);
    %   mp_params = snw_mp_param('default_moredense_a65zh81zs5_e2m2', false, 'tauchen', true, 8, 8);
    mp_params('beta') = 0.95;
    xi=0.651; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    mp_params('xi') = xi;
    mp_params('b') = b;
    
    mp_controls = snw_mp_control('default_test');
    
    Covid_checks_per_capita = 18.7255856*100/62502;
    %   Covid_checks_per_capita = 0;
    
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_find_tax_rate') = true;
    mp_controls('bl_print_find_tax_rate_verbose') = true;
    
    bl_adjust_a0 = false;
    bl_load_existing = false;
    
end

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r', 'g_n', 'g_cons', 'a2','jret'});
[theta, r, g_n, g_cons, a2,jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp', 'n_incgrid', 'inc_grid'});
[pi_unemp, n_incgrid, inc_grid] = params_group{:};

params_group = values(mp_params, {'n_welfchecksgrid', 'xi','b'});
[n_welfchecksgrid, xi, b] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, ...
    {'bl_print_find_tax_rate', 'bl_print_find_tax_rate_verbose'});
[bl_print_find_tax_rate, bl_print_find_tax_rate_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% Check if Already saved DS and VFI
mp_paths = snw_mp_path('fan');
spn_mat_path = [mp_paths('spt_simu_outputs'), 'snw_find_tax_rate.mat'];

if (isfile(spn_mat_path) && bl_load_existing)
    load(spn_mat_path, 'Phi_true', 'A_agg', 'A_per_capita', ...
        'Y_inc_agg', 'Y_inc_agg_per_capita_1',...
        'Gov_cons', 'Gov_cons_per_capita',...
        'Y_inc_agg', 'cutoffs');
else
    
    %% A. Solve VFI non COVID world
    % 2. Solve VFI and Distributon
    % Solve the Model to get V working and unemployed
    [~,ap_ss,cons_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    
    %% B1. Solve Dist non COVID world
    % units translations. There are two units. One unit of when Phi_true sums
    % to 1, meaning is is a probability density. Another case where Phi_ture
    % are population measures with measure(age18)=1. Y_inc_agg from the
    % distribution function is in population measure units, code that follow
    % from here are also in population measure units. Covid_checks, however, is
    % in per-capita units, need to be multiplyed by sum of population measure.
    
    [Phi_true,Phi_adj,A_agg,Y_inc_agg] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss);
    
    % GDP per capita, the two numbers below should give the same answer.
    Y_inc_agg_per_capita_1 = Y_inc_agg/sum(Phi_true, 'all');
    Y_inc_agg_per_capita_2 = Y_inc_agg/sum(mp_params('Pop'));
    
    % Pre-covid steady State GOV Consumption
    Gov_cons = g_cons*Y_inc_agg;
    Gov_cons_per_capita = Gov_cons/sum(Phi_true, 'all');
    
    % Savings per capita
    A_per_capita = A_agg/sum(Phi_true, 'all');
    
    %% B2. Cutoffs
    cutoffs = snw_wage_cutoffs(Phi_true, theta, epsilon, eta_H_grid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid, jret);
    
    % save to mat
    save(spn_mat_path, 'Phi_true', 'A_agg', 'A_per_capita', ...
        'Y_inc_agg', 'Y_inc_agg_per_capita_1',...
        'Gov_cons', 'Gov_cons_per_capita',...
        'Y_inc_agg', 'cutoffs');
end

% Covid Check higher scale
Covid_checks = Covid_checks_per_capita*sum(Phi_true, 'all');
Covid_checks_share_of_GDP = Covid_checks_per_capita/Y_inc_agg_per_capita_1;

disp(['Y_inc_agg=' num2str(Y_inc_agg)]);
disp(['A_agg=' num2str(A_agg)]);
disp(['Y_inc_agg_per_capita_1=' num2str(Y_inc_agg_per_capita_1)]);
disp(['A_per_capita=' num2str(A_per_capita)]);
disp(['Gov_cons_per_capita=' num2str(Gov_cons_per_capita)]);
disp(['Covid_checks_share_of_GDP=' num2str(Covid_checks_share_of_GDP)]);

%% C. Aggregate Statistics Etc
% Aggregate variables
Y_inc_agg_COVID=0;
SS_spend=0;
UI_benefits=0;

for j=1:n_jgrid
    
    % Age Timer
    if (bl_print_find_tax_rate) tm_ds_age = tic; end
    
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid
                        
                        wages=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
                        if wages<=cutoffs(1)
                            wage_ind=1;
                        elseif wages>cutoffs(1) && wages<=cutoffs(2)
                            wage_ind=2;
                        elseif wages>cutoffs(2) && wages<=cutoffs(3)
                            wage_ind=3;
                        elseif wages>cutoffs(3) && wages<=cutoffs(4)
                            wage_ind=4;
                        elseif wages>cutoffs(4)
                            wage_ind=5;
                        end
                        
                        %                         [~,earn]=individual_income(j,a,eta,educ); % What individual earnings are before we account for the drop in earnings due to unemployment
                        [~,earn]=snw_hh_individual_income(j,a,eta,educ,...
                            theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option);
                        %                         spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ)); % What average spousal earnings are before we account for the drop in earnings due to unemployment
                        spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                        
                        inc_aux=pi_unemp(j,wage_ind)*epsilon(j,educ)*theta*exp(eta_H_grid(eta))*xi+(1-pi_unemp(j,wage_ind))*epsilon(j,educ)*theta*exp(eta_H_grid(eta))+r*(agrid(a)+Bequests*(bequests_option-1)); % Income (excluding Social Security benefits) after accounting for potential earnings drop in case of unemployment
                        
                        Y_inc_agg_COVID=Y_inc_agg_COVID+Phi_true(j,a,eta,educ,married,kids)*( inc_aux+pi_unemp(j,wage_ind)*(married-1)*spouse_inc*exp(eta_S_grid(eta))*xi+(1-pi_unemp(j,wage_ind))*(married-1)*spouse_inc*exp(eta_S_grid(eta)) ); % Aggregate income (labor earnings, spousal income, interest earnings)
                        
                        SS_spend=SS_spend+sum(Phi_true(j,a,eta,educ,married,kids)*SS(j,educ)); % Total spending on Social Security
                        
                        % 2021 12/05, the first commented out line below
                        % includes spousal income, but our model only has 
%                         UI_benefits=UI_benefits+sum(Phi_true(j,a,eta,educ,married,kids))*pi_unemp(j,wage_ind)*( epsilon(j,educ)*theta*exp(eta_H_grid(eta))+(married-1)*spouse_inc*exp(eta_S_grid(eta)) )*b*(1-xi); % Total spending on unemployment insurance benefits
                        UI_benefits=UI_benefits+sum(Phi_true(j,a,eta,educ,married,kids))*pi_unemp(j,wage_ind)*(epsilon(j,educ)*theta*exp(eta_H_grid(eta)))*b*(1-xi); % Total spending on unemployment insurance benefits
                        
                    end
                end
            end
        end
    end
    
    if (bl_print_find_tax_rate)
        tm_ds_age_end = toc(tm_ds_age);
        disp(strcat(['SNW_FIND_TAX_RATE: Aggregation, Finished Age Group:' ...
            num2str(j) ' of ' num2str(n_jgrid) ...
            ', time-this-age:' num2str(tm_ds_age_end)]));
    end
    
end

%% D. Tax Adjustments to Pay for Government Bills
% Find value of a2 that balances government budget
tol = 10^-4;
tol_gap = 10^-4;
err=1;
err_last=99;
err_gap=1;
a2_guess_orig=mp_params('a2_covidyr_manna_heaven');
a2 = mp_params('a2_covidyr_manna_heaven');

% Income aggregation


% Tax parameters
a0_guess_orig=mp_params('a0');
a0 = mp_params('a0');

it=0;
while err>tol && err_gap>tol_gap
    
    it=it+1;
    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
        for a=1:n_agrid
            for eta=1:n_etagrid
                for educ=1:n_educgrid
                    for married=1:n_marriedgrid
                        for kids=1:n_kidsgrid
                            % [inc,earn]=individual_income(j,1:n_agrid,eta,educ,xi,b);
                            [inc,earn]=snw_hh_individual_income(j,a,eta,educ,...
                                theta, r, agrid, epsilon, eta_H_grid, SS, Bequests, bequests_option,...
                                xi,b);
                            
                            % spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ)); % What average spousal earnings are before we account for the drop in earnings due to unemployment
                            spouse_inc=snw_hh_spousal_income(j,educ,kids,earn,SS(j,educ), jret);
                            
                            % Tax_revenues_aux = Tax_revenues_aux+Phi_true(j,1:n_agrid,eta,educ,married,kids)*...
                            %                     max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi)))); % Tax revenues
                            Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,a,eta,educ,married,kids)*...
                                max(0,snw_tax_hh(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta))*(xi+b*(1-xi)), a2, a0)); % Tax revenues
                            
                        end
                    end
                end
            end
        end
    end
    
    fl_total_costs = UI_benefits+SS_spend+Gov_cons+Covid_checks;
    if (bl_adjust_a0)
        a0=a0*((fl_total_costs/Tax_revenues_aux)^0.75); % Find value of a0 that balances government budget
    else
        a2_last = a2;
        a2=a2*((fl_total_costs/Tax_revenues_aux)^0.75); % Find value of a2 that balances government budget
    end
    err=abs((Tax_revenues_aux/(fl_total_costs))-1);
    err_gap = abs(err-err_last);
    err_last = err;
    if (bl_print_find_tax_rate)
        st_tax_iter = strjoin(...
            ["SNW_FIND_TAX_RATE tax a2 or a0 adjustments", ...
            ['a2=' num2str(a2)] ...
            ['a0=' num2str(a0)] ...
            ['err=' num2str(err)] ...
            ['fl_total_costs=' num2str(fl_total_costs)] ...
            ['Tax_revenues_aux=' num2str(Tax_revenues_aux)] ...
            ], ";");
        disp(st_tax_iter);
    end
    
end


if (bl_print_find_tax_rate)
    
    disp(['--------------------------------']);
    disp(['--- SNW_FIND_TAX_RATE finished -']);
    disp(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    disp(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    disp(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    
    name='SNW_FIND_TAX_RATE: Number of a2-adjustments (for taxation) used to balance the government budget= ';
    name2=[name,num2str(it)];
    disp(name2);
    
    st_tax_iter = strjoin(...
        ["SNW_FIND_TAX_RATE info A", ...
        ['a2_new=' num2str(a2)] ...
        ['a2_guess_orig(mana-heaven)=' num2str(a2_guess_orig)] ...
        ['a0_new=' num2str(a0)] ...
        ['a0_guess_orig=' num2str(a0_guess_orig)] ...
        ], ";");
    disp(st_tax_iter);
    
    st_tax_iter = strjoin(...
        ["SNW_FIND_TAX_RATE info B", ...
        ['Y_inc_agg=' num2str(Y_inc_agg)] ...
        ['Y_inc_agg_COVID=' num2str(Y_inc_agg_COVID)] ...
        ['GPD_COVID_CHANGE=' num2str((Y_inc_agg-Y_inc_agg_COVID)/Y_inc_agg)] ...
        ['Y_inc_agg_per_capita=' num2str(Y_inc_agg/sum(Phi_true, 'all'))] ...
        ['Y_inc_agg_per_capita_COVID=' num2str(Y_inc_agg_COVID/sum(Phi_true, 'all'))] ...
        ], ";");
    disp(st_tax_iter);
    
    st_tax_iter = strjoin(...
        ["SNW_FIND_TAX_RATE info C", ...
        ['Covid_checks_per_capita=' num2str(Covid_checks_per_capita)] ...
        ['Covid_checks_share_of_GDP=' num2str(Covid_checks_share_of_GDP)] ...
        ], ";");
    disp(st_tax_iter);
    
    st_tax_iter = strjoin(...
        ["SNW_FIND_TAX_RATE info D", ...
        ['UI_benefits=' num2str(UI_benefits)] ...
        ['SS_spend=' num2str(SS_spend)] ...
        ['Gov_cons=' num2str(Gov_cons)] ...
        ['Covid_checks=' num2str(Covid_checks)] ...
        ['Tax_revenues_aux=' num2str(Tax_revenues_aux)] ...
        ], ";");
    disp(st_tax_iter);
    
    disp(['xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx']);
    disp(['--------------------------------']);
    
end

end
