%% SNW_EVUVW19_JMKY_ALLCHECKS 2019 Planner Values All Checks
%    Call the SNW_EVUVW19_JMKY function to find the Value at All Checks.
%    Solving for one check at a time and storing resulting jmky results one
%    check at a time. Uses pre-stored matrixes that are re-used across
%    check calculations
%
%    WARNING: DO NOT run other files that will shift
%    snw_hh_spousal_income.m if running this file, otherwise will over-ride
%    low vs high education income equation specifications. Go work on a
%    different computer for example, or some other task. multiple beta/edu
%    jobs can not be run jointly under current "dynamic" file generation
%    structure. 
%      
%    * ST_BIDEN_OR_TRUMP string with values trumpchk or bidenchk
%    * BL_LOAD_MAT boolean if true load saving mat file with VFI and
%    distributiona results if the file exists
%
%    [EV19_JMKY_ALLCHECKS, EC19_JMKY_ALLCHECKS, OUTPUT] =
%    SNW_EVUVW19_JMKY_ALLCHECKS(MP_PARAMS, MP_CONTROLS, ST_SOLU_TYPE,
%    BL_PARFOR, IT_WORKERS, BL_EXPORT, SNM_SUFFIX) provide V_SS and V_UNEMP
%    solved out elsewhere, and only get two outputs out. Solves for Trump
%    check. 
%
%    [EV19_JMKY_ALLCHECKS, EC19_JMKY_ALLCHECKS, OUTPUT] =
%    SNW_EVUVW19_JMKY_ALLCHECKS(MP_PARAMS, MP_CONTROLS, ST_BIDEN_OR_TRUMP,
%    ST_SOLU_TYPE, BL_PARFOR, IT_WORKERS, BL_EXPORT, SNM_SUFFIX) parameter
%    ST_BIDEN_OR_TRUMP is added to consider optionally either trump or
%    biden check. For biden check, will assume: (1) Trump check manna from
%    heaven (so same tax rate as under steady state; (2) Trump check fully
%    replaces income loss due to covid; (3) Double MIT shock. 
%
%    See also SNW_EVUVW19_JMKY, SNW_EVUVW19_JMKY_MASS, SNW_EVUVW19_JAEEMK_FOC,
%    SNW_EVUVW20_JAEEMK, SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jmky_allchecks(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    st_biden_or_trump = 'trumpchk';
    
    if (length(varargin)==8)
        % Original invokation structure, solves for optimal Trump check
        % policy
        [mp_params, mp_controls, st_solu_type, ...
            bl_parfor, it_workers, ...
            bl_export, bl_load_mat, snm_suffix] = varargin{:};
        
    elseif (length(varargin)==9)
        % Modified inputs to allow for solving for optimal Biden check
        % policy. The difference from prior is that the distributional code
        % will be based on optimal policy from 1 period MIT shock but also
        % given the Trump stimulus, assuming by default that there is manna
        % from heaven, and by default for year 1 of covid there is full
        % wage replacement. 
        
        % From Trump to Biden Checks For Biden policy, (2b) provides
        % mass/weights at state-space elements, and (3a) provides optimal
        % choices and corresponding MPCs at each check level. We have four
        % files with the MPC at each check for
        % income-bin/marital/kids-count (2020 info for 2021 policy), each
        % for a different beta/edu group:
        % 
        % * 1a, solve for "steady-state" policy/value function. (think
        % 2019, 2022, and after)
        % * 1b, solve for distributions given "steady-state" policy.value
        % functions.
        % * 2a. (first MIT shock), given "steady-state" value function as
        % continuation value, solve for Trump-stimulus policy/value
        % functions. (2020 choices)
        % * 2b. conditional (1b) distribution, one-period forward
        % post-trump-policy distribution with (2a) policy functions. (start
        % of 2021 distribution given 2020 choices)
        % * 3a. (second MIT shock), given "steady-state" value function as
        % continuation value, solve for Biden-stimulus policy/value
        % functions (2021 choices)
        % 
        % Following our benchmark in the prior draft and to simplify, when
        % solving for (2a) as well as (3a), making these assumptions:
        % * AS1: Assume manna-from-heaven and same tax in year-1 and year-2
        % of covid as under steady-state.
        % * AS2: Assume income loss is fully covered by in covid year-1 and
        % year-2, b=1, meaning xi does not matter
        % * AS3: For (2a) Approximating in effect income in 2019, which is
        % used to determine stimulus checks, with income in 2020, and in
        % effect ignoring one year kids transition probability as well.
        
        [mp_params, mp_controls, st_biden_or_trump, st_solu_type, ...
            bl_parfor, it_workers, ...
            bl_export, bl_load_mat, snm_suffix] = varargin{:};
        
    else
        error('Need to provide 8 parameter inputs');
    end
    
else
    clc;
    clear all;
    
    % 1a. Parfor controls
%      bl_parfor = true;
%      it_workers = 14;
     bl_parfor = true;
     it_workers = 4;
%     bl_parfor = false;
%     it_workers = 1;

    % 1b. Export Controls
    % bl_export = false;
    bl_export = true;
    bl_load_mat = false;
    

    % 1c. Solution Type
    % st_biden_or_trump = 'trumpchk';
    st_biden_or_trump = 'bidenchk';
    st_solu_type = 'bisec_vec';
        
    % 2. Set Up Parameters
    % Solve the VFI Problem and get Value Function        
    %     mp_params = snw_mp_param('default_moredense_a55zh43zs11');
    %     mp_params = snw_mp_param('default_moredense_a100zh266zs1');
    %       mp_params = snw_mp_param('default_moredense_a100zh266zs1');
    %     mp_params = snw_mp_param('default_moredense_a75zh101zs5');
    %     mp_params = snw_mp_param('default_moredense_a55z363');
        %     mp_params = snw_mp_param('default_moredense');
%       mp_params = snw_mp_param('default_dense');
%     mp_params = snw_mp_param('default_small', false, 'tauchen', false, 8, 8);
%       mp_params = snw_mp_param('default_tiny');
    
    
    % 2a. Simulation 1, e1m1, dense a and zh test
    % mp_params = snw_mp_param('default_moredense_a100zh266_e1m1');
    % 2b. Simulation 2, both edu and marriage, no spouse shock
    % mp_params = snw_mp_param('default_moredense_a100zh266_e2m2');
    % 2c. Simulation 3, both edu and marriage, 5 spouse shock
%     mp_params = snw_mp_param('default_moredense_a65zh81zs5_e2m2');
    % mp_params = snw_mp_param('default_moredense_a100zh81zs5_e2m2');   
    % mp_params = snw_mp_param('default_moredense_a65zh133zs5_e2m2');
%     mp_params = snw_mp_param('default_moredense_a65zh266zs5_e2m2', false, 'tauchen', false, 8, 8);
    mp_params = snw_mp_param('default_base', false, 'tauchen', false, 8, 8);
    
    % 3. Controls
    mp_controls = snw_mp_control('default_test');
    
    % 4. Unemployment
    % L283 LABEL B
    % set Unemployment Related Variables
    
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    TR=100/58056; % Value of a wezlfare check (can receive multiple checks). TO DO: Update with alternative values
    
    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_params('TR') = TR;
    
    % 5a. Check Count
    xi=0; % Full income loss if get covid shock
    b=1     ; % Fully unemployment insurance
    n_welfchecksgrid = 89; % 89 allows for double adults + double kids
    % n_welfchecksgrid = 245; % 245 max
    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_params('n_welfchecksgrid') = n_welfchecksgrid;    
    
    mp_params('a2_covidyr') = mp_params('a2_covidyr_manna_heaven');

    snm_suffix = ['_b1_xi0_mannna_' num2str(n_welfchecksgrid-1)];
%     snm_suffix = ['_b1_xi0_fullypay_' num2str(n_welfchecksgrid-1)];
    
    % 5b. Income Grid
    % Income Groups
    % 4 refers to 4*58056=232224 dollars in 2012USD
    % max 7 refers to 7*58056=406392 dollars in 2012USD
    % all phase out = (4400/5)*100 + 150000 = 238000
    % if 500 dollar interval, need 476 inc groups before 238000  
    % if have 85 percent of points betwen 238000, 
    fl_max_phaseout = 238000;
    fl_multiple = 58056;
    it_bin_dollar_before_phaseout = 500;
    it_bin_dollar_after_phaseout = 5000;
    fl_thres = fl_max_phaseout/fl_multiple;
    inc_grid1 = linspace(0,fl_thres,(fl_max_phaseout)/it_bin_dollar_before_phaseout);
    inc_grid2 = linspace(fl_thres, 7, (7*fl_multiple-fl_max_phaseout)/it_bin_dollar_after_phaseout); 
    inc_grid=sort(unique([inc_grid1 inc_grid2]'));    
%     n_incgrid=25; % Number of income groups
%     inc_grid=linspace(0,7,n_incgrid)';
    
    mp_params('n_incgrid') = length(inc_grid);
    mp_params('inc_grid') = inc_grid;
    
    % 6. Display Controls
    % Solve for Unemployment Values
    mp_controls('bl_print_vfi') = true;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = true;
    mp_controls('bl_print_ds_verbose') = true;
    mp_controls('bl_print_precompute') = true;
    mp_controls('bl_print_precompute_verbose') = false;
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jaeemk') = false;
    mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jmky') = false;
    mp_controls('bl_print_evuvw19_jmky_verbose') = false;    
    mp_controls('bl_print_evuvw19_jmky_mass') = false;
    mp_controls('bl_print_evuvw19_jmky_mass_verbose') = false;
            
end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_welfchecksgrid', 'n_jgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_welfchecksgrid, n_jgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'pi_unemp', 'n_incgrid', 'inc_grid'});
[pi_unemp, n_incgrid, inc_grid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, ...
    {'bl_print_evuvw19_jmky_allchecks', 'bl_print_evuvw19_jmky_allchecks_verbose'});
[bl_print_evuvw19_jmky_allchecks, bl_print_evuvw19_jmky_allchecks_verbose] = params_group{:};

%% Mat and CSV store Names

mp_path = snw_mp_path('fan');
snm_invoke_suffix = strrep(mp_params('mp_params_name'), 'default_', '');
snm_file = ['snwx_' char(st_biden_or_trump) '_' char(snm_invoke_suffix) char(snm_suffix)];
spn_csv_path = fullfile(mp_path('spt_simu_outputs'), [snm_file '.csv']);
spt_mat_path = fullfile(mp_path('spt_simu_outputs_mat'), [snm_file '.mat']);

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

% pt_all = snw_mp_path('fan');
% spt_mat_path = fullfile(pt_all('spt_simu_outputs_mat'), 'allchecks_vfi_dist.mat');
% spt_mat_path = fullfile('C:\Users\fan\Documents\temp', 'allchecks_vfi_dist_7.mat');

% if do not load mat or if the mat file does not exist
if (~bl_load_mat || ~isfile(spt_mat_path))
    %% A. Solve VFI
    % 2. Solve VFI and Distributon
    % Solve the Model to get V working and unemployed
    % solved with calibrated regular a2
    
    if strcmp(st_biden_or_trump, 'bchklock')
        invbtlock = mp_params('invbtlock');
        mp_params('invbtlock') = 1;
    end
    [V_ss,ap_ss,cons_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
    if strcmp(st_biden_or_trump, 'bchklock')
        mp_params('invbtlock') = invbtlock;
    end
    
%     mp_params('xi') = 1;
%     mp_params('b') = 0;
%     [V_ss_2020,~,cons_ss_2020,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
%     
%     % the difference should be zero if xi=1 and b=0
%     V_ss_diff = min(V_ss_2020 - V_ss, [], 'all');
%     if (V_ss_diff > 0)
%         error('V_ss_diff > 0');
%     end
    
    % 2020 if employed, same as steady state unless tax differs
    % 2020 V and C same as V_SS and cons_ss if tax the same
    if (mp_params('a2_covidyr') == mp_params('a2'))
        if strcmp(st_biden_or_trump, 'bchklock')
            % need to resolve for MIT Shock Year 1 with LOCKDOWN
            [V_ss_2020, ~, cons_ss_2020] = ...
                snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
        else
            % mana from heaven
            V_ss_2020 = V_ss;
            cons_ss_2020 = cons_ss;
        end
    else
        % change xi and b to for people without unemployment shock
        % solving for employed but 2020 tax results
        % a2_covidyr > a2, we increased tax in 2020 to pay for covid and other
        % costs resolve for both employed and unemployed
        % since resolving, if lockdown will have invbtlock used.
        xi = mp_params('xi');
        b = mp_params('b');
        mp_params('xi') = 1;
        mp_params('b') = 0;
        [V_ss_2020,~,cons_ss_2020,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
        mp_params('xi') = xi;
        mp_params('b') = b;
    end
    
    % 2020 if unemployed, can differ from steady state because xi and b
    % even if tax is the same. Solve unemployment, with three input
    % parameters, auto will use a2_covidyr as tax, similar for employed
    % call above
    [V_unemp_2020,~,cons_unemp_2020] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);
    
    %% B. Solve Dist COVID year 1 or COVID year 2 (with year one stimulus)   
       
    % Solve for steady-state distribution before COVID year start-of-year
    [Phi_true, Phi_adj_ss] = snw_ds_main_vec(mp_params, mp_controls, ap_ss, cons_ss);
    
    % One-shot MIT shock, distribution at the beginning of 2020 and in 2019
    % determined by non MIT year steady-state policy functions.
    if strcmp(st_biden_or_trump, 'trumpchk')
        % no actions needed
        % unlike in SNW_DS_MAIN_VEC_MULTITYPES, not solving for end of the
        % year post COVID distribution. Just need the policy functions
        % given lock-down which are already solved. 
        disp('Trump Check, do not need to resolve distribution')
        
    elseif (strcmp(st_biden_or_trump, 'bidenchk') || ...
            strcmp(st_biden_or_trump, 'bchklock') || ...
            strcmp(st_biden_or_trump, 'bchknoui'))
        % Use steady-state continuation value, to solve for optimal choices
        % given stimulus checks
        % Assume manna-from-heaven, (1) same tax in year-1 of covid as under
        % steady state. (2) assume income loss is fully covered by
        % "unemployment benefits". Because of this, do not have to solve
        % for unemployed and employed separately, their
        % resource-availability are identical. Approximating in effect
        % income in 2019, which is used to determine stimulus checks, with
        % income in 2020, and in effect ignoring one year kids transition
        % probability. Because of these assumptions, the effects of the
        % Trump-stimulus is to strictly increase resources availability for
        % households receiving stimulus, and no impact on those not
        % receiving stimulus. And those that do receive can save more,
        % leading to more savings than under steady-state. 
        
        disp('Biden Check, resolve for distributions given Trump check')
        
        xi = mp_params('xi');
        b = mp_params('b');
        mp_params('xi') = 0;
        mp_params('b') = 1;        
        % lockdown parameter potentially activated already here
        [~,ap_xi0b1_trumpchecks,cons_xi0b1_trumpchecks, mp_valpol_more_trumpchecks] = ...
            snw_vfi_main_bisec_vec_stimulus(mp_params, mp_controls, V_ss);
        mp_params('xi') = xi;
        mp_params('b') = b;
        
        % Distribution upon arrival in 2nd-covid year, the biden year, update Phi_true_ss
        [Phi_true] = snw_ds_main_vec(mp_params, mp_controls, ...
            ap_xi0b1_trumpchecks,cons_xi0b1_trumpchecks, ...
            mp_valpol_more_trumpchecks, ...
            Phi_adj_ss);
    else 
        error(['st_biden_or_trump=' char(st_biden_or_trump) ...
            ' is not allowed, has to be trumpchk or bidenchk or bchklock or bchknoui'])
    end
    
    %% C. Pre-compute
    % cl_st_precompute_list = {'a', ...
    %     'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid', ...
    %     'ap_idx_lower_ss', 'ap_idx_higher_ss', 'ap_idx_lower_weight_ss', ...
    %     'inc_tot_ygroup_grid', ''};
    cl_st_precompute_list = {'a', ...
        'inc', 'inc_unemp', 'spouse_inc', 'spouse_inc_unemp', 'ref_earn_wageind_grid', ...
        'inc_tot_ygroup_grid', 'ar_z_ctr_amz'};
    mp_controls('bl_print_precompute_verbose') = false;
    % Pre-Compute Matrixes and YMKY Mass
    % Pre-compute
    [mp_precompute_res] = snw_hh_precompute(mp_params, mp_controls, cl_st_precompute_list, ap_ss, Phi_true);
    
    % Load Results and Save to MAT
    ar_a_amz = mp_precompute_res('ar_a_amz');
    ar_inc_amz = mp_precompute_res('ar_inc_amz');
    ar_inc_unemp_amz = mp_precompute_res('ar_inc_unemp_amz');
    ar_spouse_inc_amz = mp_precompute_res('ar_spouse_inc_amz');
    ar_spouse_inc_unemp_amz = mp_precompute_res('ar_spouse_inc_unemp_amz');
    ref_earn_wageind_grid = mp_precompute_res('ref_earn_wageind_grid');
    inc_tot_ygroup_grid = mp_precompute_res('inc_tot_ygroup_grid');
    ar_z_ctr_amz = mp_precompute_res('ar_z_ctr_amz');
    
    %% D. YMKY Mass
    [Phi_true_jmky] = snw_evuvw19_jmky_mass(mp_params, mp_controls, Phi_true, inc_tot_ygroup_grid);
    
    %% save
    save(spt_mat_path, 'Phi_true', 'ap_ss', 'V_ss_2020', 'cons_ss_2020', 'V_unemp_2020', 'cons_unemp_2020',...
        'ar_a_amz', ...
        'ar_inc_amz', 'ar_inc_unemp_amz', 'ar_spouse_inc_amz', 'ar_spouse_inc_unemp_amz', 'ref_earn_wageind_grid', ...
        'inc_tot_ygroup_grid', 'ar_z_ctr_amz', ...
        'Phi_true_jmky');
    
else
    
%     % Load
%     load(spt_mat_path, 'Phi_true', 'ap_ss', 'V_ss_2020', 'cons_ss_2020', 'V_unemp_2020', 'cons_unemp_2020', ...
%         'ar_a_amz', ...
%         'ar_inc_amz', 'ar_inc_unemp_amz', 'ar_spouse_inc_amz', 'ar_spouse_inc_unemp_amz', 'ref_earn_wageind_grid', ...
%         'inc_tot_ygroup_grid', 'ar_z_ctr_amz', ...
%         'Phi_true_jmky');
%     
%     mp_precompute_res = containers.Map('KeyType', 'char', 'ValueType', 'any');
%     mp_precompute_res('ar_a_amz') = ar_a_amz;
%     mp_precompute_res('ar_inc_amz') = ar_inc_amz;
%     mp_precompute_res('ar_inc_unemp_amz') = ar_inc_unemp_amz;
%     mp_precompute_res('ar_spouse_inc_amz') = ar_spouse_inc_amz;
%     mp_precompute_res('ar_spouse_inc_unemp_amz') = spouse_inc_unemp;
%     mp_precompute_res('ref_earn_wageind_grid') = ref_earn_wageind_grid;
%     mp_precompute_res('inc_tot_ygroup_grid') = inc_tot_ygroup_grid;
%     mp_precompute_res('ar_z_ctr_amz') = ar_z_ctr_amz;
%     
%     clear ar_a_amz ar_inc_amz ar_inc_unemp_amz ar_spouse_inc_amz 
%     clear spouse_inc_unemp ref_earn_wageind_grid ar_z_ctr_amz
end


%% E1. Start Cluster
if (bl_parfor)
    delete(gcp('nocreate'));
    myCluster = parcluster('local');
    myCluster.NumWorkers = it_workers;
    parpool(it_workers);
end

%% E2. Start Check Loop
if (bl_print_evuvw19_jmky_allchecks)
    disp(strcat(['SNW_EVUVW19_JMKY_ALLCHECKS Start']));
end

ev19_jmky_allchecks=zeros(n_welfchecksgrid,n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);
ec19_jmky_allchecks=zeros(n_welfchecksgrid,n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

if (bl_parfor)
    
    mp_controls_parfor = mp_controls;
%     mp_controls_parfor('bl_timer') = false;
%     mp_controls_parfor('bl_print_evuvw19_jmky_allchecks') = false;    
     
    parfor i=1:n_welfchecksgrid        
        
        % check
        welf_checks=i-1;
        
        % Run
        [ev19_jmky_check, ec19_jmky_check] = ...
            snw_evuvw19_jmky_check(welf_checks, st_solu_type, ...
            mp_params, mp_controls_parfor, ...
            spt_mat_path, bl_print_evuvw19_jmky_allchecks);

        % Store Results
        ev19_jmky_allchecks(i,:,:,:,:) = ev19_jmky_check;
        ec19_jmky_allchecks(i,:,:,:,:) = ec19_jmky_check;

    end    
        
else
    
    for welf_checks=0:(n_welfchecksgrid-1)

        % Run
        [ev19_jmky_check, ec19_jmky_check] = ...
            snw_evuvw19_jmky_check(welf_checks, st_solu_type, ...
            mp_params, mp_controls, ...
            spt_mat_path, bl_print_evuvw19_jmky_allchecks);

        % Store Results
        ev19_jmky_allchecks(welf_checks+1,:,:,:,:) = ev19_jmky_check;
        ec19_jmky_allchecks(welf_checks+1,:,:,:,:) = ec19_jmky_check;

    end    
end

if (bl_parfor)
   delete(gcp('nocreate'));
end

%% F1. Mass to Probability
if exist('spt_mat_path','var')
    load(spt_mat_path, 'Phi_true_jmky');
end

Phi_true_jmky_prob = Phi_true_jmky./sum(Phi_true_jmky, 'all');

if exist('spt_mat_path','var')
    clear Phi_true_jmky
end

%% F2. Output for computing optimal allocation
Output=zeros((n_jgrid-1)*n_marriedgrid*n_kidsgrid*n_welfchecksgrid*n_incgrid,9);
counter=0;
inc_grid = [inc_grid; 10E30];
for inc_group=1:n_incgrid
    for kids=1:n_kidsgrid % Number of kids
        for married=1:n_marriedgrid % Marital status
            for j=1:(n_jgrid-1) % Age
                for welf_checks=0:(n_welfchecksgrid-1) % Number of welfare checks
                    
                    counter=counter+1;
                    
                    Output(counter,1)=17+j;
                    Output(counter,2)=married-1;
                    Output(counter,3)=kids-1;
                    Output(counter,4)=welf_checks;
                    
                    Output(counter,5)=inc_grid(inc_group);
                    Output(counter,6)=Phi_true_jmky_prob(j,married,kids,inc_group);
                    Output(counter,7)=psi(j);
                    
                end
            end
        end
    end
end
Output(:,8) = ev19_jmky_allchecks(:);
Output(:,9) = ec19_jmky_allchecks(:);

%% F2. Drop Zero Mass Rows
Output = Output(Output(:,6) > 0,:);

%% G. Save File
if (bl_export)
    writematrix(Output, spn_csv_path);
end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete_vu_vw_checks = strjoin(...
        ["Completed SNW_EVUVW19_JMKY_ALLCHECKS", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");    
    disp(st_complete_vu_vw_checks);
end

%% Print
if (bl_print_evuvw19_jmky_allchecks_verbose)
    mp_outcomes = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_outcomes('ev19_jmky_allchecks') = ev19_jmky_allchecks;
    mp_outcomes('ec19_jmky_allchecks') = ec19_jmky_allchecks;
    mp_outcomes('ev19_jmky_allchecks_posmass') = Output(:,8);
    mp_outcomes('ec19_jmky_allchecks_posmass') = Output(:,9);
    mp_outcomes('Output') = Output;
    ff_container_map_display(mp_outcomes, 10, 7);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = ev19_jmky_allchecks;
    elseif (it_k==2)
        ob_out_cur = ec19_jmky_allchecks;
    elseif (it_k==3)
        ob_out_cur = Output;
    end
    varargout{it_k} = ob_out_cur;
end

end

function [ev19_jmky_check, ec19_jmky_check] = ...
    snw_evuvw19_jmky_check(welf_checks, st_solu_type, ...
        mp_params, mp_controls, ...
        spt_mat_path, bl_print_evuvw19_jmky_allchecks) 

    if (bl_print_evuvw19_jmky_allchecks) 
        disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX');        
        tm_start_check = tic; 
    end    
        
    % Solve ev 19 JAEEMK
    [ev19_jaeemk_check, ec19_jaeemk_check] = snw_evuvw19_jaeemk_foc(...
        welf_checks, st_solu_type, ...
        mp_params, mp_controls, ...
        spt_mat_path);
%     clear V_ss_2020 ap_ss cons_ss_2020 V_unemp_2020 cons_unemp_2020
    
    % Solve ev 19 JMKY    
    load(spt_mat_path, 'Phi_true', 'inc_tot_ygroup_grid', 'Phi_true_jmky');
    [ev19_jmky_check, ec19_jmky_check] = snw_evuvw19_jmky(...
        mp_params, mp_controls, ...
        ev19_jaeemk_check, ec19_jaeemk_check, ...
        Phi_true, Phi_true_jmky, inc_tot_ygroup_grid);
    clear Phi_true Phi_true_jmky ev19_jaeemk_check ec19_jaeemk_check
    
    if (bl_print_evuvw19_jmky_allchecks)
        tm_end_check = toc(tm_start_check);
        disp(strcat([...            
            'SNW_EVUVW19_JMKY_ALLCHECKS: Finished Check ' ...
            num2str(welf_checks) ' of ' ...
            num2str(mp_params('n_welfchecksgrid')-1) ...
            ', time=' num2str(tm_end_check)]));
    end
    
end