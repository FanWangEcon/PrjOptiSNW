%% SNW_EVUVW19_JMKY_ALLCHECKS 2019 Planner Values All Checks
%    Call the SNW_EVUVW19_JMKY function to find the Value at All Checks.
%    Solving for one check at a time and storing resulting jmky results one
%    check at a time. Uses pre-stored matrixes that are re-used across
%    check calculations
%
%    [EV19_JMKY_ALLCHECKS, EC19_JMKY_ALLCHECKS, OUTPUT] =
%    SNW_EVUVW19_JMKY_ALLCHECKS(WELF_CHECKS, ST_SOLU_TYPE, MP_PARAMS,
%    MP_CONTROLS, V_SS, V_UNEMP) provide V_SS and V_UNEMP solved out
%    elsewhere, and only get two outputs out.
%
%    See also SNW_EVUVW19_JMKY, SNW_EVUVW19_JMKY_MASS, SNW_EVUVW19_JAEEMK,
%    SNW_EVUVW20_JAEEMK, SNW_HH_PRECOMPUTE
%

%%
function [varargout]=snw_evuvw19_jmky_allchecks(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==7)
        [mp_params, mp_controls, ...
            bl_parfor, it_workers, ...
            bl_export, snm_suffix] = varargin{:};
    else
        error('Need to provide 7 parameter inputs');
    end
    
else
    clc;
    clear all;

    bl_parfor = true;
    it_workers = 2;
    bl_export = true;
    
%     bl_parfor = false;
%     it_workers = 1;
%     bl_export = true;
    
    snm_suffix = '';
    
    st_solu_type = 'bisec_vec';
    
    % 1. Set Up Parameters
    % Solve the VFI Problem and get Value Function        
%      mp_params = snw_mp_param('default_moredense_a100zh266_e0m0');
%     mp_params = snw_mp_param('default_moredense_a55zh43zs11');
    mp_params = snw_mp_param('default_moredense_a100zh266zs1');
%     mp_params = snw_mp_param('default_moredense_a75zh101zs5');        
%     mp_params = snw_mp_param('default_moredense_a55z363');
    %     mp_params = snw_mp_param('default_moredense');
    %     mp_params = snw_mp_param('default_dense');
%         mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    
    % set Unemployment Related Variables
    xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
    b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
    TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
    n_welfchecksgrid = 45;
    
    mp_params('xi') = xi;
    mp_params('b') = b;
    mp_params('TR') = TR;
    mp_params('n_welfchecksgrid') = n_welfchecksgrid;
    
    % Income Groups
    % 4 refers to 4*58056=232224 dollars in 2012USD
    % max 7 refers to 7*58056=406392 dollars in 2012USD
    % all phase out = (4400/5)*100 + 150000 = 238000
    % if 500 dollar interval, need 476 inc groups before 238000  
    % if have 85 percent of points betwen 238000, 
    fl_max_phaseout = 238000;
    fl_multiple = 58056;
    it_bin_dollar_before_phaseout = 500;
    it_bin_dollar_after_phaseout = 2000;
    fl_thres = fl_max_phaseout/fl_multiple;
    inc_grid1 = linspace(0,fl_thres,(fl_max_phaseout)/it_bin_dollar_before_phaseout);
    inc_grid2 = linspace(fl_thres, 7, (7*fl_multiple-fl_max_phaseout)/it_bin_dollar_after_phaseout); 
    inc_grid=sort(unique([inc_grid1 inc_grid2]'));
    
    % n_incgrid=25; % Number of income groups
    % inc_grid=linspace(0,7,n_incgrid)';
    
    mp_params('n_incgrid') = length(inc_grid);
    mp_params('inc_grid') = inc_grid;
    
    % Solve for Unemployment Values
    mp_controls('bl_print_vfi') = false;
    mp_controls('bl_print_vfi_verbose') = false;
    mp_controls('bl_print_ds') = false;
    mp_controls('bl_print_ds_verbose') = false;
    mp_controls('bl_print_precompute') = false;
    mp_controls('bl_print_precompute_verbose') = false;
    mp_controls('bl_print_a4chk') = false;
    mp_controls('bl_print_a4chk_verbose') = false;
    mp_controls('bl_print_evuvw20_jaeemk') = false;
    mp_controls('bl_print_evuvw20_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jaeemk') = false;
    mp_controls('bl_print_evuvw19_jaeemk_verbose') = false;
    mp_controls('bl_print_evuvw19_jmky') = false;
    mp_controls('bl_print_evuvw19_jmky_verbose') = false;
        
end

%% Parse Model Parameters
params_group = values(mp_params, ...
    {'n_welfchecksgrid', 'n_jgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_welfchecksgrid, n_jgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

params_group = values(mp_params, {'n_incgrid'});
[n_incgrid] = params_group{:};

%% Parse Model Controls
% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, ...
    {'bl_print_evuvw19_jmky_allchecks', 'bl_print_evuvw19_jmky_allchecks_verbose'});
[bl_print_evuvw19_jmky_allchecks, bl_print_evuvw19_jmky_allchecks_verbose] = params_group{:};

%% Timing and Profiling Start
if (bl_timer)
    tm_start = tic;
end

%% A. Solve VFI
% 2. Solve VFI and Distributon
% Solve the Model to get V working and unemployed
[V_ss,ap_ss,cons_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
% Solve unemployment
[V_unemp,~,cons_unemp] = snw_vfi_main_bisec_vec(mp_params, mp_controls, V_ss);

%% B. Solve Dist
[Phi_true] = snw_ds_main(mp_params, mp_controls, ap_ss, cons_ss);

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
inc_tot_ygroup_grid = mp_precompute_res('inc_tot_ygroup_grid');

%% D. YMKY Mass
[Phi_true_jmky] = snw_evuvw19_jmky_mass(mp_params, mp_controls, Phi_true, inc_tot_ygroup_grid);

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
            V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, ...
            Phi_true, Phi_true_jmky, inc_tot_ygroup_grid,...
            mp_precompute_res, bl_print_evuvw19_jmky_allchecks);

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
            V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, ...
            Phi_true, Phi_true_jmky, inc_tot_ygroup_grid,...
            mp_precompute_res, bl_print_evuvw19_jmky_allchecks);

        % Store Results
        ev19_jmky_allchecks(welf_checks+1,:,:,:,:) = ev19_jmky_check;
        ec19_jmky_allchecks(welf_checks+1,:,:,:,:) = ec19_jmky_check;

    end    
end


if (bl_parfor)
   delete(gcp('nocreate'));
end

%% F1. Output for computing optimal allocation
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
                    Output(counter,6)=Phi_true_jmky(j,married,kids,inc_group);
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
    mp_path = snw_mp_path('fan');
    snm_invoke_suffix = strrep(mp_params('mp_params_name'), 'default_', '');
    snm_file_csv = ['snwx_v_planner_' char(snm_invoke_suffix) char(snm_suffix) '.csv'];
    writematrix(Output, [mp_path('spt_simu_outputs') snm_file_csv]);
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
        V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, ...
        Phi_true, Phi_true_jmky, inc_tot_ygroup_grid,...
        mp_precompute_res, bl_print_evuvw19_jmky_allchecks) 

    if (bl_print_evuvw19_jmky_allchecks) 
        disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX');        
        tm_start_check = tic; 
    end
    
    % Solve ev 19 JAEEMK
    [ev19_jaeemk_check, ec19_jaeemk_check] = snw_evuvw19_jaeemk_foc(...
        welf_checks, st_solu_type, ...
        mp_params, mp_controls, ...
        V_ss, ap_ss, cons_ss, V_unemp, cons_unemp, mp_precompute_res);
    
    % Solve ev 19 JMKY
    [ev19_jmky_check, ec19_jmky_check] = snw_evuvw19_jmky(...
        mp_params, mp_controls, ...
        ev19_jaeemk_check, ec19_jaeemk_check, ...
        Phi_true, Phi_true_jmky, inc_tot_ygroup_grid);
        
    if (bl_print_evuvw19_jmky_allchecks)
        tm_end_check = toc(tm_start_check);
        disp(strcat([...            
            'SNW_EVUVW19_JMKY_ALLCHECKS: Finished Check ' ...
            num2str(welf_checks) ' of ' ...
            num2str(mp_params('n_welfchecksgrid')-1) ...
            ', time=' num2str(tm_end_check)]));
    end
    
end