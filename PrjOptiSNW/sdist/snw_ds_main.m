%% SNW_DS_MAIN Simulates Distributions
%    Given policy functions, simulate distributions
%
%    Pref, Technology, and prices SCALARS:
%
%    * AP_SS ndarry of policy function solution from SNW_VFI_MAIN
%    * MP_PARAMS contain map of parameters
%    * MP_CONTROLS container map of parameters
%
%    [MP_DSVFI_RESULTS] = SNW_DS_MAIN() invoke model with externally set
%    parameter map and control map. Results outputed to a map containing
%    various output matrixes. 
%
%    See also SNWX_DS_MAIN, SNW_VFI_MAIN, SNW_MP_CONTROLS, SNW_MP_PARAM
%

%%
function mp_dsvfi_results=snw_ds_main(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))
    
    if (length(varargin)==1)
        [mp_params] = varargin{:};
        mp_controls = snw_mp_controls('default_base');
    elseif (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
        [v_ss, ap_ss, c_ss, exitflag_VFI] = snw_vfi_main(mp_params, mp_controls);        
    elseif (length(varargin)==3)
        [mp_params, mp_controls, ap_ss] = varargin{:};
    end
    
else
    
    mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    [v_ss, ap_ss, c_ss, exitflag_VFI] = snw_vfi_main(mp_params, mp_controls);
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global a2 g_cons agrid SS pi_eta pi_kids Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

%% Parse Model Parameters
params_group = values(mp_params, {'g_cons', 'a2'});
[g_cons, a2] = params_group{:};

params_group = values(mp_params, {'agrid'});
[agrid] = params_group{:};

params_group = values(mp_params, {'pi_eta', 'pi_kids'});
[pi_eta, pi_kids] = params_group{:};

params_group = values(mp_params, {'SS'});
[SS] = params_group{:};

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
params_group = values(mp_controls, {'bl_compute_drv_stats'});
[bl_compute_drv_stats] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_ds', 'bl_print_ds_verbose'});
[bl_print_ds, bl_print_ds_verbose] = params_group{:};

% Store Controls
params_group = values(mp_controls, {'bl_vfi_store_all', 'bl_ds_store_all'});
[bl_vfi_store_all, bl_ds_store_all] = params_group{:};

%% Aggregation
Phiss=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

% Assume everyone starts with 0 assets
% Use stationary distribution of productivity shock and initial
% distribution for educational attainment, marital status, and number of
% kids from PSID
for eta=1:n_etagrid % Productivity
    for educ=1:n_educgrid % Fixed effects
        for married=1:n_marriedgrid % Marital status
            for kids=1:n_kidsgrid % No. of kids
                Phiss(1,1,eta,educ,married,kids)=stat_distr_eta(eta)*stat_distr_educ(educ)*stat_distr_married(educ,married)*stat_distr_kids(educ,married,kids);
            end
        end
    end
end

% Use policy functions and survival probabilities to get distribution of remaining idiosyncratic states
for j=1:(n_jgrid-1) % Age
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % No. of kids
                        
                        if ap_ss(j,a,eta,educ,married,kids)==0
                            inds(1)=1;
                            inds(2)=1;
                            vals(1)=1;
                            vals(2)=0;
                            
                        elseif ap_ss(j,a,eta,educ,married,kids)==agrid(n_agrid)
                            inds(1)=n_agrid;
                            inds(2)=n_agrid;
                            vals(1)=1;
                            vals(2)=0;
                            
                        else
                            
                            ind_aux=find(agrid<=ap_ss(j,a,eta,educ,married,kids),1,'last');
                            
                            inds(1)=ind_aux;
                            inds(2)=ind_aux+1;
                            
                            % Linear interpolation
                            vals(1)=1-((ap_ss(j,a,eta,educ,married,kids)-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                            vals(2)=1-vals(1);
                            
                        end
                        
                        for etap=1:n_etagrid
                            for kidsp=1:n_kidsgrid
                                Phiss(j+1,inds(1),etap,educ,married,kidsp)=Phiss(j+1,inds(1),etap,educ,married,kidsp)+Phiss(j,a,eta,educ,married,kids)*vals(1)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                                Phiss(j+1,inds(2),etap,educ,married,kidsp)=Phiss(j+1,inds(2),etap,educ,married,kidsp)+Phiss(j,a,eta,educ,married,kids)*vals(2)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                            end
                        end
                        
                    end
                end
            end
        end
    end
end

% Normalize distribution of idiosyncratic states to sum to 1 for each age
% Phi_true = joint probability mass function all states
% Phi_adj  = conditional probability mass function given age
Phi_adj=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
Phi_true=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid
    
    dummy=sum(sum(sum(sum(sum(Phiss(j,:,:,:,:,:))))));
    
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid
                        
                        if dummy>0
                            Phi_adj(j,a,eta,educ,married,kids)=Phiss(j,a,eta,educ,married,kids)/dummy;
                        else
                            Phi_adj(j,a,eta,educ,married,kids)=0;
                        end
                        
                        Phi_true(j,a,eta,educ,married,kids)=Phi_adj(j,a,eta,educ,married,kids)*Pop(j);
                        
                    end
                end
            end
        end
    end
    
end


% Check if the upper bound on assets binds
check_asset_distr=sum(sum(sum(sum(sum(Phi_true(:,n_agrid,:,:,:,:))))));

if (bl_print_ds || bl_ds_store_all)
    name='SNW_DS_MAIN: Share of population with assets equal to upper bound on asset grid=';
    st_pop_max_asset = [name,num2str(check_asset_distr/sum(Pop))];
    disp(st_pop_max_asset);
end

% Aggregate variables
A_agg=0;
Y_inc_agg=0;
Tax_revenues=0;
SS_spend=0;

for j=1:n_jgrid
    for a=1:n_agrid
        for eta=1:n_etagrid
            for educ=1:n_educgrid
                for married=1:n_marriedgrid
                    for kids=1:n_kidsgrid
                        
                        A_agg=A_agg+Phi_true(j,a,eta,educ,married,kids)*agrid(a); % Aggregate wealth
                        
                        [inc_aux,earn]=individual_income(j,a,eta,educ);
                        spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                        
                        Y_inc_agg=Y_inc_agg+Phi_true(j,a,eta,educ,married,kids)*( inc_aux+(married-1)*spouse_inc ); % Aggregate income
                        
                        %                      inc_aux=r*agrid(a)+epsilon(j,educ)*theta*exp(eta_grid(eta))+SS(j,educ); % Income excluding spousal income (if married)
                        Tax_revenues=Tax_revenues+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc_aux,(married-1)*spouse_inc)); % Tax revenues
                        
                        SS_spend=SS_spend+Phi_true(j,a,eta,educ,married,kids)*SS(j,educ); % Total spending on Social Security
                        
                    end
                end
            end
        end
    end
end

% Update guess for a2 (determines average level of income taxation)
% Assuming government balances its budget period-by-period

tol=10^-4;
err=abs((Tax_revenues/(SS_spend+g_cons*Y_inc_agg))-1);

a2_update=a2;

while err>tol
    
    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
        for a=1:n_agrid
            for eta=1:n_etagrid
                for educ=1:n_educgrid
                    for married=1:n_marriedgrid
                        for kids=1:n_kidsgrid
                            
                            [inc_aux,earn]=individual_income(j,a,eta,educ);
                            spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                            Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,a,eta,educ,married,kids)*max(0,Tax(inc_aux,(married-1)*spouse_inc)); % Tax revenues
                            
                        end
                    end
                end
            end
        end
    end
    
    a2=a2*(((SS_spend+g_cons*Y_inc_agg)/Tax_revenues_aux)^0.5); % Find value of a2 that balances government budget
    
    err=abs((Tax_revenues_aux/(SS_spend+g_cons*Y_inc_agg))-1);
        
    if (bl_print_ds_verbose)
        disp(['SNW_DS_MAIN: err=abs((TAX/(SS+g_cons*Y))-1)=' num2str(err)])
    end
end

a2_update=[a2_update,a2];

if (bl_print_ds || bl_ds_store_all)
    name='SNW_DS_MAIN: Old and updated value of a2=';
    st_a2_update=[name,num2str(a2_update)];
    disp(st_a2_update);
end

%% Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_ds = strjoin(...
        ["Completed SNW_DS_MAIN", ...
         ['SNW_MP_PARAM=' mp_params('mp_params_name')], ...
         ['SNW_MP_CONTROL=' mp_controls('mp_params_name')] ...
        ], ";");
    disp(st_complete_ds);
end

%% Collect Results
mp_dsvfi_results = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_dsvfi_results('mp_params') = mp_params;
mp_dsvfi_results('mp_controls') = mp_controls;
mp_dsvfi_results('ap_ss') = ap_ss;

mp_dsvfi_results('Phi_true') = Phi_true;
mp_dsvfi_results('Phi_adj') = Phi_adj;
mp_dsvfi_results('A_agg') = A_agg;
mp_dsvfi_results('Y_inc_agg') = Y_inc_agg;

% Store all of DS Function Outputs
if (bl_ds_store_all)    
    % Numeric Values
    mp_dsvfi_results('dummy') = dummy;
    mp_dsvfi_results('check_asset_distr') = check_asset_distr;
    mp_dsvfi_results('Tax_revenues') = Tax_revenues;    
    mp_dsvfi_results('SS_spend') = SS_spend;    
    mp_dsvfi_results('err') = err;
    
    % Strings
    mp_dsvfi_results('st_pop_max_asset') = string(st_pop_max_asset);
    mp_dsvfi_results('st_a2_update') = string(st_a2_update);
end

% Store all of VFI Function Outputs
if (bl_vfi_store_all && (length(varargin)<3))
    mp_dsvfi_results('v_ss') = v_ss;    
    mp_dsvfi_results('c_ss') = c_ss;
    
    % Household Asset State Space Value
    a_ss = zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    for a=1:n_agrid
        a_ss(:,a,:,:,:,:) = agrid(a);
    end
    mp_dsvfi_results('a_ss') = a_ss;
    
    % Households Population At State-Space
    n_ss = zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
    % add kids
    for it_kids=2:n_kidsgrid
        n_ss(:,:,:,:,:,it_kids) = it_kids-1;
    end
    % marry
    n_ss(:,:,:,:,2,:) = n_ss(:,:,:,:,2,:) + 1;
    % add own
    n_ss = n_ss + 1;
    mp_dsvfi_results('n_ss') = n_ss;

    mp_dsvfi_results('exitflag_VFI') = exitflag_VFI;
end

%% Results Basic Print and Verbose print
if (bl_print_ds_verbose)
    ff_container_map_display(mp_params);
    ff_container_map_display(mp_controls);
    ff_container_map_display(mp_dsvfi_results);
end

if (bl_compute_drv_stats && (length(varargin)<3))
    % Array Inputs
    mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
    mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss(:), zeros(1)};
    mp_cl_ar_xyz_of_s('a_ss') = {a_ss(:), zeros(1)};
    mp_cl_ar_xyz_of_s('c_ss') = {c_ss(:), zeros(1)};
    mp_cl_ar_xyz_of_s('v_ss') = {v_ss(:), zeros(1)};    
    mp_cl_ar_xyz_of_s('n_ss') = {n_ss(:), zeros(1)};
    mp_cl_ar_xyz_of_s('ar_st_y_name') = ["a_ss", "ap_ss", "c_ss", "v_ss", "n_ss"];

    % controls
    mp_support = containers.Map('KeyType','char', 'ValueType','any');
    mp_support('ar_fl_percentiles') = [0.01 0.1 1 5 10 20 25 30 40 50 60 70 75 80 90 95 99 99.9 99.99];
    if (bl_print_ds)
        mp_support('bl_display_final') = true;
    else
        mp_support('bl_display_final') = false;
    end
    mp_support('bl_display_detail') = false;    
    mp_support('bl_display_drvm2outcomes') = false;  
    mp_support('bl_display_drvstats') = false;
    mp_support('bl_display_drvm2covcor') = false;

    % Call Function     
    mp_cl_mt_xyz_of_s = ff_simu_stats(Phi_true(:)/sum(Phi_true,'all'), mp_cl_ar_xyz_of_s, mp_support);
    if (bl_vfi_store_all)
        mp_dsvfi_results('mp_cl_mt_xyz_of_s') = mp_cl_mt_xyz_of_s;
    end
end

end