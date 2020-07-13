%% SNW_DS_MAIN_GRID_SEARCH Simulates Distributions
%    Given policy functions, simulate distributions. When VFI results are
%    not provided, will solve for VFI using SNW_VFI_MAIN_BISEC_VEC.
%
%    Pref, Technology, and prices SCALARS:
%
%    * AP_SS ndarry of policy function solution from SNW_VFI_MAIN
%    * MP_PARAMS contain map of parameters
%    * MP_CONTROLS container map of parameters
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] = SNW_DS_MAIN()
%    invoke model with externally set parameter map and control map.
%    Results outputed to a map containing various output matrixes in
%    mp_dsvfi_results, and also distributional matrixes.
%
%    See also SNWX_DS_MAIN, SNW_VFI_MAIN, SNW_MP_CONTROLS, SNW_MP_PARAM
%

%%
function varargout=snw_ds_main_grid_search(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))

    if (length(varargin)==2)        
        [mp_params, mp_controls] = varargin{:};                
    elseif (length(varargin)==4)        
        % This will not produce extra statistics outputs
        [mp_params, mp_controls, ap_ss, cons_ss] = varargin{:};        
    elseif (length(varargin)==5)        
        % This will not produce extra statistics outputs
        [mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss] = varargin{:};        
    else
        error('Need to provide 2/4 parameter inputs');
    end
    
    bl_ds_store_all = false;
    
else

    mp_params = snw_mp_param('default_small');
    mp_controls = snw_mp_control('default_test');
    [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_grid_search(mp_params, mp_controls);
    
    bl_ds_store_all = true;
    
end

%% Solve VFI if not Provided

% output more
if (nargout >= 6)
    bl_ds_store_all = true;
end

% requesting all stats outputs, but does not have vfi more info
if (bl_ds_store_all && ~exist('mp_valpol_more_ss','var'))
    [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_grid_search(mp_params, mp_controls);
end

% requesting subset outputs if don't have policy function, and does not
% request extra outputs
if (~exist('ap_ss','var'))
    [v_ss, ap_ss, cons_ss] = snw_vfi_main_grid_search(mp_params, mp_controls);
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global a2 g_cons agrid SS pi_eta pi_kids Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid 
% Added
global epsilon theta eta_H_grid eta_S_grid r g_n psi Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r', 'g_n', 'g_cons', 'a2'});
[theta, r, g_n, g_cons, a2] = params_group{:};

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

%% Initialize 6D Distributional Array and Initial Distribution

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

                       for etap=1:n_etagrid
                           for kidsp=1:n_kidsgrid
                               Phiss(j+1,ap_ss(j,a,eta,educ,married,kids),etap,educ,married,kidsp)=Phiss(j+1,ap_ss(j,a,eta,educ,married,kids),etap,educ,married,kidsp)+Phiss(j,a,eta,educ,married,kids)*pi_eta(eta,etap)*pi_kids(kids,kidsp,j,educ,married);
                           end
                       end

                   end
               end
           end
       end
   end 
end

% Normalize distribution of idiosyncratic states to sum to 1 for each age
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
if (bl_print_ds)
    disp(strcat(['SNW_DS_MAIN_GRID_SEARCH: Share of population with assets equal to upper bound on asset grid:' ...
        num2str(check_asset_distr/sum(Pop))]));
end

% Aggregate variables
A_agg=0;
Aprime_agg=0;
C_agg=0;
Y_inc_agg=0;
Tax_revenues=0;
SS_spend=0;
Bequests_aux=0;

for j=1:n_jgrid
   for eta=1:n_etagrid
       for educ=1:n_educgrid
           for married=1:n_marriedgrid
               for kids=1:n_kidsgrid

                   A_agg=A_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*agrid(1:n_agrid); % Aggregate wealth
                   Aprime_agg=Aprime_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*agrid(ap_ss(j,1:n_agrid,eta,educ,married,kids)); % Aggregate saving
                   C_agg=C_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*cons_ss(j,1:n_agrid,eta,educ,married,kids)'; % Aggregate consumption

                   [inc,earn]=individual_income(j,1:n_agrid,eta,educ);
                   spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));

                   inc_aux=r*(agrid(1:n_agrid)+Bequests*(bequests_option-1))+epsilon(j,educ)*theta*exp(eta_H_grid(eta)); % Income (excluding Social Security benefits)

                   Y_inc_agg=Y_inc_agg+Phi_true(j,1:n_agrid,eta,educ,married,kids)*( inc_aux+(married-1)*spouse_inc*exp(eta_S_grid(eta)) ); % Aggregate income (labor earnings, spousal income, interest earnings)

                   Tax_revenues=Tax_revenues+Phi_true(j,1:n_agrid,eta,educ,married,kids)*max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)))); % Tax revenues

                   SS_spend=SS_spend+sum(Phi_true(j,1:n_agrid,eta,educ,married,kids)*SS(j,educ)); % Total spending on Social Security

                   Bequests_aux=Bequests_aux+Phi_true(j,1:n_agrid,eta,educ,married,kids)*agrid(ap_ss(j,1:n_agrid,eta,educ,married,kids))*(1-psi(j)); % Accidental Bequests*(bequests_option-1)

               end                       
           end
       end
   end
end

% Alternative 1: Suppose accidental bequests go to the government
% Alternative 2: Allocate accidental bequests uniformly across the population
if bequests_option==1
    if throw_in_ocean==0
        if (bl_print_ds)
            disp('SNW_DS_MAIN: Accidental bequests are part of government revenues');
        end
        Tax_revenues=Tax_revenues+Bequests_aux*(1+r);
    elseif throw_in_ocean==1
        if (bl_print_ds)
            disp('SNW_DS_MAIN: Accidental bequests are thrown in the ocean');
        end
    end
elseif bequests_option==2
    if (bl_print_ds)
        disp('SNW_DS_MAIN: Accidental bequests are uniformly distributed across the population');
    end
    Bequests=Bequests_aux/(sum(Pop)*(1+g_n));
end

% Update guess for a2 (determines average level of income taxation)
% Assuming government balances its budget period-by-period

tol=5*10^-4; % 10^-4;
err=abs((Tax_revenues/(SS_spend+g_cons*Y_inc_agg))-1);

a2_update=a2;

it=0;

while err>tol
    
    it=it+1;
    
    Tax_revenues_aux=0;
    
    for j=1:n_jgrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid

                       [inc,earn]=individual_income(j,1:n_agrid,eta,educ);
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                       Tax_revenues_aux=Tax_revenues_aux+Phi_true(j,1:n_agrid,eta,educ,married,kids)*max(0,Tax(inc,(married-1)*spouse_inc*exp(eta_S_grid(eta)))); % Tax revenues

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
    
    a2=a2*(((SS_spend+g_cons*Y_inc_agg)/Tax_revenues_aux)^0.25); % Find value of a2 that balances government budget
    
    err=abs((Tax_revenues_aux/(SS_spend+g_cons*Y_inc_agg))-1);
    
    if (bl_print_ds)
        st_tax_iter = strjoin(...
            ["SNW_DS_MAIN tax and spend", ...
             ['it=' num2str(it)], ...
             ['err=' num2str(err)] ...
            ], ";");        
        disp(st_tax_iter);
    end    
end

a2=a2*0.75+a2_update*0.25;

if (bl_print_ds)
    name='SNW_DS_MAIN_GRID_SEARCH: Number of a2-adjustments (for taxation) used to balance the government budget= ';
    name2=[name,num2str(it)];
    disp(name2);

    a2_update=[a2_update,a2];

    name='SNW_DS_MAIN_GRID_SEARCH: Old and updated value of a2=';
    name2=[name,num2str(a2_update)];
    disp(name2);

    name='SNW_DS_MAIN_GRID_SEARCH: Aggregates: Cons., Gov. cons., Save, Assets, Income, Bequests ';
    aggregates=[C_agg,g_cons*Y_inc_agg,A_agg,Aprime_agg,Y_inc_agg,Bequests_aux];
    name2=[name,num2str(aggregates)];
    disp(name2);

    name='SNW_DS_MAIN_GRID_SEARCH: Resource constraint: C_t+A_{t+1}+G_t=A_t+Y_t ';
    name2=[name,num2str([C_agg+g_cons*Y_inc_agg+Aprime_agg,A_agg+Y_inc_agg])];
    disp(name2);    
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

%% Collect Results if Require Extra Outputs
if (bl_ds_store_all)
    mp_dsvfi_results = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_dsvfi_results('mp_params') = mp_params;
    mp_dsvfi_results('mp_controls') = mp_controls;
    
    % Store Policy Functions as Outputs
    mp_dsvfi_results('v_ss') = v_ss;
    mp_dsvfi_results('ap_ss') = ap_ss;
    mp_dsvfi_results('cons_ss') = cons_ss;
    
    y_inc_ss = mp_valpol_more_ss('inc_VFI');
    y_earn_ss = mp_valpol_more_ss('earn_VFI');
    y_spouse_inc_ss = mp_valpol_more_ss('spouse_inc_VFI');
    SS_ss = mp_valpol_more_ss('SS_VFI');
    tax_ss = mp_valpol_more_ss('tax_VFI');
    mp_dsvfi_results('y_head_inc_ss') = y_inc_ss;
    mp_dsvfi_results('y_head_earn_ss') = y_earn_ss;
    mp_dsvfi_results('y_spouse_inc_ss') = y_spouse_inc_ss;
    mp_dsvfi_results('SS_ss') = SS_ss;
    mp_dsvfi_results('tax_ss') = tax_ss;       
    
    % Additional aggregate Statistics
    mp_dsvfi_results('Aprime_agg') = Aprime_agg;
    mp_dsvfi_results('C_agg') = C_agg;
    mp_dsvfi_results('Tax_revenues') = Tax_revenues;
    mp_dsvfi_results('SS_spend') = SS_spend;
    mp_dsvfi_results('Bequests_aux') = Bequests_aux;
   
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

    % Capital/interest Earning of Income
    yshr_interest = (a_ss*r)./(y_inc_ss+y_spouse_inc_ss);
    % Wage Share of Income
    yshr_wage = (y_earn_ss+y_spouse_inc_ss)./(y_inc_ss+y_spouse_inc_ss);
    % SS Share of Income
    yshr_SS = (SS_ss)./(y_inc_ss+y_spouse_inc_ss);
    % Tax Rate
    yshr_tax = (tax_ss)./(y_inc_ss+y_spouse_inc_ss);
    % Net Transfer Share
    yshr_nettxss = (tax_ss - SS_ss)./(y_inc_ss+y_spouse_inc_ss);

    mp_dsvfi_results('yshr_interest_ss') = yshr_interest;
    mp_dsvfi_results('yshr_wage_ss') = yshr_wage;
    mp_dsvfi_results('yshr_SS_ss') = yshr_SS;
    mp_dsvfi_results('yshr_tax_ss') = yshr_tax;
    mp_dsvfi_results('yshr_nttxss_ss') = yshr_nettxss;
    
end

%% Results Basic Print and Verbose print
if (bl_print_ds_verbose)
    ff_container_map_display(mp_params);
    ff_container_map_display(mp_controls);
    ff_container_map_display(mp_dsvfi_results);
end

if (bl_compute_drv_stats)
    % Array Inputs
    mp_cl_ar_xyz_of_s = containers.Map('KeyType','char', 'ValueType','any');
    mp_cl_ar_xyz_of_s('ap_ss') = {ap_ss(:), zeros(1)};
    mp_cl_ar_xyz_of_s('cons_ss') = {cons_ss(:), zeros(1)};
    
    if (bl_ds_store_all)
        
        mp_cl_ar_xyz_of_s('a_ss') = {a_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('v_ss') = {v_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('n_ss') = {n_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_head_inc') = {y_inc_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_head_earn') = {y_earn_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_spouse_inc') = {y_spouse_inc_ss(:), zeros(1)};
        mp_cl_ar_xyz_of_s('y_all') = {y_inc_ss(:) + y_spouse_inc_ss(:), zeros(1)};
        
        % Add to Map
        mp_cl_ar_xyz_of_s('yshr_interest') = {yshr_interest(:), zeros(1)};
        mp_cl_ar_xyz_of_s('yshr_wage') = {yshr_wage(:), zeros(1)};
        mp_cl_ar_xyz_of_s('yshr_SS') = {yshr_SS(:), zeros(1)};
        mp_cl_ar_xyz_of_s('yshr_tax') = {yshr_tax(:), zeros(1)};
        mp_cl_ar_xyz_of_s('yshr_nttxss') = {yshr_nettxss(:), zeros(1)};        
        
        % Add Names to list
        mp_cl_ar_xyz_of_s('ar_st_y_name') = ...
            ["a_ss", "ap_ss", "cons_ss", "v_ss", "n_ss", ...
            "y_all", "y_head_inc", "y_head_earn", "y_spouse_inc", ...
            "yshr_interest", "yshr_wage", "yshr_SS", ...
            "yshr_tax", "yshr_nttxss"];
        
    else
        mp_cl_ar_xyz_of_s('ar_st_y_name') = ["ap_ss", "c_ss"];
    end

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
    if (bl_ds_store_all)
        mp_dsvfi_results('mp_cl_mt_xyz_of_s') = mp_cl_mt_xyz_of_s;
    end
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = Phi_true;
    elseif (it_k==2)
        ob_out_cur = Phi_adj;
    elseif (it_k==3)
        ob_out_cur = A_agg;
    elseif (it_k==4)
        ob_out_cur = Y_inc_agg;
    elseif (it_k==5)
        ob_out_cur = it;        
    elseif (it_k==6)
        ob_out_cur = mp_dsvfi_results;        
    end
    varargout{it_k} = ob_out_cur;
end

end
