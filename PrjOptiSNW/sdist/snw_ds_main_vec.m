%% SNW_DS_MAIN_VEC Simulates Distributions Exact Savings Vectorized
%    Given policy functions, simulate distributions. When VFI results are
%    not provided, will solve for VFI using SNW_VFI_MAIN_BISEC_VEC.
%
%    This is the vectorized version of SNW_DS_MAIN with significant speed
%    gains.
%
%    [Phi_true,Phi_adj,A_agg,Y_inc_agg,it,mp_dsvfi_results] =
%    SNW_DS_MAIN_VEC() invoke model with externally set parameter map and
%    control map. Results outputed to a map containing various output
%    matrixes in mp_dsvfi_results, and also distributional matrixes.
%
%    See also SNW_DS_MAIN, SNW_DS_GRID_SEARCH
%

%%
function varargout=snw_ds_main_vec(varargin)

%% Default and Parse Inputs
if (~isempty(varargin))

    bl_loop_aggregate = true;

    if (length(varargin)==2)
        [mp_params, mp_controls] = varargin{:};
    elseif (length(varargin)==4)
        % This will not produce extra statistics outputs
        [mp_params, mp_controls, ap_ss, cons_ss] = varargin{:};
    elseif (length(varargin)==5)
        % This will produce extra statistics outputs
        [mp_params, mp_controls, ap_ss, cons_ss, mp_valpol_more_ss] = varargin{:};
    else
        error('Need to provide 2/4 parameter inputs');
    end

    bl_ds_store_all = false;

else

%     clc;
    clear all;
%     mp_params = snw_mp_param('default_docdense');
    mp_params = snw_mp_param('default_small');
%     mp_params = snw_mp_param('default_tiny');
    mp_controls = snw_mp_control('default_test');
    [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);

    bl_ds_store_all = true;
    bl_loop_aggregate = true;

end

%% Solve VFI if not Provided

% output more
if (nargout >= 6)
    bl_ds_store_all = true;
end

% requesting all stats outputs, but does not have vfi more info
if (bl_ds_store_all && ~exist('mp_valpol_more_ss','var'))
    [v_ss, ap_ss, cons_ss, mp_valpol_more_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
end

% requesting subset outputs if don't have policy function, and does not
% request extra outputs
if (~exist('ap_ss','var'))
    [v_ss, ap_ss, cons_ss] = snw_vfi_main_bisec_vec(mp_params, mp_controls);
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
% global a2 g_cons agrid SS pi_eta pi_kids Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Added
% global epsilon theta eta_H_grid eta_S_grid r g_n psi Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r', 'g_n', 'g_cons', 'a2','jret'});
[theta, r, g_n, g_cons, a2,jret] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, ...
    {'pi_eta', 'pi_kids', 'cl_mt_pi_jem_kidseta', 'psi'});
[pi_eta, pi_kids, cl_mt_pi_jem_kidseta, psi] = params_group{:};

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

%% Timer

if (bl_timer)
    tm_start = tic;
end

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

%% Vectorized Life-Cycle Probability Mass Calculations
% Vectorized method to compute probability mass: see
% https://fanwangecon.github.io/MEconTools/MEconTools/doc/ds/htmlpdfm/fx_ds_az_cts_vec.html
% for a simpler example of the algorithm below

% A. initialize age/shocks sub-matrix
mt_dist_az_zeros = zeros(n_agrid,n_etagrid*n_kidsgrid);

% B1. Loop over Age
for j=1:(n_jgrid-1) % Age

    % Age Timer
    if (bl_print_ds) tm_ds_age = tic; end

    % B2. Loop over Permanent Heterogeneity Education Levels
    for educ=1:n_educgrid

        % B3. Loop over Permanent Heterogeneity Marital Status
        for married=1:n_marriedgrid

            % C1. Get P(S'|S), S = [eta x kids] by [eta x kids] transition matrix
            % transition matrix is specific to jem:age/edu/marry
            mt_pi_jem_kidseta = kron(pi_kids(:,:,j,educ,married), pi_eta);

            % C2. Choice and Current/Last Age Probability Mass 6D to 2D
            % transformed matrix rows are savings, columns are shocks
            % column shocks include productivity as well kids shocks
            mt_ap_ss_jem = reshape(ap_ss(j,:,:,educ,married,:), n_agrid, []);
            mt_Phiss_jem = reshape(Phiss(j,:,:,educ,married,:), n_agrid, []);

            % C3a. lower/upper ratio and lower index for savings choices
            % These determine where do mass from age j land at at j+1
            % Rows are savings states, columns are shocks for mt
            ar_it_ap_near_lower_idx  = sum(agrid' <= mt_ap_ss_jem(:), 2);
            mt_it_ap_near_lower_idx = reshape(ar_it_ap_near_lower_idx, size(mt_ap_ss_jem));
            % find higher index
            mt_it_ap_near_higher_idx  = mt_it_ap_near_lower_idx + 1;
            mt_it_ap_near_higher_idx(mt_it_ap_near_higher_idx >= n_agrid) = n_agrid;
            % reshape to apindex matrix
            mt_fl_ap_near_lower = reshape(agrid(mt_it_ap_near_lower_idx), size(mt_ap_ss_jem));
            mt_fl_ap_near_higher = reshape(agrid(mt_it_ap_near_higher_idx), size(mt_ap_ss_jem));

            % C3b. Chance of landing oat just above and just below points
            mt_fl_lower_idx_share = 1 - (mt_ap_ss_jem - mt_fl_ap_near_lower)./(mt_fl_ap_near_higher - mt_fl_ap_near_lower);
            mt_fl_lower_idx_share(mt_it_ap_near_higher_idx == mt_it_ap_near_lower_idx) = 1;

            % C4. Initialize jem specific az distributional mass
            % initialize empty
            mt_dist_az = mt_dist_az_zeros;

            % C5. Loop over all possible shock states
            % household head x spousal x kids shocks
            for it_z_i=1:size(mt_pi_jem_kidseta,2)

                % C6a. current probablity at f(:,z) for all a
                ar_cur_za_prob = mt_Phiss_jem(:, it_z_i);

                % C6b. only consider positive mass elements
                ar_it_cur_za_prob_nonzero = find(ar_cur_za_prob);
                ar_cur_za_prob = ar_cur_za_prob(ar_it_cur_za_prob_nonzero);

                % C6c. Proportion of Current Mass to send to the left closest
                ar_aprime_lower_share = mt_fl_lower_idx_share(ar_it_cur_za_prob_nonzero, it_z_i);

                % C6d. f(z'|z) transition for all z' given current z
                ar_ztoz_trans = mt_pi_jem_kidseta(it_z_i, :);

                % C6e. P(a',z'|a,z)*P(a,z) for where P(a,z)>0, given z
                mt_zfromza_lower = (ar_aprime_lower_share.*ar_cur_za_prob)*ar_ztoz_trans;
                mt_zfromza_higher = ((1-ar_aprime_lower_share).*ar_cur_za_prob)*ar_ztoz_trans;

                % C7. Loop over positive P(a,z)>0 to add mass to (a',z')
                it_pos_counter = 0;
                for it_a_j = ar_it_cur_za_prob_nonzero'

                    % C8a. Counter because mt_zfromza_lower index by
                    % nonzero P(a,z)>0, no usint it_a_j index.
                    it_pos_counter = it_pos_counter + 1;

                    % C8b. Use it_a_j index, which are full all a index
                    it_aprime_lower = mt_it_ap_near_lower_idx(it_a_j, it_z_i);
                    it_aprime_higher = mt_it_ap_near_higher_idx(it_a_j, it_z_i);

                    % C8c. assign mass to (a',z') for all z' at once
                    mt_dist_az(it_aprime_lower, :) = mt_dist_az(it_aprime_lower, :) + mt_zfromza_lower(it_pos_counter, :);
                    mt_dist_az(it_aprime_higher, :) = mt_dist_az(it_aprime_higher, :) + mt_zfromza_higher(it_pos_counter, :);
                end

            end

            % D. Assign Mass to overall Probability Mass 6D full Structure
            mn_dist_az = reshape(mt_dist_az, [n_agrid, n_etagrid, n_kidsgrid]);
            Phiss(j+1,:, :, educ, married, :) = mn_dist_az;

       end
   end

   if (bl_print_ds)
       tm_ds_age_end = toc(tm_ds_age);
       disp(strcat(['SNW_DS_MAIN_VEC ACUMU MASS: Finished Age Group:' ...
           num2str(j) ' of ' num2str(n_jgrid-1) ...
           ', time-this-age:' num2str(tm_ds_age_end)]));
   end

end

%% Normalize distribution of idiosyncratic states to sum to 1 for each age
Phi_adj=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
Phi_true=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid

    % Age Timer
    if (bl_print_ds) tm_ds_age = tic; end

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

   if (bl_print_ds)
       tm_ds_age_end = toc(tm_ds_age);
       disp(strcat(['SNW_DS_MAIN NORMALIZE MASS: Finished Age Group:' ...
           num2str(j) ' of ' num2str(n_jgrid-1) ...
           ', time-this-age:' num2str(tm_ds_age_end)]));
   end

end

% Check if the upper bound on assets binds
check_asset_distr=sum(sum(sum(sum(sum(Phi_true(:,n_agrid,:,:,:,:))))));
if (bl_print_ds)
    disp(strcat(['SNW_DS_MAIN: Share of population with assets equal to upper bound on asset grid:' ...
        num2str(check_asset_distr/sum(Pop))]));
end

%% Explicit Looped Aggregate variables Calculation
if (bl_loop_aggregate)

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

    tol=10^-4; %10^-3; %5*10^-4;
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

        err=abs((Tax_revenues_aux/(SS_spend+g_cons*Y_inc_agg))-1);

        if (bl_print_ds)
            st_tax_iter = strjoin(...
                ["SNW_DS_MAIN_VEC tax and spend", ...
                ['it=' num2str(it)], ...
                ['err=' num2str(err)] ...
                ], ";");
            disp(st_tax_iter);
        end
    end


    if (bl_print_ds)
        name='SNW_DS_MAIN_VEC: Number of a2-adjustments (for taxation) used to balance the government budget= ';
        name2=[name,num2str(it)];
        disp(name2);

        a2_update=[a2_update,a2];

        name='SNW_DS_MAIN_VEC: Old and updated value of a2=';
        name2=[name,num2str(a2_update)];
        disp(name2);

        name='SNW_DS_MAIN_VEC: Aggregates: Cons., Gov. cons., Save, Assets, Income, Bequests ';
        aggregates=[C_agg,g_cons*Y_inc_agg,A_agg,Aprime_agg,Y_inc_agg,Bequests_aux];
        name2=[name,num2str(aggregates)];
        disp(name2);

        name='SNW_DS_MAIN_VEC: Resource constraint: C_t+A_{t+1}+G_t=A_t+Y_t ';
        name2=[name,num2str([C_agg+g_cons*Y_inc_agg+Aprime_agg,A_agg+Y_inc_agg])];
        disp(name2);
    end
end

%% Timing and Profiling End
if (bl_timer)
    tm_end = toc(tm_start);
    st_complete = strjoin(...
        ["Completed SNW_DS_MAIN_VEC", ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))], ...
         ['time=' num2str(tm_end)] ...
        ], ";");
    disp(st_complete);
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

    % Additional aggregate Statistics (Duplicate Initial Looped Method)
    if (bl_loop_aggregate)
        mp_dsvfi_results('A_agg') = A_agg;
        mp_dsvfi_results('Aprime_agg') = Aprime_agg;
        mp_dsvfi_results('Y_inc_agg') = Y_inc_agg;
        mp_dsvfi_results('C_agg') = C_agg;
        mp_dsvfi_results('Tax_revenues') = Tax_revenues;
        mp_dsvfi_results('SS_spend') = SS_spend;
        mp_dsvfi_results('Bequests_aux') = Bequests_aux;

        mp_dsvfi_results('A_agg_perhh') = A_agg/sum(Pop);
        mp_dsvfi_results('Aprime_agg_perhh') = Aprime_agg/sum(Pop);
        mp_dsvfi_results('Y_inc_agg_perhh') = Y_inc_agg/sum(Pop);

        mp_dsvfi_results('C_agg_perhh') = C_agg/sum(Pop);
        mp_dsvfi_results('Tax_revenues_perhh') = Tax_revenues/sum(Pop);
        mp_dsvfi_results('SS_spend_perhh') = SS_spend/sum(Pop);
        mp_dsvfi_results('Bequests_aux_perhh') = Bequests_aux/sum(Pop);
    end

end

%% Results Basic Print and Verbose print
if (bl_print_ds_verbose)
%     ff_container_map_display(mp_params);
%     ff_container_map_display(mp_controls);
    if exist('mp_dsvfi_results', 'var')
        ff_container_map_display(mp_dsvfi_results);
    end
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
        mp_cl_ar_xyz_of_s('y_all') = {y_inc_ss(:) + y_spouse_inc_ss(:) - SS_ss(:), zeros(1)};

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
        mp_cl_ar_xyz_of_s('ar_st_y_name') = ["ap_ss", "cons_ss"];
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
    elseif (it_k==7)
        % update the tax parameter
        ob_out_cur = a2;
    end
    varargout{it_k} = ob_out_cur;
end

end
