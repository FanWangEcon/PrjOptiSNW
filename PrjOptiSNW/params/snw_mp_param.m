%% SNW_MP_PARAM Organizes and Sets Various Model Input Scalar and Array Parameters
%    SNW_MP_PARAM sets several main default parameter structures. There are
%    four types of parameters: scalar (pref, tech, prices), state space
%    arrays, Transition Matrixes (Exogenous), Permanent Type (Life-cycle
%    Arrays)
%
%    ST_PARAM_GROUP options:
%
%    * "default_base": Main
%    * "default_small": Quick Testing
%    * "default_tiny": Quicker Testing
%
%    Pref, Technology, and prices SCALARS:
%
%    * BETA discount
%    * THETA total factor productivity normalizer
%    * R interest rate
%
%    Vectorized State Space ARRAYS:
%
%    * AGRID asset grid
%    * ETA_GRID productivity shock grid
%
%    Transition Matrixes ARRAYS:
%
%    * PI_ETA shock productivity transition
%    * PI_KIDS shock kids count transition
%    * PSI shock survival probability
%
%    Permanent Education Type Heterogeneity ARRAYS:
%
%    * EPSILON perfect-foresight education type transition
%    * SS Social Security
%
%    MP_PARAMS = SNW_MP_PARAM() get default parameters all in the same
%    container map
%
%    MP_PARAMS = SNW_MP_PARAM(ST_PARAM_GROUP) generates default parameters
%    for the type ST_PARAM_GROUP. ST_PARAM_GROUP groups include:
%    "default_base", "default_small", etc.
%
%    MP_PARAMS = SNW_MP_PARAM(ST_PARAM_GROUP, BL_PRINT_MP_PARAMS) generates
%    default parameters for the type ST_PARAM_GROUP, display parameter map
%    details if BL_PRINT_MP_PARAMS is true.
%
%    MP_PARAMS = SNW_MP_PARAM(ST_PARAM_GROUP, BL_PRINT_MP_PARAMS,
%    IT_ROW_N_KEEP, IT_COL_N_KEEP) Control for output matrixes how many
%    rows and columns to print out.
%
%    [MP_PARAMS, MP_PARAMS_PREFTECHPRICEGOV, MP_PARAMS_STATESGRID,
%    MP_PARAMS_EXOTRANS, MP_PARAMS_TYPELIFE, MP_PARAMS_INTLEN] =
%    SNW_MP_PARAM(ST_PARAM_GROUP) generates default parameters for the type
%    ST_PARAM_GROUP, and output all parameters in one MP_PARAMS, but also
%    parameters oraganized in submaps by parameter types
%
%    See also SNWX_MP_PARAM
%

%%
function varargout = snw_mp_param(varargin)
%% Parse Main Inputs and Set Defaults
if (~isempty(varargin))

    [it_row_n_keep, it_col_n_keep] = deal(8, 8);

    if (length(varargin)==1)
        st_param_group = varargin{:};
        bl_print_mp_params = false;
    elseif (length(varargin)==2)
        [st_param_group, bl_print_mp_params] = varargin{:};
    elseif (length(varargin)==4)
        [st_param_group, bl_print_mp_params, it_row_n_keep, it_col_n_keep] = varargin{:};
    end

else

    st_param_group = 'default_base';
    st_param_group = 'default_tiny';
    bl_print_mp_params = true;
    [it_row_n_keep, it_col_n_keep] = deal(8, 8);

end

%% Parametesr Grid Points
% Number of grid points
n_educgrid=2; % No. of grid points for educational attainment (college vs. non-college)
n_marriedgrid=2; % No. of grid points for marital status
n_kidsgrid=6; % No. of grid points for children (0 to 5+ children)
if(strcmp(st_param_group, "default_base"))
    n_jgrid  =42; % Age runs from 18 to 100 (a period is 2 years)
    it_retire=25;
    n_agrid  =40; % No. of grid points for assets
    n_etagrid =7; % No. of grid points for persistent labor productivity shocks
    n_kidsgrid=6; % No. of grid points for children (0 to 5+ children)
elseif(strcmp(st_param_group, "default_tiny"))
    n_jgrid   =7; % Age runs from 18 to 100 (5 periods of 16 years + terminal)
    it_retire =5;
    n_agrid   =10; % No. of grid points for assets
    n_etagrid =5; % No. of grid points for persistent labor productivity shocks
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
elseif(strcmp(st_param_group, "default_small"))
    n_jgrid   =18; % Age runs from 18 to 100 (16 periods of 5 years + terminal)
    it_retire =13;
    n_agrid   =20; % No. of grid points for assets
    n_etagrid =5; % No. of grid points for persistent labor productivity shocks
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
end

%% Preferences, Technologies, etc.
% Non-calibrated parameters
if(strcmp(st_param_group, "default_base"))
    gamma=1; % Risk aversion parameter, log kids don't matter
else
    gamma=2; % Risk aversion parameter, change to 2 so kids matter
end
rho_eta=0.98; % Ozkan (2017)
sigma_eta=0.02; % Ozkan (2017)

if(strcmp(st_param_group, "default_base"))
    g_n=(1.01^2)-1; % Annual population growth of 1.1 percent
    r=(1.04^2)-1; % Annual real interest rate of 4.0 percent from McGrattan and Prescott
else
    g_n=(1.01^(80/(n_jgrid-2)))-1; 
    r=(1.04^(80/(n_jgrid-2)))-1;
end

% Government budget constraint parameters
g_cons=0.17575574; % Government consumption expenditures to GDP (BEA: Average 2015-2019)
a2=3.664; % Initial guess for average income tax burden (if we use GS)

% Calibrated parameters
beta=0.96077; % Discount factor
theta=0.42315; % TFP parameter to normalize units such that GDP per capita equals 1 in the model: Real GDP/capita in 2019: $58,056

% Consumption allocation rule (1=uniform; 2=square root)
cons_allocation_rule=2;

% Social Security benefits
SS=zeros(n_jgrid,2);
SS(it_retire:end,1)=0.24433; % Average SS non-college 2005-2009 as a share of GDP per capita
SS(it_retire:end,2)=0.29263; % Average SS college 2005-2009 as a share of GDP per capita

%% PARAM Mortality
% Assume MORT_PROB will have 80 rows. Average based on Resulting Dataframes

load('Mortality_prob_by_age_20_99.mat','mort_prob') % Age-specific mortality probabilities (20-99 year-olds)
psi_aux=1-mort_prob;

% Convert to two-year survival probabilities
if(strcmp(st_param_group, "default_base"))
    % A3. File is 81 by 2, uses the 1st to the 80th row, avg every two
    bl_loop_trans = false;
    if (bl_loop_trans)
        psi=NaN(40,1);
        for i=1:40
            psi(i)=prod(psi_aux((2*i-1):2*i));
        end
    else
        %This is kept here to make sure it generates the same results as
        %the original looped code, below tiny and small use the same
        %structure for averaging.
        psi = prod(reshape(psi_aux, 2, []), 1)';
    end
else
    psi = prod(reshape(psi_aux, 80/(n_jgrid-2), []), 1)';
end

% Averaging

psi=[psi(1);psi]; % Let survival probability of 18-year-olds be the same as that for 20-year-olds
psi=[psi;0]; % Maximum lifespan=100 (survival probability at age 100=0)

clear mort_prob psi_aux

%% PARAM Productivity by Education Types
% Generate epsilon matrix
% A1. load external
load('Life_cycle_prod_by_educ.mat','life_cycle_prod_by_educ') % Life-cycle labor productivity for 20-100 year-olds by education (non-college vs. college)
epsilon_aux=life_cycle_prod_by_educ;
% A2. Initialize
epsilon=NaN(n_jgrid-1,2);
if(strcmp(st_param_group, "default_base"))
    % A3. File is 81 by 2, uses the 1st to the 80th row, avg every two
    bl_loop_trans = false;
    if (bl_loop_trans)
        for i=1:40
            for e=1:2
                epsilon(i,e)=sum(epsilon_aux((2*i-1):2*i,e))/2;
            end
        end
        epsilon(end,:)=epsilon_aux(end,:);
    else
        %This is kept here to make sure it generates the same results as
        %the original looped code, below tiny and small use the same
        %structure for averaging.
        for e=1:2
            epsilon(1:(end-1),e) = mean(reshape(epsilon_aux(1:80, e), 2, 40))';
        end
        epsilon(end, :) = epsilon_aux(end,:);
    end
else
    for e=1:2
        epsilon(1:(end-1),e) = mean(reshape(epsilon_aux(1:80, e), 80/(n_jgrid-2), (n_jgrid-2)))';
    end
    epsilon(end, :) = epsilon_aux(end,:);
end
% A5. init ages
epsilon=[epsilon(1,:);epsilon]; % Let life-cycle labor productivity of 18-year-olds be the same as that for 20-year-olds
% A6. Older
epsilon(it_retire:end,:)=0; % Assume zero labor productivity for 65+ year-olds (exogenous retirement)

clear life_cycle_prod_by_educ epsilon_aux

%% PARAM Kids Transition
% Transition probabilities for number of children (0, 1, 2, 3, 4, or 5) (stored in the following order:
% Number of children in year 1, age, educational status, marital status. Each column refers to the number of children
% in year 2)
%
% PI_KIDS_TRANS_PROB has 6*2*41 = 492 rows. It is assuming
% strcmp(st_param_group, "default_base"). For tiny and small grid, select
% appropriate subset of rows to match appropriate age groups from PI_KIDS.
%
% pi_kids_trans_prob estimation might change in the future if there are
% more age groups added. To make sure this code here works, use
% N_JGRID_PKTP as the N_JGRID count for PI_KIDS_TRANS_PROB. But assume
% there will always be up to 6 kids, and 2 marriage/notmarry groups.

load('pi_kids_trans_prob','pi_kids_trans_prob')
n_kidsgrid_kptp = 6;
n_marriedgrid_kptp = 2;
n_jgrid_pktp = size(pi_kids_trans_prob, 1)/(n_kidsgrid_kptp*n_marriedgrid_kptp);
pi_kids_pktp =NaN(...
    n_kidsgrid_kptp,...
    n_kidsgrid_kptp,...
    n_jgrid_pktp,...
    n_marriedgrid_kptp);

for kidsp=1:n_kidsgrid_kptp % No. of kids in year 2
    counter=0;
    for kids=1:n_kidsgrid_kptp % No. of kids in year 1
        for j=1:(n_jgrid_pktp-1) % Age in year 1
            for married=1:n_marriedgrid_kptp % Marital status
                counter=counter+1;
                pi_kids_pktp(kids,kidsp,j,married)=pi_kids_trans_prob(counter,kidsp);
            end
        end
    end
end

pi_kids_pktp(:,:,n_jgrid_pktp,:)=0;
pi_kids_pktp(:,1,n_jgrid_pktp,:)=1;

% Subset Selection
if(strcmp(st_param_group, "default_base"))
    pi_kids = pi_kids_pktp;
else
    pi_kids = pi_kids_pktp(...
        1:n_kidsgrid, ...
        1:n_kidsgrid, ...
        round(linspace(1, n_jgrid_pktp, n_jgrid)), ...
        : ...
        );
end

% Ensure that all rows sum to 1 in case of rounding error
for kids=1:n_kidsgrid % No. of kids in year 1
    for j=1:(n_jgrid-1) % Age in year 1
        for married=1:n_marriedgrid % Marital status
            aux_sum=sum(pi_kids(kids,:,j,married));
            pi_kids(kids,:,j,married)=pi_kids(kids,:,j,married)/aux_sum;
        end
    end
end

%

clear aux_sum counter pi_kids_trans_prob pi_kids_pktp

%% PARAM Specifying asset grid (use non-linear spacing with minimum value of 0)
curv=3; % Governs how the grid points are allocated
scale_a=50; % Maximum value of assets (NOTE: Verifying that it does not bind in Aggregation.m)

agrid=zeros(n_agrid,1);

for i=2:n_agrid
	agrid(i)=scale_a*((i-1)/(n_agrid-1))^curv;
end

%% PARAM Derive transition probabilities and stationary distribution for productivity shock

% Discretize process for persistent productivity shocks and derive stationary distribution
[eta_grid,pi_eta]=rouwenhorst(rho_eta,sqrt(sigma_eta),n_etagrid);

stat_distr_eta=NaN(1,n_etagrid);
x0=(1/n_etagrid)*ones(1,n_etagrid);
err=1;
tol=10^-12;
while err>tol
    x1=x0*pi_eta(:,:);
    err=max(abs(x1-x0));
    if err>tol
       x0=x1; 
    end
end
stat_distr_eta(1,:)=x0;

%% Initial conditions for marital status, college attainment, and number of kids
% Distribution of educational attainment from PSID
% tab Rcollege if RAGE>=18 & RAGE!=. [aweight=WEIGHT]
stat_distr_educ(1,1)=0.6970; % No college
stat_distr_educ(1,2)=0.3030; % College

% Distribution of marital status from PSID
% tab Rmarried if RAGE>=18 & RAGE!=. [aweight=WEIGHT]
stat_distr_married(1,1)=0.5242; % Not married
stat_distr_married(1,2)=0.4758; % Married

% Stationary distribution of children at age 20 from PSID
% Not married and no college
% tab kids if Rmarried==0 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn6(1,1,1)=0.7333;
stat_distr_kids_kn6(1,1,2)=0.1513;
stat_distr_kids_kn6(1,1,3)=0.0828;
stat_distr_kids_kn6(1,1,4)=0.0236;
stat_distr_kids_kn6(1,1,5)=0.0059;
stat_distr_kids_kn6(1,1,6)=0.0030;

% Not married but college-educated
% tab kids if Rmarried==0 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn6(2,1,1)=0.9752;
stat_distr_kids_kn6(2,1,2)=0.0236;
stat_distr_kids_kn6(2,1,3)=0.0001;
stat_distr_kids_kn6(2,1,4)=0.0011;
stat_distr_kids_kn6(2,1,5)=0;
stat_distr_kids_kn6(2,1,6)=0;

% Married and no college
% tab kids if Rmarried==1 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn6(1,2,1)=0.4143;
stat_distr_kids_kn6(1,2,2)=0.2958;
stat_distr_kids_kn6(1,2,3)=0.2131;
stat_distr_kids_kn6(1,2,4)=0.0569;
stat_distr_kids_kn6(1,2,5)=0.0199;
stat_distr_kids_kn6(1,2,6)=0;


% Married and college-educated
% tab kids if Rmarried==1 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn6(2,2,1)=0.7534;
stat_distr_kids_kn6(2,2,2)=0.2153;
stat_distr_kids_kn6(2,2,3)=0.0221;
stat_distr_kids_kn6(2,2,4)=0.0091;
stat_distr_kids_kn6(2,2,5)=0;
stat_distr_kids_kn6(2,2,6)=0;

% Subset Kids for Smaller Solutions
if(strcmp(st_param_group, "default_base"))
    stat_distr_kids = stat_distr_kids_kn6;
else
    stat_distr_kids = stat_distr_kids_kn6(:,:,1:n_kidsgrid);
end

% Reweight
bl_sum_onebyone = true;
if (bl_sum_onebyone)
    aux=sum(stat_distr_kids(1,1,:));
    stat_distr_kids(1,1,:)=stat_distr_kids(1,1,:)/aux;
    aux=sum(stat_distr_kids(2,1,:));
    stat_distr_kids(2,1,:)=stat_distr_kids(2,1,:)/aux;
    aux=sum(stat_distr_kids(1,2,:));
    stat_distr_kids(1,2,:)=stat_distr_kids(1,2,:)/aux;
    aux=sum(stat_distr_kids(2,2,:));
    stat_distr_kids(2,2,:)=stat_distr_kids(2,2,:)/aux;
else
mt_aux = repmat(sum(stat_distr_kids, 3), ...
    [1,1,size(stat_distr_kids,3)]);
    stat_distr_kids(2,2,:)=stat_distr_kids./mt_aux;
end

clear aux

%% Population distribution
% Normalize mass of 18-year-olds to 1
Pop=zeros(n_jgrid,1);
Pop(1)=1;
for j=2:n_jgrid
    Pop(j)=Pop(j-1)*psi(j-1)/(1+g_n);
end

name='Old-age dependency ratio=';
st_old_age_depend =[name,num2str(sum(Pop(it_retire:end))/sum(Pop(1:(it_retire-1))))];
% disp(st_old_age_depend);

%% Set Parameter Maps
mp_params_preftechpricegov = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_preftechpricegov('gamma') = gamma;
mp_params_preftechpricegov('beta') = beta;
mp_params_preftechpricegov('theta') = theta;
mp_params_preftechpricegov('cons_allocation_rule') = cons_allocation_rule;
mp_params_preftechpricegov('r') = r;
mp_params_preftechpricegov('g_n') = g_n;
mp_params_preftechpricegov('g_cons') = g_cons;
mp_params_preftechpricegov('a2') = a2;

mp_params_statesgrid = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_statesgrid('agrid') = agrid;
mp_params_statesgrid('eta_grid') = eta_grid;

mp_params_exotrans = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_exotrans('pi_eta') = pi_eta;
mp_params_exotrans('pi_kids') = pi_kids;
mp_params_exotrans('psi') = psi;

mp_params_typelife = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_typelife('epsilon') = epsilon;
mp_params_typelife('SS') = SS;

mp_params_intlen = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_intlen('n_jgrid') = n_jgrid;
mp_params_intlen('n_agrid') = n_agrid;
mp_params_intlen('n_etagrid') = n_etagrid;
mp_params_intlen('n_educgrid') = n_educgrid;
mp_params_intlen('n_marriedgrid') = n_marriedgrid;
mp_params_intlen('n_kidsgrid') = n_kidsgrid;

mp_params_stat = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_stat('stat_distr_eta') = stat_distr_eta;
mp_params_stat('stat_distr_educ') = stat_distr_educ;
mp_params_stat('stat_distr_married') = stat_distr_married;
mp_params_stat('stat_distr_kids') = stat_distr_kids;
mp_params_stat('Pop') = Pop;
mp_params_stat('st_old_age_depend') = string(st_old_age_depend);

%% Combine Maps
mp_params = [mp_params_preftechpricegov; mp_params_statesgrid ; ...
    mp_params_exotrans; mp_params_typelife; mp_params_intlen ; ...
    mp_params_stat];

mp_params('mp_params_name') = string(st_param_group);

% MP_PARAMS = [MP_PARAMS_PREFTECHPRICE; MP_PARAMS_STATESGRID ; ...
%     MP_PARAMS_EXOTRANS; MP_PARAMS_TYPELIFE; MP_PARAMS_INTLEN];

%% Print
if (bl_print_mp_params)
    ff_container_map_display(mp_params_preftechpricegov);
    ff_container_map_display(mp_params_intlen);
    ff_container_map_display(mp_params_statesgrid, it_row_n_keep, it_col_n_keep);
    ff_container_map_display(mp_params_exotrans, it_row_n_keep, it_col_n_keep);
    ff_container_map_display(mp_params_typelife, it_row_n_keep, it_col_n_keep);
    ff_container_map_display(mp_params_stat, it_row_n_keep, it_col_n_keep);
end

%% Return
if (nargout==1)
    varargout = cell(nargout,0);
    varargout{1} = mp_params;
elseif (nargout==6)
    varargout = cell(nargout,0);
    varargout{1} = mp_params;
    varargout{2} = mp_params_preftechpricegov;
    varargout{3} = mp_params_statesgrid;
    varargout{4} = mp_params_exotrans;
    varargout{5} = mp_params_typelife;
    varargout{6} = mp_params_intlen;
end

end
