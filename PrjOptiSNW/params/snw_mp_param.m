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

%     st_param_group = 'default_base';
    st_param_group = 'default_verydense';
%     st_param_group = 'default_dense';
%     st_param_group = 'default_tiny';
    bl_print_mp_params = true;
    [it_row_n_keep, it_col_n_keep] = deal(20, 8);

end

%% Parametesr Grid Points
% Number of grid points
n_educgrid=2; % No. of grid points for educational attainment (college vs. non-college)
n_marriedgrid=2; % No. of grid points for marital status
n_kidsgrid=5; % No. of grid points for children (0 to 5+ children)
if(strcmp(st_param_group, "default_verydense"))
    n_jgrid  =83; % Age runs from 18 to 100 (a period is 2 years)
    jret     =48;
    n_agrid  =251; % No. of grid points for assets
    n_eta_H_grid=9; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=5; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=5; % No. of grid points for children (0 to 4+ children)
elseif(strcmp(st_param_group, "default_moredense"))
    n_jgrid  =83; % Age runs from 18 to 100 (a period is 2 years)
    jret     =48;
    n_agrid  =151; % No. of grid points for assets
    n_eta_H_grid=9; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=5; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=5; % No. of grid points for children (0 to 4+ children)
elseif(strcmp(st_param_group, "default_dense"))
    n_jgrid  =83; % Age runs from 18 to 100 (a period is 2 years)
    jret     =48;
    n_agrid  =55; % No. of grid points for assets
    n_eta_H_grid=7; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=3; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=5; % No. of grid points for children (0 to 4+ children)
elseif(strcmp(st_param_group, "default_base"))
    n_jgrid  =42; % Age runs from 18 to 100 (a period is 2 years)
    jret=25;
    n_agrid  =40; % No. of grid points for assets
    n_eta_H_grid=5; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=3; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=5; % No. of grid points for children (0 to 4+ children)
elseif(strcmp(st_param_group, "default_small53"))
    n_jgrid   =18; % Age runs from 18 to 100 (16 periods of 5 years + terminal)
    jret      =13;
    n_agrid   =25;
    n_eta_H_grid=5;
    n_eta_S_grid=3;
    n_kidsgrid=3;
elseif(strcmp(st_param_group, "default_smalla151"))
    n_jgrid   =18; % Age runs from 18 to 100 (16 periods of 5 years + terminal)
    jret      =13;
    n_agrid   =151; % No. of grid points for assets
    n_eta_H_grid=5; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=1; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
elseif(strcmp(st_param_group, "default_small"))
    n_jgrid   =18; % Age runs from 18 to 100 (16 periods of 5 years + terminal)
    jret      =13;
    n_agrid   =25; % No. of grid points for assets
    n_eta_H_grid=5; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=1; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
elseif(strcmp(st_param_group, "default_tiny53"))
    n_jgrid   =7; % Age runs from 18 to 100 (5 periods of 16 years + terminal)
    jret =5;
    n_agrid   =10; % No. of grid points for assets
    n_eta_H_grid=5; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=3; % 3; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
elseif(strcmp(st_param_group, "default_tiny"))
    n_jgrid   =7; % Age runs from 18 to 100 (5 periods of 16 years + terminal)
    jret =5;
    n_agrid   =10; % No. of grid points for assets
    n_eta_H_grid=5; % 9; % No. of grid points for persistent labor productivity shocks
    n_eta_S_grid=1; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to no spousal shocks)
    n_kidsgrid=3; % No. of grid points for children (0 to 5+ children)
end
n_etagrid=n_eta_H_grid*n_eta_S_grid; % Length of productivity shock grid

% Years Per Period
if (n_jgrid == 83)
    it_yrs_per_period = 1;
else
    it_yrs_per_period = (80/(n_jgrid-2));
end

%% Planning and Unemployment Parameters
% ADDITIONAL PARAMETERS/VARIABLES for planner problem are given in the following sections:
% "Compute value of employment and unemployment in 2020 conditional on number of welfare checks";
% "Probability of unemployment"
% "Solve planner's problem"
% xi: Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
% b: Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
% TR: Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
% n_welfchecksgrid: Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars
% pi_unemp: Probability of unemployment
% n_incgrid: Number of income groups
% inc_grid: Grid for income groups

% omega=0.0135; % Total spending on welfare checks as a share of aggregate income (estimated to cost $290 billion, or 1.35 percent of GDP in 2019. IRS has paid out $218 billion as of May 11)

xi=0.75; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
scaleconvertor = 58056;
TR=100/scaleconvertor; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
n_welfchecksgrid=51; % Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars

% Probability of unemployment
% pi_j=[0.22;0.175;0.16;0.165;0.22]; % Probability of unemployment in 2020 by age groups from Cajner et al. (2020, NBER)
% pi_w=[0.360;0.22;0.17;0.14;0.09]; % Probability of unemployment in 2020 by wage quintiles from Cajner et al. (2020, NBER)

% Columns are wage groups; rows are age groups.
pi_unemp=zeros(n_jgrid,5);
aux_mat=[0.08027790 0.05170647 0.04150239 0.03537994 0.02517586];
pi_unemp(1:13,:)=repmat(aux_mat,13,1);
aux_mat=[0.07070343 0.04213200 0.03192792 0.02580547 0.01560139];
pi_unemp(14:23,:)=repmat(aux_mat,10,1);
aux_mat=[0.06751194 0.03894051 0.02873643 0.02261398 0.01240990];
pi_unemp(24:33,:)=repmat(aux_mat,10,1);
aux_mat=[0.06857577 0.04000434 0.02980026 0.02367781 0.01347373];
pi_unemp(34:43,:)=repmat(aux_mat,10,1);
aux_mat=[0.08027790 0.05170647 0.04150239 0.03537994 0.02517586];
pi_unemp(44:48,:)=repmat(aux_mat,5,1);

clear aux_mat

% Solve planner's problem
n_incgrid=201; % Number of income groups
inc_grid=linspace(0,7,n_incgrid)'; % 7 refers to 7*58056=406392 dollars in 2012USD

%% Preferences, Technologies, etc.

% Non-calibrated parameters
gamma=2; % Risk aversion parameter
if(contains(st_param_group, "dense"))
    rho_eta=0.98; % Persistence of AR(1) productivity shocks
    sigma_eta=0.018; % Variance of AR(1) productivity shocks
    g_n=0.01; % Annual population growth of 1.1 percent
    r=0.04; % Annual real interest rate of 4.0 percent from McGrattan and Prescott
    beta=0.969010845568155; % 0.97068903873305; % Discount factor
else
    rho_eta=0.98^it_yrs_per_period;
    sigma_eta=sqrt(0.018^2*sum((0.98.^(0:(it_yrs_per_period-1))).^2));
    g_n=(1.01^it_yrs_per_period)-1;
    r=(1.04^it_yrs_per_period)-1;
    beta=0.969010845568155^it_yrs_per_period;
end

% Spousal Shocks
rho_eta_spouse=0; % Persistence of spousal AR(1) productivity shocks
sigma_eta_spouse=1.040654^2; % Variance of spousal AR(1) productivity shocks (standard deviation of residual from spousal income regression for 18-65 year-old household heads. See spousal_income.m for regression specification details)

% Bequest
% Bequests allocation rule (=1: accidental bequests go to the government; =2: accidental bequests uniformly across the population)
% Bequests=Bequests_aux/((1+g_n)*sum(sum(sum(sum(sum(sum(Phi_true(:,:,:,:,:,:))))))));
bequests_option=1;
Bequests=0.05826*(bequests_option-1);
throw_in_ocean=1; % If bequests go to the government, a value of 1 for throw_in_ocean means that all accidental bequests are "thrown in the ocean", whereas a value of 0 means the full amount goes to the government

if bequests_option==2
    a2=1.575;
elseif bequests_option==1
    if throw_in_ocean==0
        a2=0.7027;
    elseif throw_in_ocean==1
        a2=1.6081061107858;
    end
end


% Government budget constraint parameters
g_cons=0.17575574; % Government consumption expenditures to GDP (BEA: Average 2015-2019)

% Calibrated parameters
theta=0.66947637485374; % TFP parameter to normalize units such that average household income relative to GDP per capita is consistent with the data (the latter is normalized to 1): Real GDP/capita in 2019: $58,056

% Consumption allocation rule (1=uniform; 2=square root)
cons_allocation_rule=2;

% Social Security benefits
SS=zeros(n_jgrid,2);
SS(jret:end,1)=0.24433; % Average SS non-college 2005-2009 as a share of GDP per capita
SS(jret:end,2)=0.29263; % Average SS college 2005-2009 as a share of GDP per capita

%% PARAM Mortality
% Assume MORT_PROB will have 80 rows. Average based on Resulting Dataframes
load('Mortality_prob_by_age_18_99.mat','mort_prob')
psi_full=1-mort_prob;
psi_full=[psi_full;0]; % Maximum lifespan=100 (survival probability at age 100=0)

if(contains(st_param_group, "dense"))
    psi = psi_full;
else
    psi = NaN(n_jgrid,1);
    psi(1) = mean(1-mort_prob(1:2));
    psi(2:(end-1)) = prod(reshape(1-mort_prob(3:end), it_yrs_per_period, []), 1)';
    psi(end) = 0;
end

clear mort_prob psi_full

%% PARAM Productivity by Education Types
% Generate epsilon matrix
% A1. load external
load('Life_cycle_prod_by_educ.mat','life_cycle_prod_by_educ') % Life-cycle labor productivity for 20-100 year-olds by education (non-college vs. college)
% Set Annual
epsilon_full=NaN(83,2);
epsilon_full(3:end,:)=life_cycle_prod_by_educ(:,:);

epsilon_full(1,:)=epsilon_full(3,:); % Let life-cycle labor productivity of 18- and 19-year-olds be the same as that for 20-year-olds
epsilon_full(2,:)=epsilon_full(3,:);

if(contains(st_param_group, "dense"))
    epsilon = epsilon_full;
else
    epsilon=NaN(n_jgrid,2);
    epsilon(1,:) = epsilon_full(3,:);
    for e=1:2
        epsilon(2:end-1,e) = mean(reshape(life_cycle_prod_by_educ(1:80, e), it_yrs_per_period, (n_jgrid-2)))';
    end
end

epsilon(jret:end,:)=0; % Assume zero labor productivity for 65+ year-olds (exogenous retirement)

clear life_cycle_prod_by_educ epsilon_aux

%% PARAM Kids Transition
% Transition probabilities for number of children (0, 1, 2, 3, 4, or 5) (stored in the following order:
% Number of children in year 1, age, marital status, college attainment. Each column refers to
% the number of children in year 2)
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
n_kidsgrid_kptp = 5;
n_educgrid_kptp = 2;
n_marriedgrid_kptp = 2;
n_jgrid_pktp = 83;
pi_kids_pktp =NaN(...
    n_kidsgrid_kptp,...
    n_kidsgrid_kptp,...
    n_jgrid_pktp,...
    n_educgrid_kptp,...
    n_marriedgrid_kptp);

for kidsp=1:n_kidsgrid_kptp % No. of kids in year 2
    counter=0;
    for kids=1:n_kidsgrid_kptp % No. of kids in year 1
        for j=1:(n_jgrid_pktp-1) % Age in year 1
            for married=1:n_marriedgrid_kptp % Marital status
                for educ=1:n_educgrid_kptp % Educational level
                    counter=counter+1;
                    pi_kids_pktp(kids,kidsp,j,educ,married)=pi_kids_trans_prob(counter,kidsp);
                end
            end
        end
    end
end

pi_kids_pktp(:,:,n_jgrid_pktp,:,:)=0;
pi_kids_pktp(:,1,n_jgrid_pktp,:,:)=1;

% Subset Selection
if(contains(st_param_group, "dense"))
    pi_kids = pi_kids_pktp;
else
    pi_kids = pi_kids_pktp(...
        1:n_kidsgrid, ...
        1:n_kidsgrid, ...
        round(linspace(1, n_jgrid_pktp, n_jgrid)), ...
        :, : ...
        );
end

% Ensure that all rows sum to 1 in case of rounding error
for kids=1:n_kidsgrid % No. of kids in year 1
    for j=1:n_jgrid % Age in year 1
        for educ=1:n_educgrid % Educational level
            for married=1:n_marriedgrid % Marital status
                aux_sum=sum(pi_kids(kids,:,j,educ,married));
                pi_kids(kids,:,j,educ,married)=pi_kids(kids,:,j,educ,married)/aux_sum;
            end
        end
    end
end

%
clear aux_sum counter pi_kids_trans_prob pi_kids_pktp

%% PARAM Specifying asset grid (use non-linear spacing with minimum value of 0)
curv=3; % Governs how the grid points are allocated
scale_a=190; % Maximum value of assets (NOTE: Verifying that it does not bind in Aggregation.m) grid=1.2451e-11
scale_a=135; % Maximum value of assets (NOTE: Verifying that it does not bind in Aggregation.m) grid=0

agrid=zeros(n_agrid,1);

for i=2:n_agrid
	agrid(i)=scale_a*((i-1)/(n_agrid-1))^curv;
end

%% Derive transition probabilities and stationary distribution for productivity shock
% Discretize process for persistent productivity shocks and derive stationary distribution
[eta_H_grid_aux,pi_H_eta]=rouwenhorst(rho_eta,sqrt(sigma_eta),n_eta_H_grid);
[eta_S_grid_aux,pi_S_eta]=rouwenhorst(rho_eta_spouse,sqrt(sigma_eta_spouse),n_eta_S_grid);

pi_eta=NaN(n_etagrid,n_etagrid);
counter=0;
for eta_S=1:n_eta_S_grid
    for eta_H=1:n_eta_H_grid
        counter=counter+1;

        counterp=0;
        for eta_Sp=1:n_eta_S_grid
            for eta_Hp=1:n_eta_H_grid
                counterp=counterp+1;
                pi_eta(counter,counterp)=pi_H_eta(eta_H,eta_Hp)*pi_S_eta(eta_S,eta_Sp);
            end
        end
    end
end

eta_H_grid=repmat(eta_H_grid_aux,n_eta_S_grid,1);
eta_S_grid=sort(repmat(eta_S_grid_aux,n_eta_H_grid,1));

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

clear counter counterp eta_H_grid_aux eta_S_grid_aux

%% Initial conditions for marital status, college attainment, and number of kids
% Distribution of educational attainment from PSID
% tab Rcollege if RAGE>=18 & RAGE!=. [aweight=WEIGHT]
stat_distr_educ(1,1)=0.6970; % No college
stat_distr_educ(1,2)=0.3030; % College

% Distribution of marital status conditional on college attainment from PSID
% tab  Rmarried if RAGE>=18 & RAGE!=. & Rcollege==0 [aweight=WEIGHT]
stat_distr_married(1,1)=0.5635; % Not married
stat_distr_married(1,2)=0.4365; % Married

% tab  Rmarried if RAGE>=18 & RAGE!=. & Rcollege==1 [aweight=WEIGHT]
stat_distr_married(2,1)=0.4364; % Not married
stat_distr_married(2,2)=0.5636; % Married

% Stationary distribution of children at age 20 from PSID
% Not married and no college
% tab kids if Rmarried==0 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn5(1,1,1)=0.7333;
stat_distr_kids_kn5(1,1,2)=0.1513;
stat_distr_kids_kn5(1,1,3)=0.0828;
stat_distr_kids_kn5(1,1,4)=0.0236;
stat_distr_kids_kn5(1,1,5)=0.0059;

% Not married but college-educated
% tab kids if Rmarried==0 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn5(2,1,1)=0.9752;
stat_distr_kids_kn5(2,1,2)=0.0236;
stat_distr_kids_kn5(2,1,3)=0.0001;
stat_distr_kids_kn5(2,1,4)=0.0011;
stat_distr_kids_kn5(2,1,5)=0;

% Married and no college
% tab kids if Rmarried==1 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn5(1,2,1)=0.4143;
stat_distr_kids_kn5(1,2,2)=0.2958;
stat_distr_kids_kn5(1,2,3)=0.2131;
stat_distr_kids_kn5(1,2,4)=0.0569;
stat_distr_kids_kn5(1,2,5)=0.0199;

% Married and college-educated
% tab kids if Rmarried==1 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids_kn5(2,2,1)=0.7534;
stat_distr_kids_kn5(2,2,2)=0.2153;
stat_distr_kids_kn5(2,2,3)=0.0221;
stat_distr_kids_kn5(2,2,4)=0.0091;
stat_distr_kids_kn5(2,2,5)=0;

% Subset Kids for Smaller Solutions
if(strcmp(st_param_group, "default_base"))
    stat_distr_kids = stat_distr_kids_kn5;
else
    stat_distr_kids = stat_distr_kids_kn5(:,:,1:n_kidsgrid);
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

name='Old-age dependency ratio (ratio of 65+/(18-64))=';
st_old_age_depend =[name,num2str(sum(Pop(jret:end))/sum(Pop(1:(jret-1))))];
% disp(st_old_age_depend);

%% Set Parameter Maps
mp_params_covid_unemploy = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_covid_unemploy('xi') = xi;
mp_params_covid_unemploy('b') = b;
mp_params_covid_unemploy('scaleconvertor') = scaleconvertor;
mp_params_covid_unemploy('TR') = TR;
mp_params_covid_unemploy('n_welfchecksgrid') = n_welfchecksgrid;
mp_params_covid_unemploy('pi_unemp') = pi_unemp;
mp_params_covid_unemploy('n_incgrid') = n_incgrid;
mp_params_covid_unemploy('inc_grid') = inc_grid;

mp_params_preftechpricegov = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_preftechpricegov('gamma') = gamma;
mp_params_preftechpricegov('beta') = beta;
mp_params_preftechpricegov('theta') = theta;
mp_params_preftechpricegov('cons_allocation_rule') = cons_allocation_rule;
mp_params_preftechpricegov('r') = r;
mp_params_preftechpricegov('g_n') = g_n;
mp_params_preftechpricegov('g_cons') = g_cons;
mp_params_preftechpricegov('a2') = a2;
mp_params_preftechpricegov('jret') = jret;
mp_params_preftechpricegov('Bequests') = Bequests;
mp_params_preftechpricegov('bequests_option') = bequests_option;
mp_params_preftechpricegov('throw_in_ocean') = throw_in_ocean;

mp_params_statesgrid = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_statesgrid('agrid') = agrid;
mp_params_statesgrid('eta_H_grid') = eta_H_grid;
mp_params_statesgrid('eta_S_grid') = eta_S_grid;

mp_params_exotrans = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_exotrans('pi_eta') = pi_eta;
mp_params_exotrans('pi_H_eta') = pi_H_eta;
mp_params_exotrans('pi_S_eta') = pi_S_eta;
mp_params_exotrans('pi_kids') = pi_kids;
mp_params_exotrans('psi') = psi;

mp_params_typelife = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_typelife('epsilon') = epsilon;
mp_params_typelife('SS') = SS;

mp_params_intlen = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_params_intlen('n_jgrid') = n_jgrid;
mp_params_intlen('n_agrid') = n_agrid;
mp_params_intlen('n_etagrid') = n_etagrid;
mp_params_intlen('n_eta_S_grid') = n_eta_S_grid;
mp_params_intlen('n_eta_H_grid') = n_eta_H_grid;
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
mp_params = [...
    mp_params_covid_unemploy; ...
    mp_params_preftechpricegov; ...
    mp_params_statesgrid ; ...
    mp_params_exotrans; ...
    mp_params_typelife; ...
    mp_params_intlen ; ...
    mp_params_stat];

mp_params('mp_params_name') = string(st_param_group);

% MP_PARAMS = [MP_PARAMS_PREFTECHPRICE; MP_PARAMS_STATESGRID ; ...
%     MP_PARAMS_EXOTRANS; MP_PARAMS_TYPELIFE; MP_PARAMS_INTLEN];

%% Print
if (bl_print_mp_params)
    ff_container_map_display(mp_params_preftechpricegov);
    ff_container_map_display(mp_params_intlen);
    ff_container_map_display(mp_params_covid_unemploy, it_row_n_keep, it_col_n_keep);
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
    varargout{2} = mp_params_covid_unemploy;
    varargout{3} = mp_params_preftechpricegov;
    varargout{4} = mp_params_statesgrid;
    varargout{5} = mp_params_exotrans;
    varargout{6} = mp_params_typelife;
    varargout{7} = mp_params_intlen;
end

end
