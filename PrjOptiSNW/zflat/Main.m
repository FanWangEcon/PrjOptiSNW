%% This is the main file that solves the household's optimization problem,
% aggregates, calibrates, and computes the planner's problem

clear
close all
clc

format longg

% Last updated: 07/20/2020

global beta gamma g_cons a2 rho_eta sigma_eta theta g_n cons_allocation_rule r agrid epsilon eta_H_grid eta_S_grid SS pi_eta pi_kids pi_unemp psi Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid jret Bequests bequests_option throw_in_ocean

%% Parameters
run_calibration=0; % Whether or not to calibrate theta and beta. run_calibration=1 recalibrates the parameters

% Non-calibrated parameters
gamma=2; % Risk aversion parameter
rho_eta=0.98; % Persistence of AR(1) productivity shocks
sigma_eta=0.018; % Variance of AR(1) productivity shocks
g_n=0.01; % Annual population growth of 1.1 percent
r=0.04; % Annual real interest rate of 4.0 percent from McGrattan and Prescott

rho_eta_spouse=0; % Persistence of spousal AR(1) productivity shocks
sigma_eta_spouse=1.040654^2; % Variance of spousal AR(1) productivity shocks (standard deviation of residual from spousal income regression. See spousal_income.m for regression specification details)

% Calibrated parameters
beta=0.969010845568155; % Discount factor
theta=0.66947637485374; % TFP parameter to normalize units such that average household income relative to GDP per capita is consistent with the data (the latter is normalized to 1): Real GDP/capita in 2019: $58,056

% Consumption allocation rule (1=uniform; 2=square root; 0=equivalent to only 1 household member independent of family size)
cons_allocation_rule=2;

% Bequests allocation rule (=1: accidental bequests go to the government; =2: accidental bequests uniformly across the population)
bequests_option=1;
Bequests=0.05826*(bequests_option-1);
throw_in_ocean=1; % If bequests go to the government, a value of 1 for throw_in_ocean means that all accidental bequests are "thrown in the ocean", whereas a value of 0 means the full amount goes to the government

% Government budget constraint parameters
g_cons=0.17575574; % Government consumption expenditures to GDP (BEA: Average 2015-2019)

% How accidental bequests are handled
if bequests_option==2
    a2=1.575;
elseif bequests_option==1
    if throw_in_ocean==0
        a2=0.7027;
    elseif throw_in_ocean==1
        a2=1.6081061107858;
    end
end

% Number of grid points
n_jgrid=83; % Age runs from 18 to 100 (a period is 1 year)
n_agrid=151; % 351; % No. of grid points for assets
n_eta_H_grid=9; % No. of grid points for persistent labor productivity shocks
n_eta_S_grid=5; % 1; % No. of grid points for spousal labor productivity shocks (=1 corresponds to deterministic spousal income)
n_etagrid=n_eta_H_grid*n_eta_S_grid; % Length of productivity shock grid
n_educgrid=2; % No. of grid points for educational attainment (college vs. non-college)
n_marriedgrid=2; % No. of grid points for marital status
n_kidsgrid=5; % No. of grid points for children (0 to 4+ children)

% Retirement age (age 65)
jret=48;

% Social Security benefits
SS=zeros(n_jgrid,2);
SS(jret:end,1)=0.24433; % Average SS non-college 2005-2009 as a share of GDP per capita
SS(jret:end,2)=0.29263; % Average SS college 2005-2009 as a share of GDP per capita

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

%% Data
load('Mortality_prob_by_age_18_99.mat','mort_prob') % Age-specific mortality probabilities (18-99 year-olds)
psi=1-mort_prob;
psi=[psi;0]; % Maximum lifespan=100 (survival probability at age 100=0)

%clear mort_prob psi_aux
clear mort_prob

load('Life_cycle_prod_by_educ.mat','life_cycle_prod_by_educ') % Life-cycle labor productivity for 20-100 year-olds by education (non-college vs. college)
epsilon=NaN(83,2);
epsilon(3:end,:)=life_cycle_prod_by_educ(:,:);

epsilon(1,:)=epsilon(3,:); % Let life-cycle labor productivity of 18- and 19-year-olds be the same as that for 20-year-olds
epsilon(2,:)=epsilon(3,:);

epsilon(jret:end,:)=0; % Assume zero labor productivity for 65+ year-olds (exogenous retirement)

%clear life_cycle_prod_by_educ epsilon_aux
clear life_cycle_prod_by_educ

% Transition probabilities for number of children (0, 1, 2, 3, 4, or 5) (stored in the following order:
% Number of children in year 1, age, marital status, college attainment. Each column refers to the
% number of children in year 2)
%load('pi_kids_trans_prob','pi_kids_trans_prob_2year_intervals')
load('pi_kids_trans_prob','pi_kids_trans_prob')
pi_kids=NaN(n_kidsgrid,n_kidsgrid,n_jgrid,n_educgrid,n_marriedgrid);

for kidsp=1:n_kidsgrid % No. of kids in year 2

    counter=0;

    for kids=1:n_kidsgrid % No. of kids in year 1
        for j=1:(n_jgrid-1) % Age in year 1
            for married=1:n_marriedgrid % Marital status
                for educ=1:n_educgrid % Educational level

                    counter=counter+1;
                    pi_kids(kids,kidsp,j,educ,married)=pi_kids_trans_prob(counter,kidsp);

                end
            end
        end
    end

end

pi_kids(:,:,n_jgrid,:,:)=0;
pi_kids(:,1,n_jgrid,:,:)=1;

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

clear aux_sum counter pi_kids_trans_prob

%% Specifying asset grid (use non-linear spacing with minimum value of 0)
curv=3; % Governs how the grid points are allocated
scale_a=135; % Maximum value of assets (NOTE: Verifying that it does not bind in Aggregation.m) grid=0

agrid=zeros(n_agrid,1);

for i=2:n_agrid
	agrid(i)=scale_a*((i-1)/(n_agrid-1))^curv;
end

%% For minimization problem
%amin=0;
%amax=agrid(end);
A_aux=[];
B_aux=[];
Aeq=[];
Beq=[];

nonlcon=[];
options=optimoptions('fmincon','Display', 'off');
options2=optimoptions('fsolve','Display','off');

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
stat_distr_kids(1,1,1)=0.7333;
stat_distr_kids(1,1,2)=0.1513;
stat_distr_kids(1,1,3)=0.0828;
stat_distr_kids(1,1,4)=0.0236;
stat_distr_kids(1,1,5)=0.0090;

aux=sum(stat_distr_kids(1,1,:));
stat_distr_kids(1,1,:)=stat_distr_kids(1,1,:)/aux;

% Not married but college-educated
% tab kids if Rmarried==0 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(2,1,1)=0.9752;
stat_distr_kids(2,1,2)=0.0236;
stat_distr_kids(2,1,3)=0.0001;
stat_distr_kids(2,1,4)=0.0011;
stat_distr_kids(2,1,5)=0;

aux=sum(stat_distr_kids(2,1,:));
stat_distr_kids(2,1,:)=stat_distr_kids(2,1,:)/aux;

% Married and no college
% tab kids if Rmarried==1 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(1,2,1)=0.4143;
stat_distr_kids(1,2,2)=0.2958;
stat_distr_kids(1,2,3)=0.2131;
stat_distr_kids(1,2,4)=0.0569;
stat_distr_kids(1,2,5)=0.0199;

aux=sum(stat_distr_kids(1,2,:));
stat_distr_kids(1,2,:)=stat_distr_kids(1,2,:)/aux;

% Married and college-educated
% tab kids if Rmarried==1 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(2,2,1)=0.7534;
stat_distr_kids(2,2,2)=0.2153;
stat_distr_kids(2,2,3)=0.0221;
stat_distr_kids(2,2,4)=0.0091;
stat_distr_kids(2,2,5)=0;

aux=sum(stat_distr_kids(2,2,:));
stat_distr_kids(2,2,:)=stat_distr_kids(2,2,:)/aux;

clear aux

%% Population distribution
% Normalize mass of 18-year-olds to 1
Pop=zeros(n_jgrid,1);
Pop(1)=1;
for j=2:n_jgrid
    Pop(j)=Pop(j-1)*psi(j-1)/(1+g_n);
end

name='Old-age dependency ratio (ratio of 65+/(18-64))=';
name2=[name,num2str(sum(Pop(48:end))/sum(Pop(1:47)))];
disp(name2);

%% Calibration
err=1;
tol=0.005;

if run_calibration==1

    disp('Start calibration')

    while err>tol

        it=1;

        while it>0

            % Solve optimization problem
            % Uncomment for continuous choice for ap
            tic;
            [V,ap,cons,exitflag]=VFI(A_aux,B_aux,Aeq,Beq,nonlcon,options);
            toc;
            % Uncomment for grid search method for ap rather than continuous choice
%             tic;
%             [V,ap,cons,exitflag]=VFI_grid_search;
%             toc;

            % Aggregation
             [Phi_true,Phi_adj,A_agg,Y_inc_agg,it]=Aggregation(ap,cons,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids);
%              [Phi_true,Phi_adj,A_agg,Y_inc_agg,it]=Aggregation_grid_search(ap,cons,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids);

        end

        name='Average household income (target=1.38)=';
        name2=[name,num2str(Y_inc_agg/sum(Pop))];
        disp(name2);
        name='Aggregate wealth to aggregate income (target=3.0)=';
        name2=[name,num2str(A_agg/Y_inc_agg)];
        disp(name2);

        err1=abs((Y_inc_agg/sum(Pop))-1.38); % Target: Average household income relative to income per capita (latter is normalized to 1 in the model)
        err2=abs((A_agg/Y_inc_agg)-3.0); % Target: Annual capital/income ratio of 3

        err=max(err1,err2);

        param_update=[theta;beta];

        if err>tol

            theta=theta*((1.38/(Y_inc_agg/sum(Pop)))^0.2); % Normalize theta such that income per capita equals 1
            beta=beta*((3.0/(A_agg/Y_inc_agg))^0.025); % Calibrate beta such that annual capital/income ratio equals 3

        end

        param_update=[param_update(1,1),theta;param_update(2,1),beta];

        name='Old and updated value of theta=';
        name2=[name,num2str(param_update(1,:))];
        disp(name2);

        name='Old and updated value of beta=';
        name2=[name,num2str(param_update(2,:))];
        disp(name2);

        disp([err1,err2])

    end

elseif run_calibration==0

    disp('Not calibrating the parameters')

    it=1;

    while it>0

        % Solve optimization problem
        % Uncomment for continuous choice for ap
        tic;
        [V,ap,cons,exitflag]=VFI(A_aux,B_aux,Aeq,Beq,nonlcon,options);
        toc;
        % Uncomment for grid search method for ap rather than continuous choice
%         tic;
%         [V,ap,cons,exitflag]=VFI_grid_search;
%         toc;

        % Aggregation
         [Phi_true,Phi_adj,A_agg,Y_inc_agg,it]=Aggregation(ap,cons,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids);
%          [Phi_true,Phi_adj,A_agg,Y_inc_agg,it]=Aggregation_grid_search(ap,cons,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids);

    end

end

disp('Done with calibration')

name='Calibrated parameters: beta,theta,a2=';
name2=[name,num2str([beta,theta,a2])];
disp(name2);

clear name name2

%% Save value and policy functions from stationary distribution
%load('Value_and_policy_functions_ss','V','ap','cons');
%save('Value_and_policy_functions_ss','V','ap','cons');

% Uncomment to store value and policy functions in csv-format
% Output=NaN(n_jgrid*n_agrid*n_etagrid*n_educgrid*n_marriedgrid*n_kidsgrid,10);
%
% counter=0;
%
% for j=1:n_jgrid % Age
%    for a=1:n_agrid % Assets
%        for eta=1:n_etagrid % Productivity
%            for educ=1:n_educgrid % Educational level
%                for married=1:n_marriedgrid % Marital status
%                    for kids=1:n_kidsgrid % Number of kids
%
%                        counter=counter+1;
%
%                        Output(counter,1)=17+j;
%                        Output(counter,2)=agrid(a);
%                        Output(counter,3)=eta_H_grid(eta);
%                        Output(counter,4)=educ-1;
%                        Output(counter,5)=married-1;
%                        Output(counter,6)=kids-1;
%                        Output(counter,7)=V(j,a,eta,educ,married,kids);
%                        Output(counter,8)=cons(j,a,eta,educ,married,kids);
%                        Output(counter,9)=agrid(ap(j,a,eta,educ,married,kids));
% %                        Output(counter,9)=ap(j,a,eta,educ,married,kids);
%
%                        if Phi_true(j,a,eta,educ,married,kids)>0
%                            Output(counter,10)=Phi_true(j,a,eta,educ,married,kids)/sum(sum(sum(sum(sum(sum(Phi_true))))));
%                        else
%                            Output(counter,10)=0;
%                        end
%
%                    end
%                end
%            end
%        end
%    end
% end
%
% % Save output for computation of optimal allocation
% disp('Save value and policy functions')
% writematrix(Output,'Value_and_policy_functions_steady_state.csv')
%
% clear Output

%% Asset distribution
% asset_distr=zeros(n_agrid,2);
% asset_distr(:,1)=agrid;
% for a=1:n_agrid
%    asset_distr(a,2)=sum(sum(sum(sum(sum(Phi_true(:,a,:,:,:,:))))))/sum(Pop);
% end

%% Age profiles
% Uncomment to compute age profiles
% assets_avg=zeros(n_jgrid,1);
% cons_avg=zeros(n_jgrid,1);
% inc_avg=zeros(n_jgrid,1);
% nr_of_kids=zeros(n_jgrid,n_kidsgrid);
%
% for j=1:n_jgrid % Age
%    for a=1:n_agrid % Assets
%        for eta=1:n_etagrid % Productivity
%            for educ=1:n_educgrid % Educational level
%                for married=1:n_marriedgrid % Marital status
%                    for kids=1:n_kidsgrid % No. of kids
%
%                        assets_avg(j)=assets_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*agrid(a);
%                        cons_avg(j)=cons_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*cons(j,a,eta,educ,married,kids);
%
%                        [inc,earn]=individual_income(j,a,eta,educ);
%                        spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
%
%                        inc_avg(j)=inc_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*( inc+(married-1)*spouse_inc*exp(eta_S_grid(eta)) );
%
%                        nr_of_kids(j,kids)=nr_of_kids(j,kids)+Phi_adj(j,a,eta,educ,married,kids);
%                    end
%                end
%            end
%        end
%    end
% end
%
% % Age profiles by marital status
% Phi_adj2=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
% for j=1:n_jgrid % Age
%     for married=1:n_marriedgrid % Marital status
%         dummy=sum(sum(sum(sum(Phi_true(j,:,:,:,married,:)))));
%         for a=1:n_agrid % Assets
%            for eta=1:n_etagrid % Productivity
%                for educ=1:n_educgrid % Educational level
%                    for kids=1:n_kidsgrid % No. of kids
%                        if dummy>0
%                            Phi_adj2(j,a,eta,educ,married,kids)=Phi_true(j,a,eta,educ,married,kids)/dummy;
%                        else
%                            Phi_adj2(j,a,eta,educ,married,kids)=0;
%                        end
%                    end
%                end
%            end
%         end
%     end
% end
%
% assets_avg_marr=zeros(n_jgrid,2);
% cons_avg_marr=zeros(n_jgrid,2);
% inc_avg_marr=zeros(n_jgrid,2);
% nr_of_kids_marr=zeros(n_jgrid,n_kidsgrid,2);
%
% for j=1:n_jgrid % Age
%    for a=1:n_agrid % Assets
%        for eta=1:n_etagrid % Productivity
%            for educ=1:n_educgrid % Educational level
%                for married=1:n_marriedgrid % Marital status
%                    for kids=1:n_kidsgrid % No. of kids
%
%                        assets_avg_marr(j,married)=assets_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*agrid(a);
%                        cons_avg_marr(j,married)=cons_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*cons(j,a,eta,educ,married,kids);
%
%                        [inc,earn]=individual_income(j,a,eta,educ);
%                        spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
%
%                        inc_avg_marr(j,married)=inc_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*( inc+(married-1)*spouse_inc*exp(eta_S_grid(eta)) );
%
%                        nr_of_kids_marr(j,kids,married)=nr_of_kids_marr(j,kids,married)+Phi_adj2(j,a,eta,educ,married,kids);
%
%                    end
%                end
%            end
%        end
%    end
% end
%
% clear dummy Phi_adj2

%% Compute value of employment and unemployment in 2020 conditional on number of welfare checks: "Manna-from-heaven" where taxes do not change
xi=0.75; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=0; % 1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)

% Compute policy functions in the event of unemployment. Required to compute V_U in the Planner's problem
disp('Compute value function and policy functions in the event of unemployment')
tic;
[V_unemp,~,cons_unemp,~]=VFI_unemp(A_aux,B_aux,Aeq,Beq,nonlcon,options,V,xi,b);
% [V_unemp,~,cons_unemp,~]=VFI_unemp_grid_search(V,xi,b);
toc;

TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
name='Value of each welfare check (in 2012USD)= ';
name2=[name,num2str(TR*58056)];
disp(name2)

clear name name2

n_welfchecksgrid=51; % Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
V_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);

C_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
C_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);

disp('Solve for V_W and V_U for different number of welfare checks')
for welf_checks=0:(n_welfchecksgrid-1)
    tic;
    [V_W(:,:,:,:,:,:,welf_checks+1),C_W(:,:,:,:,:,:,welf_checks+1),~]=V_working_proxy(welf_checks,TR,V,cons,options2);
    [V_U(:,:,:,:,:,:,welf_checks+1),C_U(:,:,:,:,:,:,welf_checks+1),~]=V_unemp_proxy(welf_checks,TR,xi,b,V_unemp,cons_unemp,options2);

    name='Welfare checks= ';
    name2=[name,num2str(welf_checks)];
    disp(name2)
    toc;
end

%% Probability of unemployment
disp('Update this section to allow pi_unemp() to be continuous in age and wage')
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

cutoffs=wage_cutoffs(Phi_true);

%% Compute value of employment and unemployment in 2020 conditional on number of welfare checks: Taxes are fully adjusted in 2020 to balance the government budget
% disp('This section must be updated!')
% xi=0.75; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
% b=1; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement. b=1 in benchmark model to capture that the generosity of UI benefits was increased under the CARES Act with the aim of replacing 100 percent of average lost earnings)
%
% % Find tax rate that balances government budget given total spending on
% % unemployment bebenfits and welfare checks
% omega=0.0135; % Total spending on welfare checks as a share of aggregate income (estimated to cost $290 billion, or 1.35 percent of GDP in 2019. IRS has paid out $218 billion as of May 11)
% a2_guess=a2; % Initial guess for a2
%
% a2_COVID=find_tax_rate(a2_guess,Phi_true,omega,xi,b,cutoffs);
%
% % Compute policy functions in the event of working and unemployment when the government adjusts taxes to balances the budget. Required to compute V_U and V_W in the Planner's problem
% disp('Compute policy functions when taxes adjust in the event of working')
% [V_working_tax,~,~,~]=VFI_working_tax(A_aux,B_aux,Aeq,Beq,nonlcon,options,V,a2_COVID);
% disp('Compute policy functions when taxes adjust in the event of unemployment')
% [V_unemp_tax,~,~,~]=VFI_unemp_tax(A_aux,B_aux,Aeq,Beq,nonlcon,options,V,xi,b,a2_COVID);
%
% TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values
%
% n_welfchecksgrid=51; % Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars
%
% V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
% V_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
%
% C_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
% C_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
%
% disp('Solve for V_W and V_U for different number of welfare checks')
% for welf_checks=0:(n_welfchecksgrid-1)
%     [V_W(:,:,:,:,:,:,welf_checks+1),C_W(:,:,:,:,:,:,welf_checks+1),~]=V_working_proxy_tax(welf_checks,TR,V_working_tax,cons_working_tax,a2_COVID,options2);
%     [V_U(:,:,:,:,:,:,welf_checks+1),C_U(:,:,:,:,:,:,welf_checks+1),~]=V_unemp_proxy_tax(welf_checks,TR,xi,b,V_unemp_tax,cons_unemp_tax,a2_COVID,options2);
%
%     name='Welfare checks=';
%     name2=[name,num2str(welf_checks)];
%     disp(name2)
% end

%% Solve planner's problem
n_incgrid=201; % Number of income groups
n_incgrid_aux=round(0.75*n_incgrid);
inc_grid1=linspace(0,4,n_incgrid_aux)'; % 4 refers to 4*58056=232224 dollars in 2012USD
inc_grid=[inc_grid1;linspace(4+((7-4)/(n_incgrid-n_incgrid_aux)),7,n_incgrid-n_incgrid_aux)']; % 7 refers to 7*58056=406392 dollars in 2012USD

clear inc_grid1 n_incgrid_aux

% Income and earning grids used to speed up the code
ref_earn_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
inc_tot_grid=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
Phi_true_1=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid
   for a=1:n_agrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid
                       [inc,earn]=individual_income(j,a,eta,educ);
                       ref_earn_grid(j,a,eta,educ,married,kids)=earn;
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                       inc_tot_grid(j,a,eta,educ,married,kids)=inc+(married-1)*spouse_inc*exp(eta_S_grid(eta));
                       if Phi_true(j,a,eta,educ,married,kids)>0
                           Phi_true_1(j,a,eta,educ,married,kids)=Phi_true(j,a,eta,educ,married,kids)/sum(Phi_true,'all');
                       end
                   end
               end
           end
       end
   end
end

EV=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,5);
for j=1:n_jgrid
   for a=1:n_agrid
       for eta=1:n_etagrid
           for educ=1:n_educgrid
               for married=1:n_marriedgrid
                   for kids=1:n_kidsgrid
                       for welf_checks=1:n_welfchecksgrid
                            for wage_ind=1:5
                                EV(j,a,eta,educ,married,kids,welf_checks,wage_ind)=pi_unemp(j,wage_ind)*V_U(j,a,eta,educ,married,kids,welf_checks)+(1-pi_unemp(j,wage_ind))*V_W(j,a,eta,educ,married,kids,welf_checks);
                            end
                       end
                   end
               end
           end
       end
   end
end

V_planner=NaN(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid);
Phi_mass=NaN(n_jgrid-1,n_marriedgrid,n_kidsgrid,n_incgrid);

disp('Solve the problem of the Planner')
tic;
for j=1:(n_jgrid-1) % Age
   for married=1:n_marriedgrid % Marital status
       for kids=1:n_kidsgrid % Number of kids
           for welf_checks=0:(n_welfchecksgrid-1)
               for inc_group=1:n_incgrid

                   if inc_group<n_incgrid
                      [V_planner(j,married,kids,welf_checks+1,inc_group),Phi_mass(j,married,kids,inc_group)]=Planner(Phi_true_1,j,married,kids,welf_checks,inc_grid(inc_group),inc_grid(inc_group+1),ap,cutoffs,ref_earn_grid,inc_tot_grid,EV,pi_eta,pi_kids,agrid,n_agrid,n_etagrid,n_educgrid,n_kidsgrid);
%                        [V_planner(j,married,kids,welf_checks+1,inc_group),Phi_mass(j,married,kids,inc_group)]=Planner_grid_search(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),inc_grid(inc_group+1),V_U,V_W,ap,cutoffs,ref_inc_grid,ref_earn_grid,spouse_inc_grid);

                   elseif inc_group==n_incgrid
                      [V_planner(j,married,kids,welf_checks+1,inc_group),Phi_mass(j,married,kids,inc_group)]=Planner(Phi_true_1,j,married,kids,welf_checks,inc_grid(inc_group),10E30,ap,cutoffs,ref_earn_grid,inc_tot_grid,EV,pi_eta,pi_kids,agrid,n_agrid,n_etagrid,n_educgrid,n_kidsgrid);
%                        [V_planner(j,married,kids,welf_checks+1,inc_group),Phi_mass(j,married,kids,inc_group)]=Planner_grid_search(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),10E30,V_U,V_W,ap,cutoffs,ref_inc_grid,ref_earn_grid,spouse_inc_grid);

                   end

               end
           end
       end
   end

   disp(j)

end
toc;

% Output for computing optimal allocation
Output=NaN((n_jgrid-1)*n_marriedgrid*n_kidsgrid*n_welfchecksgrid*inc_group,9);

counter=0;

for j=1:(n_jgrid-1) % Age
   for married=1:n_marriedgrid % Marital status
       for kids=1:n_kidsgrid % Number of kids
           for welf_checks=0:(n_welfchecksgrid-1) % Number of welfare checks
               for inc_group=1:n_incgrid

                   counter=counter+1;

                   Output(counter,1)=17+j;
                   Output(counter,2)=married-1;
                   Output(counter,3)=kids-1;
                   Output(counter,4)=welf_checks;

                   if inc_group<n_incgrid
                       Output(counter,5)=inc_grid(inc_group)*58056;
                       Output(counter,6)=inc_grid(inc_group+1)*58056;
                       Output(counter,7)=Phi_mass(j,married,kids,inc_group,inc_group+1);
                       Output(counter,9)=V_planner(j,married,kids,welf_checks+1,inc_group,inc_group+1);
                   elseif inc_group==n_incgrid
                       Output(counter,5)=inc_grid(inc_group)*58056;
                       Output(counter,6)=10E30;
                       Output(counter,7)=Phi_mass(j,married,kids,inc_group,inc_group);
                       Output(counter,9)=V_planner(j,married,kids,welf_checks+1,inc_group,inc_group);
                   end

                   Output(counter,8)=psi(j);

               end
           end
       end
   end

end

% Save output for computation of optimal allocation
%disp('Save output for computation of optimal allocation')
%writematrix(Output,'Output.csv')
