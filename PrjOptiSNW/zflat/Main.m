%% This is the main file that solves the household's optimization problem,
% aggregates, calibrates, and computes the planner's problem

clear
close all
clc

format longg

% Last updated: 06/18/2020

global beta gamma g_cons a2 rho_eta sigma_eta theta g_n cons_allocation_rule r agrid epsilon eta_grid SS pi_eta pi_kids pi_unemp psi Pop n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

%% Parameters
% Non-calibrated parameters
gamma=1; % Risk aversion parameter
rho_eta=0.98; % Ozkan (2017)
sigma_eta=0.02; % Ozkan (2017)
g_n=(1.01^2)-1; % Annual population growth of 1.1 percent
r=(1.04^2)-1; % Annual real interest rate of 4.0 percent from McGrattan and Prescott

% Government budget constraint parameters
g_cons=0.17575574; % Government consumption expenditures to GDP (BEA: Average 2015-2019) 
a2=3.664; % Initial guess for average income tax burden (if we use GS)

% Calibrated parameters
beta=0.96077; % Discount factor
theta=0.42315; % TFP parameter to normalize units such that GDP per capita equals 1 in the model: Real GDP/capita in 2019: $58,056

% Consumption allocation rule (1=uniform; 2=square root)
cons_allocation_rule=2;

% Number of grid points
n_jgrid=42; % Age runs from 18 to 100 (a period is 2 years)
n_agrid=40; % No. of grid points for assets
n_etagrid=7; % No. of grid points for persistent labor productivity shocks
n_educgrid=2; % No. of grid points for educational attainment (college vs. non-college)
n_marriedgrid=2; % No. of grid points for marital status
n_kidsgrid=6; % No. of grid points for children (0 to 5+ children)

% Social Security benefits
SS=zeros(n_jgrid,2);
SS(25:end,1)=0.24433; % Average SS non-college 2005-2009 as a share of GDP per capita
SS(25:end,2)=0.29263; % Average SS college 2005-2009 as a share of GDP per capita

%% Data
load('Mortality_prob_by_age_20_99.mat','mort_prob') % Age-specific mortality probabilities (20-99 year-olds)
psi_aux=1-mort_prob;

% Convert to two-year survival probabilities
psi=NaN(40,1);
for i=1:40
    psi(i)=prod(psi_aux((2*i-1):2*i));
end

psi=[psi(1);psi]; % Let survival probability of 18-year-olds be the same as that for 20-year-olds
psi=[psi;0]; % Maximum lifespan=100 (survival probability at age 100=0)

clear mort_prob psi_aux

load('Life_cycle_prod_by_educ.mat','life_cycle_prod_by_educ') % Life-cycle labor productivity for 20-100 year-olds by education (non-college vs. college)
epsilon_aux=life_cycle_prod_by_educ;

epsilon=NaN(41,2);
for i=1:40
    for e=1:2
        epsilon(i,e)=sum(epsilon_aux((2*i-1):2*i,e))/2;
    end
end

epsilon(end,:)=epsilon_aux(end,:);

epsilon=[epsilon(1,:);epsilon]; % Let life-cycle labor productivity of 18-year-olds be the same as that for 20-year-olds
epsilon(25:end,:)=0; % Assume zero labor productivity for 65+ year-olds (exogenous retirement)

clear life_cycle_prod_by_educ epsilon_aux

% Transition probabilities for number of children (0, 1, 2, 3, 4, or 5) (stored in the following order:
% Number of children in year 1, age, educational status, marital status. Each column refers to the number of children
% in year 2)
load('pi_kids_trans_prob','pi_kids_trans_prob')
pi_kids=NaN(n_kidsgrid,n_kidsgrid,n_jgrid,n_marriedgrid);

for kidsp=1:n_kidsgrid % No. of kids in year 2
    
    counter=0;
    
    for kids=1:n_kidsgrid % No. of kids in year 1
        for j=1:(n_jgrid-1) % Age in year 1
            for married=1:n_marriedgrid % Marital status

                counter=counter+1;
                pi_kids(kids,kidsp,j,married)=pi_kids_trans_prob(counter,kidsp);

            end
        end
    end
    
end

pi_kids(:,:,n_jgrid,:)=0;
pi_kids(:,1,n_jgrid,:)=1;

% Ensure that all rows sum to 1 in case of rounding error
for kids=1:n_kidsgrid % No. of kids in year 1
    for j=1:(n_jgrid-1) % Age in year 1
        for married=1:n_marriedgrid % Marital status
            aux_sum=sum(pi_kids(kids,:,j,married));
            pi_kids(kids,:,j,married)=pi_kids(kids,:,j,married)/aux_sum;
        end
    end
end

clear aux_sum counter pi_kids_trans_prob

%% Specifying asset grid (use non-linear spacing with minimum value of 0)
curv=3; % Governs how the grid points are allocated
scale_a=50; % Maximum value of assets (NOTE: Verifying that it does not bind in Aggregation.m)

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
stat_distr_kids(1,1,1)=0.7333;
stat_distr_kids(1,1,2)=0.1513;
stat_distr_kids(1,1,3)=0.0828;
stat_distr_kids(1,1,4)=0.0236;
stat_distr_kids(1,1,5)=0.0059;
stat_distr_kids(1,1,6)=0.0030;

aux=sum(stat_distr_kids(1,1,:));
stat_distr_kids(1,1,:)=stat_distr_kids(1,1,:)/aux;

% Not married but college-educated
% tab kids if Rmarried==0 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(2,1,1)=0.9752;
stat_distr_kids(2,1,2)=0.0236;
stat_distr_kids(2,1,3)=0.0001;
stat_distr_kids(2,1,4)=0.0011;
stat_distr_kids(2,1,5)=0;
stat_distr_kids(2,1,6)=0;

aux=sum(stat_distr_kids(2,1,:));
stat_distr_kids(2,1,:)=stat_distr_kids(2,1,:)/aux;

% Married and no college
% tab kids if Rmarried==1 & Rcollege==0 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(1,2,1)=0.4143;
stat_distr_kids(1,2,2)=0.2958;
stat_distr_kids(1,2,3)=0.2131;
stat_distr_kids(1,2,4)=0.0569;
stat_distr_kids(1,2,5)=0.0199;
stat_distr_kids(1,2,6)=0;

aux=sum(stat_distr_kids(1,2,:));
stat_distr_kids(1,2,:)=stat_distr_kids(1,2,:)/aux;

% Married and college-educated
% tab kids if Rmarried==1 & Rcollege==1 & inrange(RAGE,18,25) [aweight=WEIGHT]
stat_distr_kids(2,2,1)=0.7534;
stat_distr_kids(2,2,2)=0.2153;
stat_distr_kids(2,2,3)=0.0221;
stat_distr_kids(2,2,4)=0.0091;
stat_distr_kids(2,2,5)=0;
stat_distr_kids(2,2,6)=0;

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

name='Old-age dependency ratio=';
name2=[name,num2str(sum(Pop(25:end))/sum(Pop(1:24)))];
disp(name2);

%% Calibration
err=1;
tol=0.005;

while err>tol
    
    disp('Start calibration')

    % Solve optimization problem
    [V,ap,cons,exitflag]=VFI(A_aux,B_aux,Aeq,Beq,nonlcon,options);
    % Uncomment for grid search method for ap rather than continuous choice
    %[V_VFI,ap_VFI,cons_VFI,exitflag_VFI]=VFI_grid_search;

    % Aggregation
    [Phi_true,Phi_adj,A_agg,Y_inc_agg]=Aggregation(ap,stat_distr_eta,stat_distr_educ,stat_distr_married,stat_distr_kids);
    
    name='Income per capita (target=1)=';
    name2=[name,num2str(Y_inc_agg/sum(Pop))];
    disp(name2);
    name='Aggregate wealth to aggregate income (target=1.5)=';
    name2=[name,num2str(A_agg/Y_inc_agg)];
    disp(name2);
    
    err1=abs((Y_inc_agg/sum(Pop))-1); % Target: Income per capita equals 1
    err2=abs((A_agg/Y_inc_agg)-1.5); % Target: Annual capital/income ratio of 3
        
    err=max(err1,err2);
    
    param_update=[theta;beta];
    
    if err>tol
   
        theta=theta*((sum(Pop)/Y_inc_agg)^0.1); % Normalize theta such that income per capita equals 1
        beta=beta*((1.5/(A_agg/Y_inc_agg))^0.1); % Calibrate beta such that annual capital/income ratio equals 3
        
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

disp('Done with calibration')

%% Save value and policy functions from stationary distribution
%load('Value_and_policy_functions_ss','V','ap','cons');
save('Value_and_policy_functions_ss','V','ap','cons');

Output=NaN(n_jgrid*n_agrid*n_etagrid*n_educgrid*n_marriedgrid*n_kidsgrid,10);

counter=0;

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % Number of kids
                   
                       counter=counter+1;

                       Output(counter,1)=16+(2*j);
                       Output(counter,2)=agrid(a);
                       Output(counter,3)=eta_grid(eta);
                       Output(counter,4)=educ-1;
                       Output(counter,5)=married-1;
                       Output(counter,6)=kids-1;
                       Output(counter,7)=V(j,a,eta,educ,married,kids);
                       Output(counter,8)=cons(j,a,eta,educ,married,kids);
                       Output(counter,9)=ap(j,a,eta,educ,married,kids);
                       
                       if Phi_true(j,a,eta,educ,married,kids)>0
                           Output(counter,10)=Phi_true(j,a,eta,educ,married,kids)/sum(sum(sum(sum(sum(sum(Phi_true))))));
                       else
                           Output(counter,10)=0;
                       end
                       
                   end
               end
           end
       end
   end
end
                   
% Save output for computation of optimal allocation
disp('Save value and policy functions')
writematrix(Output,'Value_and_policy_functions_steady_state.csv')

clear Output

%% Asset distribution
asset_distr=zeros(n_agrid,2);
asset_distr(:,1)=agrid;
for a=1:n_agrid
   asset_distr(a,2)=sum(sum(sum(sum(sum(Phi_true(:,a,:,:,:,:))))))/sum(Pop); 
end

%% Age profiles
assets_avg=zeros(n_jgrid,1);
cons_avg=zeros(n_jgrid,1);
inc_avg=zeros(n_jgrid,1);
nr_of_kids=zeros(n_jgrid,n_kidsgrid);

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids

                       assets_avg(j)=assets_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*agrid(a);
                       cons_avg(j)=cons_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*cons(j,a,eta,educ,married,kids);
                       
                       [inc,earn]=individual_income(j,a,eta,educ);
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                       
                       inc_avg(j)=inc_avg(j)+Phi_adj(j,a,eta,educ,married,kids)*( inc+(married-1)*spouse_inc );
                       
                       nr_of_kids(j,kids)=nr_of_kids(j,kids)+Phi_adj(j,a,eta,educ,married,kids);
                   end
               end
           end
       end
   end
end

% Age profiles by marital status
Phi_adj2=zeros(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
for j=1:n_jgrid % Age
    for married=1:n_marriedgrid % Marital status
        dummy=sum(sum(sum(sum(Phi_true(j,:,:,:,married,:)))));
        for a=1:n_agrid % Assets
           for eta=1:n_etagrid % Productivity
               for educ=1:n_educgrid % Educational level
                   for kids=1:n_kidsgrid % No. of kids
                       if dummy>0
                           Phi_adj2(j,a,eta,educ,married,kids)=Phi_true(j,a,eta,educ,married,kids)/dummy;
                       else
                           Phi_adj2(j,a,eta,educ,married,kids)=0;
                       end
                   end
               end
           end
        end
    end
end

assets_avg_marr=zeros(n_jgrid,2);
cons_avg_marr=zeros(n_jgrid,2);
inc_avg_marr=zeros(n_jgrid,2);
nr_of_kids_marr=zeros(n_jgrid,n_kidsgrid,2);

for j=1:n_jgrid % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids

                       assets_avg_marr(j,married)=assets_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*agrid(a);
                       cons_avg_marr(j,married)=cons_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*cons(j,a,eta,educ,married,kids);
                       
                       [inc,earn]=individual_income(j,a,eta,educ);
                       spouse_inc=spousal_income(j,educ,kids,earn,SS(j,educ));
                       
                       inc_avg_marr(j,married)=inc_avg_marr(j,married)+Phi_adj2(j,a,eta,educ,married,kids)*( inc+(married-1)*spouse_inc );
                       
                       nr_of_kids_marr(j,kids,married)=nr_of_kids_marr(j,kids,married)+Phi_adj2(j,a,eta,educ,married,kids);
                       
                   end
               end
           end
       end
   end
end

clear dummy Phi_adj2

%% Probability of unemployment
pi_j=[0.22;0.175;0.16;0.165;0.22]; % Probability of unemployment in 2020 by age groups from Cajner et al. (2020, NBER)
pi_w=[0.360;0.22;0.17;0.14;0.09]; % Probability of unemployment in 2020 by wage quintiles from Cajner et al. (2020, NBER)
[age_factor,cutoffs]=pi_unemp_calibration(Phi_true,pi_j,pi_w);

pi_unemp=zeros(n_jgrid,5);
for i=1:5
    pi_unemp(1:7,i)=age_factor(1)*pi_w(i);
    pi_unemp(8:12,i)=age_factor(2)*pi_w(i);
    pi_unemp(13:17,i)=age_factor(3)*pi_w(i);
    pi_unemp(18:22,i)=age_factor(4)*pi_w(i);
    pi_unemp(23:24,i)=age_factor(5)*pi_w(i);
end

%% Compute value of employment and unemployment in 2020 conditional on number of welfare checks
% "Manna-from-heaven" where taxes do not change
xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)

% Compute policy functions in the event of unemployment. Required to compute V_U in the Planner's problem
disp('Compute policy functions in the event of unemployment')
[V_unemp,~,~,~]=VFI_unemp(A_aux,B_aux,Aeq,Beq,nonlcon,options,V,xi,b); 

TR=100/58056; % Value of a welfare check (can receive multiple checks). TO DO: Update with alternative values

n_welfchecksgrid=51; % Number of welfare checks. 0 refers to 0 dollars; 51 refers to 5000 dollars

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);
V_U=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid);

disp('Solve for V_W and V_U for different number of welfare checks')
for welf_checks=0:(n_welfchecksgrid-1)
    [V_W(:,:,:,:,:,:,welf_checks+1),~]=V_working_proxy(welf_checks,TR,V,options2);
    [V_U(:,:,:,:,:,:,welf_checks+1),~]=V_unemp_proxy(welf_checks,TR,xi,b,V_unemp,options2);

    name='Welfare checks=';
    name2=[name,num2str(welf_checks)];
    disp(name2)
end

%% Compute value of employment and unemployment in 2020 conditional on number of welfare checks
% Taxes are fully adjusted in 2020 to balance the government budget
% xi=0.5; % Proportional reduction in income due to unemployment (xi=0 refers to 0 labor income; xi=1 refers to no drop in labor income)
% b=0; % Unemployment insurance replacement rate (b=0 refers to no UI benefits; b=1 refers to 100 percent labor income replacement)
% 
% % Find tax rate that balances government budget given total spending on
% % unemployment bebenfits and welfare checks
% omega=0.025; % Total spending on welfare checks as a share of aggregate income
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
% disp('Solve for V_W and V_U for different number of welfare checks')
% for welf_checks=0:(n_welfchecksgrid-1)
%     [V_W(:,:,:,:,:,:,welf_checks+1),~]=V_working_proxy_tax(welf_checks,TR,V_working_tax,a2_COVID,options2);
%     [V_U(:,:,:,:,:,:,welf_checks+1),~]=V_unemp_proxy_tax(welf_checks,TR,xi,b,V_unemp_tax,a2_COVID,options2);
% 
%     name='Welfare checks=';
%     name2=[name,num2str(welf_checks)];
%     disp(name2)
% end

%% Solve planner's problem
n_incgrid=101; % Number of income groups
inc_grid=linspace(0,5,n_incgrid)'; % 5 refers to 5*58056=290280 dollars in 2012USD

V_planner=NaN(n_jgrid,n_marriedgrid,n_kidsgrid,n_welfchecksgrid,n_incgrid,n_incgrid);
Phi_mass=NaN(n_jgrid,n_marriedgrid,n_kidsgrid,n_incgrid,n_incgrid);

disp('Solve the problem of the Planner')
for j=1:(n_jgrid-1) % Age
   for married=1:n_marriedgrid % Marital status
       for kids=1:n_kidsgrid % Number of kids
           for welf_checks=0:(n_welfchecksgrid-1)
               for inc_group=1:n_incgrid
                   
                   if inc_group<n_incgrid
                       [V_planner(j,married,kids,welf_checks+1,inc_group,inc_group+1),Phi_mass(j,married,kids,inc_group,inc_group+1)]=Planner(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),inc_grid(inc_group+1),V_U,V_W,ap,cutoffs);
                   elseif inc_group==n_incgrid
                       [V_planner(j,married,kids,welf_checks+1,inc_group,inc_group),Phi_mass(j,married,kids,inc_group,inc_group)]=Planner(Phi_true,j,married,kids,welf_checks,inc_grid(inc_group),10E30,V_U,V_W,ap,cutoffs);
                   end
                   
               end
           end
       end
   end
   
   disp(j)
   
end

% Output for computing optimal allocation
Output=NaN((n_jgrid-1)*n_marriedgrid*n_kidsgrid*n_welfchecksgrid*inc_group,9);

counter=0;

for j=1:(n_jgrid-1) % Age
   for married=1:n_marriedgrid % Marital status
       for kids=1:n_kidsgrid % Number of kids
           for welf_checks=0:(n_welfchecksgrid-1)
               for inc_group=1:n_incgrid
                   
                   counter=counter+1;
                   
                   Output(counter,1)=16+(2*j);
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
disp('Save output for computation of optimal allocation')
writematrix(Output,'Output.csv')



%% Plot policy functions for next period's assets (TO DO: Needs to be updated)
% figure(1);
% hold on;
% plot(agrid(1:end),agrid(1:end),'k--',agrid(1:end),ap(21,:,1,1,1,1),'r-',agrid(1:end),ap(21,:,4,1,1,1),'b-',agrid(1:end),ap(21,:,7,1,1,1),'g');
% l1 = legend('45deg','\eta=1','\eta=4','\eta=7','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('policy function for tomorrows capital');
% title('Policy function for ap(j=1,:,\eta,1,1,1)');
% hold off;
% 
% figure(2);
% hold on;
% plot(agrid(1:end),ap(1,:,1,1,1,1),'ro',agrid(1:end),ap(1,:,4,1,1,1),'b+',agrid(1:end),ap(1,:,7,1,1,1),'k-');
% l1 = legend('\eta=1','\eta=4','\eta=7','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('policy function for tomorrows capital');
% title('Policy function for ap(j=1,:,\eta,1,1,1)');
% hold off;
% 
% figure(3);
% hold on;
% plot(agrid(1:end),agrid(1:end),'k--',agrid(1:end),ap(1,:,4,1,1,1),'r-',agrid(1:end),ap(5,:,4,1,1,1),'b-',agrid(1:end),ap(11,:,4,1,1,1),'g');
% l1 = legend('45deg','j=1','j=5','j=11','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('policy function for tomorrows capital');
% title('Policy function for ap(j,:,4,1,1,1)');
% hold off;

% figure(1);
% hold on;
% plot(agrid(1:end),V(1,:,1,1,1,1),'r-',agrid(1:end),V(1,:,4,1,1,1),'b-',agrid(1:end),V(1,:,7,1,1,1),'g');
% l1 = legend('\eta=1','\eta=4','\eta=7','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('Value function for tomorrows capital');
% title('Value function for ap(j=1,:,\eta,1,1,1)');
% hold off;
% 
% figure(1);
% hold on;
% plot(agrid(1:end),V(1,:,4,1,1,1),'r-',agrid(1:end),V(1,:,4,1,2,1));
% l1 = legend('married=0','married=1','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('Value function for tomorrows capital');
% title('Value function for (j=1,:,4,1,m,1)');
% hold off;
% 
% figure(2);
% hold on;
% plot(agrid(1:end),cons(1,:,4,1,1,1),'r-',agrid(1:end),cons(1,:,4,1,2,1));
% l1 = legend('married=0','married=1','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('Policy function for consumption');
% title('Policy function for cons(j=1,:,4,1,m,1)');
% hold off;
%
% figure(4);
% hold on;
% plot(agrid(1:end),ap(21,:,4,1,1,1),'ro',agrid(1:end),ap(21,:,4,2,1,1),'b+',agrid(1:end),ap(21,:,4,1,2,1),'k-',agrid(1:end),ap(21,:,4,2,2,1),'g*');
% l1 = legend('e=1&m=1','e=2&m=1','e=1&m=2','e=2&m=2','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('policy function for tomorrows capital');
% title('Policy function for ap(21,:,4,e,m,1)');
% hold off;
% 
% figure(5);
% hold on;
% plot(agrid(1:end),ap(21,:,4,1,1,1),'ro',agrid(1:end),ap(21,:,4,1,1,2),'b+',agrid(1:end),ap(21,:,4,1,1,3),'k-',agrid(1:end),ap(21,:,4,1,1,4),'g*',agrid(1:end),ap(21,:,4,1,1,5),'y',agrid(1:end),ap(21,:,4,1,1,6),'b');
% l1 = legend('k=1','k=2','k=3','k=4','k=5','k=6','Location','NorthWest');
% set(l1,'box','off');
% xlabel('capital grid');
% ylabel('policy function for tomorrows capital');
% title('Policy function for ap(21,:,4,1,1,m)');
% hold off;