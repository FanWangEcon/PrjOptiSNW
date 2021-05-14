%% The Stimulus Check Planning Problem
% The planner chooses the amount of stimulus checks for each group, where groups 
% are defined by marital status, number of children, income, and age in 2019.
%% 2019 Information Planning Problem
% Given the expected outcomes we computed conditional on 2019 information, we 
% can solve the planning problem. We have a number of different planning problems 
% that we solve given different individual level constraints and what the planner 
% can condition allocations on.
% 
% For FEASIBLE allocation, there are *970=5*2*97* types/cells of households:
%% 
% * 5 children groups
% * 2 spousal groups
% * 97 income bins: the allocation planner sees approximately $2500 income bins 
% between $0 and $238,800, and 1 bin after $238,800. There are 97 bins
%% 
% for OPTIMAL G4 (4 age groups 18 to 64) allocation, there are *3880=5*2*97*4* 
% types/cells of households:
%% 
% * 5 children groups
% * 2 spousal groups
% * 4 age groups
% * 97 income bins
%% 
% for OPTIMAL G47 (47 age groups) allocation, there are *45590=5*2*97*47* types/cells 
% of households:
%% 
% * 5 children groups
% * 2 spousal groups
% * 47 age groups
% * 97 income bins
%% 
% Optimal G4 has a + 1 version where we allocate for a fifth age group of individuals 
% older than 64 years of age. Optimal G47 has a + 35 version where optimal allocation 
% for all age groups are determined. 
%% Allocation Functions
% Functions in the <https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/alloc_discrete_fun_R 
% AllocateR/alloc_discrete_fun_R> folder of the project repository page is responsible 
% for feeding the dynamic programming results into the allocation functions. The 
% functions in this folder call the <https://fanwangecon.github.io/PrjOptiAlloc/reference/ffp_snw_process_inputs.html 
% ffp_snw_process_inputs> function to solve the allocation problems and compute 
% REV, and call the <https://fanwangecon.github.io/PrjOptiAlloc/reference/ffp_snw_graph_feasible.html 
% ffp_snw_graph_feasible> function to generate allocation graphs. These two functions 
% are a part of the <https://fanwangecon.github.io/PrjOptiAlloc/ PrjOptiAlloc> 
% package.