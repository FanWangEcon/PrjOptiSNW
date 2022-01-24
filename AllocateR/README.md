# README

Work flow information and some optimal allocation information.

## Folders

- [alloc_discrete_fun_R](https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/alloc_discrete_fun_R): Optimal allocation files that call functions from the R [PrjOptiAlloc](https://fanwangecon.github.io/PrjOptiAlloc/) package and that take as inputs results from the dynamic programming simulation functions from Matlab.
- [alloc_discrete_paper_graphs](https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/alloc_discrete_paper_graphs): Given optimal allocation results, visualize for the paper manuscript.
- [jointbetedu](https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/jointbetedu): Files that merge CSV files solved for different beta and education types together, to be used as inputs for the optimal allocation files. 
- [sandbox](https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/sandbox): Various RMD files testing functions later included in [alloc_discrete_fun_R](https://github.com/FanWangEcon/PrjOptiSNW/tree/master/AllocateR/alloc_discrete_fun_R). Also tested code that became a part of functions [ffp_snw_process_inputs](https://fanwangecon.github.io/PrjOptiAlloc/reference/ffp_snw_process_inputs.html), [ffp_snw_graph_feasible](https://fanwangecon.github.io/PrjOptiAlloc/reference/ffp_snw_graph_feasible.html) and other functions in the [PrjOptiAlloc](https://fanwangecon.github.io/PrjOptiAlloc/) package.

## Work Flow

Work Flow is about how to use the various functions in the *AllocateR* folder to generate optimal allocation results.

### Support Function

First, optimal allocation relies on the allocation package [PrjOptiAlloc](https://fanwangecon.github.io/PrjOptiAlloc/), which needs to be separately installed and ran.

Second, run support function where path functionality is centralized.

- Code Path: PrjOptiSNW/AllocateR/fs_opti_support_202111.R

### Generate Main Allocation Inputs

**Step 1**

Solve the dynamic programming problem for several four beta and education groups separately.

- Code Path: PrjOptiSNW/invoke/202111/snw_res_b1_manna_mixture_arp.m
- Results Path: Dropbox/PrjNygaardSorensenWang/Output202111/

**Step 2**

Merge skill (education) and beta specific groups together:

- Code Path: AllocateR/jointbetedu/fs_wgt_merge_betaedutypes_202111.R
- Results Path: Dropbox/PrjNygaardSorensenWang/Output202111/

**Step 3**

Generate MPC and APC Files and MPV and APV Files

- Code Path: AllocateR/alloc_discrete_fun_R/fs_mpc_tables_increments_202111.R
- Results Path: dropbox/PrjNygaardSorensenWang/Results202111/apc_mpc/

Graph Jointly A and alpha.

- Code: AllocateR/alloc_discrete_paper_graphs/fs_A_alpha_scatter_202111.R
- Results Path: dropbox/PrjNygaardSorensenWang/Results202111/apc_mpc/

### Optimal Allocation Results Across Planners

#### Optimal Allocation Across Upper Bounds and Planner Types

**Step 4a**

Allocation Results, solve allocation problem across planner types and over different allocation bounds

- Code Path: AllocateR/alloc_discrete_fun_R/fs_opti_202111_bush.R
    - first run test with: *it_run_collection <- 1*
    - second run full problem with: *it_run_collection <- 2*
- Results Path: for example,
    - for rho = 1: dropbox/PrjNygaardSorensenWang/Results202111/capi_rh1/
    - for rho = -1: dropbox/PrjNygaardSorensenWang/Results202111/capi_rhn1/
    - for rho = -99: dropbox/PrjNygaardSorensenWang/Results202111/capi_rhn99/

A number of graphs are automatically generated, for each planner type and allocation bound.

**Step 4b**

Graphs to show joint allocation results across multiple planners.

- Code Path: AllocateR/alloc_discrete_paper_graphs/fs_opti_graph_multiplanner_202111.R
- Results Path: dropbox/PrjNygaardSorensenWang/Results202111/capi_rhmultiple/

#### Welfare Implication of Optimal Allocation Results

We have already generated welfare results in various folders in step *4a*, the task here is to combine results from various files together

**Step 5**

Visualize aggregate "REV" type results, how "welfare" in resource units change.

- Code Path: AllocateR/alloc_discrete_paper_graphs/fs_rev_202111.R
- Results Path: Dropbox/PrjNygaardSorensenWang/Results202111/capi_rhomultiple_rev/

### Robustness Exercises Optimal Allocaitons

Follow the steps above, first in *6a* combine beta/edu sub files together in different ways, then solve for optimal allocation, then combine to generate graphs in *6b*

**Step 6a**

Combine together various results for different robustness exercises

- Code Path: AllocateR/jointbetedu/fs_wgt_merge_betaedutypes_202111.R
- First: different weights of low beta and high beta
  > st_biden_or_trump <- 'bushchck'
  > ls_it_mixture_grp <- c(2,3,4,5)
- Second: Lower savings interest at r = 0.02
  > st_biden_or_trump <- 'bchklkr2'
  > ls_it_mixture_grp <- c(1)
- Third: no UI
  > st_biden_or_trump <- 'bcklknou'
  > ls_it_mixture_grp <- c(1)

Solve allocation problems from the First/Second/Third problems

- Code Path: AllocateR/alloc_discrete_fun_R/fs_opti_202111.R

> ls_st_file_jedc_2021_robust <- c(ls_st_file_jedc_2021_bchklock_mix, ls_st_file_jedc_2021_bchklock_ui_r)
> it_run_collection <- 4

**Step 6b**

Visualize robustness results.

- Code Path: /AllocateR/alloc_discrete_paper_graphs/fs_opti_graph_multisolufiles_2021.R

### Results in Perturb Folder

This is a special robustness exercises, where result for perturbing things are stored in separate folder.

**Step 7a**

Perturb allocation results, solve allocation problem for one planner type with perturbed A values, do it 500 times.

- Code Path (same as *4a*): AllocateR/alloc_discrete_fun_R/fs_opti_202111_bush.R
    - first run a simple test with: *it_run_collection <- 11*
    - second run 100 set of draws with graphs: *it_run_collection <- 12*
    - third run 400 sets of draws: *it_run_collection <- 13*
- Results Path: Dropbox/Results202111_perturb/

**Step 7b**

Graphs to show joint allocation results across multiple perturbation results for one planner, draw confidence interval.

- Code Path: AllocateR/alloc_discrete_paper_graphs/fs_opti_graph_multiperturb_202111.R
- Results Path: Dropbox/PrjNygaardSorensenWang/Results202111_perturb/capi_perturb/

### Results in Beta Variation Folder

We solve separately along a finer array of beta values for high and low education, and generate APC and MPC in aggregate across check amount to see economy-wide aggregate MPC/APC values. For the JEDC revision, we do this for the bush check problem

**Step 8a1**

- Code Path: PrjOptiSNW/invoke/202111/snw_res_b1_manna_mixture_bush.m
  > ls_fl_beta_val = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.97, 0.99];
- Results Path: Dropbox/PrjNygaardSorensenWang/Output202111/

**Step 8a2**

This is separate from *8a1* Generate low beta, high beta, low edu, high edu versions of the main file, not the beta-specific files.

- Code Path: AllocateR/jointbetedu/fs_wgt_merge_betaedutypes_202111.R
  > it_solve_jedc == 4

**Step 8b**

Generate MPC and APC Files and MPV and APV Files by household income bins, kids count, and marital status.

- Code Path: AllocateR/alloc_discrete_fun_R/fs_mpc_tables_increments_202111.R
- Results Path (in subfolders)
    + dropbox/PrjNygaardSorensenWang/Results202111/apc_mpc/
    + dropbox/PrjNygaardSorensenWang/Results202111_betaheter/apc_mpc/

> Run: ls_st_file_suffix_bushchck
> Run: ls_st_file_suffix_bushchck_betaedu

**Step 8c**

Generate MPC and APC Files and MPV and APV Files by household income bins, kids count, and marital status.

- Code Path: AllocateR/alloc_discrete_fun_R/fs_mpc_tables_aggregator_202111.R
  > ls_st_file_collection <- c('bushchck_main', 'bushchck_betaedu')
- Results Path
    + dropbox/PrjNygaardSorensenWang/Results202111/apc_mpc/
    + dropbox/PrjNygaardSorensenWang/Results202111_betaheter/apc_mpc/


## Allocation problem dimensions

This section contains some README information related to the dimensionality of the COVID related optimal allocation problems.

### The state-space of the Underlying Model

We solve for the value of 244 checks, each at $100, at these state-space points. We computed the marginal effects on consumption and value of 1st to 244th check for each element of the state-space:

- 2 marital groups: Marital Status 0 or 1
- 5 kids groups: Children 0 to 4
- 2 education groups: college or not
- 83 age groups: Age 18 to 100
- 65 savings groups
- 266 household head productivity shock.
- 5 spousal productivity shocks if married

Overall state-space: 86,104,200
  + married state-space: 5*2*83*65*266*5 = 71,753,500
  + single state-space (no spousal shock): 5*2*83*65*266*1 = 14,350,700

### 2020 Unemployment Shock COVID

There is a covid unemployment shock in 2020, that lowers (proportional) income realization.

So we solve the 86,104,200 twice once with the unemployment shock, once without

- State-space points total with COVID shock: 172,208,400

Since we have 244 checks, we solve for the consumption increase and value increase from $100 increment of the check 244 times, + 1 time for without check.

So overall, we have 42.2 billion MPCs: (244+1)*172208400=42,191,058,000

### How Many Cells We Solved for from 2019 Perspective for the COVID planner

The planner does not observe the full state-space so we intergrate from 2019 perspective given 2019 to 2020 COVID transition probability (conditional on state-space) and kids and shock transitions planner value.

- 5 kids groups: Children 0 to 4
- 2 marital groups: Marital Status 0 or 1
- 83 age groups: Age 18 to 100
- 509 Income Groups:
    + 476 bins below max actual phaseout: solved at $500 intervals between $0 and $238,000
    + 33 bins after max actual phaseout: solved at $5000 interval after $238,000, where the 33rd final bin is between $401,130 and Maximum
    + The income group is composed of individuals of with different 2019 productivity shocks and savings levels.

Together, these are:

- 5*2*83*509 = 422,470 groups/bins

We have (244+1) marginal average consumption gain, and value gain for each of the groups, so from 2019 planner perspective, we have 103.5 million MPCs

- 103,505,150 = 422470*(244+1)

### How Many Cells do we consider for the Allocation Problem

For the optimal allocation problem, we have several versions, how many individual types are there in each version?

- for FEASIBLE allocation, there are 970=5*2*97 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 97 income bins: the allocation planner sees approximately $2500 income bins between $0 and $238,000, and 1 bin after $238,000. There are 97 bins
- for OPTIMAL G4 (4 age groups 18 to 64) allocation, there are 3880=5*2*97*4 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 4 age groups
    + 97 income bins
- for OPTIMAL G47 (47 age groups) allocation, there are 45590=5*2*97*47 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 47 age groups
    + 97 income bins -->
