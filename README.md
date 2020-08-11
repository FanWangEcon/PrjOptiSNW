[![HitCount](http://hits.dwyl.io/fanwangecon/PrjOptiSNW.svg)](https://github.com/FanWangEcon/PrjOptiSNW)  [![Star](https://img.shields.io/github/stars/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/stargazers) [![Fork](https://img.shields.io/github/forks/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/network/members) [![Star](https://img.shields.io/github/watchers/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/watchers)

This is a work-in-progress Matlab package consisting of functions that solve the dynamic life cycle model in Nygaard, Sorensen and Wang (2020). The code companion presents solutions to the dynamic life-cycle problem, and methods for evaluating the marginal gains from allocating additional welfare checks. Tested with [Matlab](https://www.mathworks.com/products/matlab.html) 2019a [@matlab].

All functions are parts of a matlab toolbox that can be installed:

> Download and install the Matlab toolbox: [PrjOptiSNW.mltbx](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW.mltbx)

This bookdown file is a collection of mlx based vignettes for functions that are available from [PrjOptiSNW](https://github.com/FanWangEcon/PrjOptiSNW). Each Vignette file contains various examples for invoking each function. The goal of this repository is to make it easier to find/re-use codes produced for various projects.

The package relies on [MEconTools](https://fanwangecon.github.io/MEconTools/), which needs to be installed first. The package does not include allocation functions, only simulation code to generate the value of each welfare check increments for households. Allocation functions rely the R optimal allocation package [PrjOptiAlloc](https://fanwangecon.github.io/PrjOptiAlloc).

From other repositories: For dynamic borrowing and savings problems, see [Dynamic Asset Repository](https://fanwangecon.github.io/CodeDynaAsset/); For code examples, see also [R Example Code](https://fanwangecon.github.io/R4Econ/), [Matlab Example Code](https://fanwangecon.github.io/M4Econ/), and [Stata Example Code](https://fanwangecon.github.io/Stata4Econ/); For intro stat with R, see [Intro Statistics for Undergraduates](https://fanwangecon.github.io/Stat4Econ/), and intro Math with Matlab, see [Intro Mathematics for Economists](https://fanwangecon.github.io/Math4Econ/). See [here](https://github.com/FanWangEcon) for all of [Fan](https://fanwangecon.github.io/)'s public repositories.

Please contact [FanWangEcon](https://fanwangecon.github.io/) for issues or problems.

[![](https://img.shields.io/github/last-commit/fanwangecon/PrjOptiSNW)](https://github.com/FanWangEcon/PrjOptiSNW/commits/master) [![](https://img.shields.io/github/commit-activity/m/fanwangecon/PrjOptiSNW)](https://github.com/FanWangEcon/PrjOptiSNW/graphs/commit-activity) [![](https://img.shields.io/github/issues/fanwangecon/PrjOptiSNW)](https://github.com/FanWangEcon/PrjOptiSNW/issues) [![](https://img.shields.io/github/issues-pr/fanwangecon/PrjOptiSNW)](https://github.com/FanWangEcon/PrjOptiSNW/pulls)

# Installation

In addition to downloading and installing [PrjOptiSNW.mltbx](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW.mltbx), can also:

```
# Clone Package from Git Bash
cd "C:/Downloads"
git clone https://github.com/fanwangecon/PrjOptiSNW.git
```

Install the Package from inside Matlab:

```
# Install Matlab Toolbox PrjOptiSNW
toolboxFile = 'C:/Downloads/PrjOptiSNW/PrjOptiSNW.mltbx';
# toolboxFile = 'C:/Users/fan/PrjOptiSNW/PrjOptiSNW.mltbx';
agreeToLicense = true;
installedToolbox = matlab.addons.toolbox.installToolbox(toolboxFile, agreeToLicense)
```

# 1  Parameters

1. [Model Parameters](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_param.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/snwx_mp_param.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_param.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_param.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_param.html)
	+ Model parameters, transition matrices, permanent heterogeneities.
	+ **PrjOptiSNW**: *[snw_mp_param()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/params/snw_mp_param.m)*
2. [Model Controls Parameters](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_control.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/snwx_mp_control.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_control.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_control.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/params/htmlpdfm/snwx_mp_control.html)
	+ Parameters to control display options etc.
	+ **PrjOptiSNW**: *[snw_mp_control()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/params/snw_mp_control.m)*

# 2  Solving the Dynamic Life Cycle Problem

1. [Policy and Value Functions Dynamic Life Cycle Vectorized Bisection](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpol/htmlpdfm/snwx_vfi_bisec_vec.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpol/snwx_vfi_bisec_vec.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpol/htmlpdfm/snwx_vfi_bisec_vec.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpol/htmlpdfm/snwx_vfi_bisec_vec.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpol/htmlpdfm/snwx_vfi_bisec_vec.html)
	+ Solving for policy and value functions from 18 to 100 years of age, at 1 year interval.
	+ Households face persistent productivity shocks for household heads, stochastic shocks for spousal income, exogenous children under age 17 transition probability, and age-specific household-head survival probabilities.
	+ The household can have up to four children under age 17, and has permanent heterogeneity in marital status and education types.
	+ Problem solved for exact savings choices using [vectorized bisection](https://fanwangecon.github.io/MEconTools/MEconTools/doc/vfi/htmlpdfm/fx_vfi_az_bisec_vec.html) from from [MEconTools](https://fanwangecon.github.io/MEconTools/).
	+ **PrjOptiSNW**: *[snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m)*

# 3  Alternative Value Function Solution Testing

1. [Small Test Looped Minimizer Routine to Solve for Exact Savings Choices](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/snwx_vfi_test.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test.html)
	+ Solve for the exact savings choices using matlab minimizer in an iterative loop.
	+ The code demonstrates the solution structure. We use [snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m) with [vectorized bisection](https://fanwangecon.github.io/MEconTools/MEconTools/doc/vfi/htmlpdfm/fx_vfi_az_bisec_vec.html) for working implementations.
	+ Due to speed, only show testing results at small grid without spousal shocks.
	+ **PrjOptiSNW**: *[snw_vfi_main()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main.m)*
2. [Small Test Looped over States Grid Search Solution](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_grid_search.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/snwx_vfi_test_grid_search.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_grid_search.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_grid_search.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_grid_search.html)
	+ The savings choice grid is the same as the savings states grid. Solve for optimal savings choices using grid-search. Loop over the state space, at each state-space point, vectorized optimization.
	+ Our problem requires very high precision to solve for the marginal gains to households from each increment of welfare checks. We rely on the exact solution method from [snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m) for the working code.
	+ Due to speed, only show testing results at small grid without spousal shocks.
	+ **PrjOptiSNW**: *[snw_vfi_main_grid_search()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_grid_search.m)*
3. [Small Test Vectorized Bisection Solve for Exact Savings Choices](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/snwx_vfi_test_bisec_vec.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec.html)
	+ Vectorized bisection exact solution code tested with small grid to compare to alternative solutiom methods.
	+ Small grid without spousal shocks.
	+ **PrjOptiSNW**: *[snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m)*
4. [Small Test Spousal Shocks  Test Vectorized Bisection Solve for Exact Savings Choices](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec_spousalshock.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/snwx_vfi_test_bisec_vec_spousalshock.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec_spousalshock.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec_spousalshock.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolsmall/htmlpdfm/snwx_vfi_test_bisec_vec_spousalshock.html)
	+ Vectorized bisection exact solution code tested with small grid to compare to alternative solutiom methods.
	+ Small grid with spousal shocks. There are three shocks: persistent household head income shock, i.i.d. spousal income shock, and persistent kids count transition shocks.
	+ **PrjOptiSNW**: *[snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m)*

# 4  Solution with Unemployment

1. [Policy and Value Functions Dynamic Life Cycle if Unemployed](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolunemploy/htmlpdfm/snwx_vfi_unemp_bisec_vec.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolunemploy/snwx_vfi_unemp_bisec_vec.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolunemploy/htmlpdfm/snwx_vfi_unemp_bisec_vec.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/svalpolunemploy/htmlpdfm/snwx_vfi_unemp_bisec_vec.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/svalpolunemploy/htmlpdfm/snwx_vfi_unemp_bisec_vec.html)
	+ Solving the dynamic programming problem conditional on having an one period unemployment shock.
	+ There is an unemployment shock in 2020. We first solve for the policy and value functions without the unemployment shock.
	+ Using the value function from the world without the 2020 covid unemployment shock as future values, we solve for optimal choices in 2020 given a COVID unemployment shock.
	+ The COVID shock lowers the realization of household's stochastic income process proportionally, but the lost income might be replenished by unemployment benefits up to 100 percent. Unemployment benefits have to be paid for by taxes.
	+ **PrjOptiSNW**: *[snwx_vfi_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/svalpol/snw_vfi_main_bisec_vec.m)*

# 5  Household Life Cycle Distribution

1. [Assets and Demographic Distributions with Continuous Exact Savings Choices](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_bisec_vec_loop.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/snwx_ds_bisec_vec_loop.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_bisec_vec_loop.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_bisec_vec_loop.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_bisec_vec_loop.html)
	+ Simulate the life cycle distribution of assets, consumptions, and demographic patterns up to age 100, given exogenous initial distributions at age 18.
	+ Solves for budget clearing tax rates given distributional results.
	+ Uses vectorized bisection to solve for exact savings choices, looped distribution code.
	+ **PrjOptiSNW**: *[snw_ds_main()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main.m)*
2. [Assets and Demographic Distributions with Grid Search](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_grid_search.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/snwx_ds_grid_search.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_grid_search.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_grid_search.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/sdist/htmlpdfm/snwx_ds_grid_search.html)
	+ Grid search solution using grid search for savings choices, the savings state-space grid is the same as the savings choice-grid.
	+ Exact choice solution from [snw_ds_main()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main.m) generates significantly smoother distributions.
	+ **PrjOptiSNW**: *[snw_ds_main_grid_search()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/sdist/snw_ds_main_grid_search.m)*

# 6  Value of Each Check

1. [Marginal Gain Per Check 2020 Employed](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_wrk_bisec_vec.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/snwx_a4chk_wrk_bisec_vec.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_wrk_bisec_vec.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_wrk_bisec_vec.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_wrk_bisec_vec.html)
	+ Evaluate the marginal gain per check in 2020 if household head is employed.
	+ Solve for the increase in savings that is equivalent to the impact of an additional check on a household's resource available in 2020, given tax and interest rates considerations.
	+ **PrjOptiSNW**: *[snw_a4chk_wrk_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanner/snw_a4chk_wrk_bisec_vec.m)*
2. [Marginal Gain Per Check 2020 Unemployed](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_unemp_bisec_vec.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/snwx_a4chk_unemp_bisec_vec.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_unemp_bisec_vec.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_unemp_bisec_vec.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannercheckval/htmlpdfm/snwx_a4chk_unemp_bisec_vec.html)
	+ Evaluate the marginal gain per check in 2020 if household head is unemployed.
	+ Solve for the increase in savings that is equivalent to the impact of an additional check on a household's resource available in 2020, given tax and interest rates considerations.
	+ **PrjOptiSNW**: *[snw_a4chk_unemp_bisec_vec()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanner/snw_a4chk_unemp_bisec_vec.m)*

# 7  2020 Outcomes Full State Space with Savings, Shocks and Education

1. [Value in 2020 Given Age, Savings, Shocks, Kids, Educaiton and Marriage](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw20_jaeemk.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/snwx_evuvw20_jaeemk.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw20_jaeemk.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw20_jaeemk.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw20_jaeemk.html)
	+ Expected value and expected consumption from 2020 for a household given at a particular age (18-100), with a particular savings level, at a particular combination of household head and spouse income shocks, with 0 to 4 children, high or low Education status, and married or not married.
	+ This uses the unemployment probability and generates the average value given the probability of the unemployment state that is dependent on the state-space.
	+ **PrjOptiSNW**: *[snw_evuvw20_jaeemk()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw20_jaeemk.m)*
2. [Expected Value in 2019 Given Age, Savings, Shocks, Kids, Educaiton and Marriage](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw19_jaeemk.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/snwx_evuvw19_jaeemk.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw19_jaeemk.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw19_jaeemk.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjaeemk/htmlpdfm/snwx_evuvw19_jaeemk.html)
	+ Expected value and expected consumption from 2019 for a household at a particular age (18-99), savings level, shocks combinations, kids/education/marriage status, given 2019 optimal savings choices, income shock transition probability as well as household children count transition probabilities.
	+ **PrjOptiSNW**: *[snw_evuvw19_jaeemk()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw19_jaeemk.m)*

# 8  2019 Expectations Given Income, Age, Kids and Marital Status

1. [Expected Value from 2019 Given Age, Kids, Income and Marriage](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/snwx_evuvw19_jmky.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky.html)
	+ Expected Value from 2019 Given Age, Kids, Income and Marriage.
	+ Each 2019 income group consists of individuals with varying productivity shocks, savings, and from lower and higher education groups.
	+ **PrjOptiSNW**: *[snw_evuvw19_jmky()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw19_jmky.m)*
2. [Expected Value from 2019 Given Age, Kids, Income and Marriage for All Checks](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky_allchecks.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/snwx_evuvw19_jmky_allchecks.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky_allchecks.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky_allchecks.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/splannerjmky/htmlpdfm/snwx_evuvw19_jmky_allchecks.html)
	+ Expected Value from 2019 Given Age, Kids, Income and Marriage for All Checks.
	+ This is the gateway function that solves policy functions, derive distributions, computes value in 2020 with and without unemployment shocks with varying check levels, derives 2019 planner expected values given household optimization and shocks, and finds the mass of individuals in different income/age/marital-status bins, and saves the simulated value of check results for the planner.
	+ **PrjOptiSNW**: *[snw_evuvw19_jmky_allchecks()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/splanneralt/snw_evuvw19_jmky_allchecks.m)*

# 9  Calibration

1. [Calibrate Discount Factor and Normalize GDP](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/calibrate/htmlpdfm/snwx_calibrate_beta_norm_gdp.html): [**mlx**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/calibrate/snwx_calibrate_beta_norm_gdp.mlx) \| [**m**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/calibrate/htmlpdfm/snwx_calibrate_beta_norm_gdp.m) \| [**pdf**](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/doc/calibrate/htmlpdfm/snwx_calibrate_beta_norm_gdp.pdf) \| [**html**](https://fanwangecon.github.io/PrjOptiSNW/PrjOptiSNW/doc/calibrate/htmlpdfm/snwx_calibrate_beta_norm_gdp.html)
	+ We calibrate the model so that the Asset/Savings/Capital to GDP/Income ratio is 3.
	+ We normalize the model so that median household income is equal to 1 in the model.
	+ **PrjOptiSNW**: *[snw_calibrate_beta_norm_gdp()](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/calibrate/snw_calibrate_beta_norm_gdp.m)*

----
Please contact [![](https://img.shields.io/github/followers/fanwangecon?label=FanWangEcon&style=social)](https://github.com/FanWangEcon) [![](https://img.shields.io/twitter/follow/fanwangecon?label=%20&style=social)](https://twitter.com/fanwangecon) for issues or problems.

![RepoSize](https://img.shields.io/github/repo-size/fanwangecon/PrjOptiSNW)
![CodeSize](https://img.shields.io/github/languages/code-size/fanwangecon/PrjOptiSNW)
![Language](https://img.shields.io/github/languages/top/fanwangecon/PrjOptiSNW)
![Release](https://img.shields.io/github/downloads/fanwangecon/PrjOptiSNW/total)
![License](https://img.shields.io/github/license/fanwangecon/PrjOptiSNW)

