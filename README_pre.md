[![HitCount](http://hits.dwyl.io/fanwangecon/PrjOptiSNW.svg)](https://github.com/FanWangEcon/PrjOptiSNW)  [![Star](https://img.shields.io/github/stars/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/stargazers) [![Fork](https://img.shields.io/github/forks/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/network/members) [![Star](https://img.shields.io/github/watchers/fanwangecon/PrjOptiSNW?style=social)](https://github.com/FanWangEcon/PrjOptiSNW/watchers)

This is a work-in-progress Matlab package consisting of functions that solve the dynamic life cycle model in [Nygård](https://sites.google.com/site/vegardmokleivnygaard/), [Sørensen](https://uh.edu/~bsorense/) and [Wang](https://fanwangecon.github.io/) (2020). The code companion presents solutions to the dynamic life-cycle problem, and methods for evaluating the marginal gains from allocating additional welfare checks. Tested with [Matlab](https://www.mathworks.com/products/matlab.html) 2019a.

All functions are parts of a matlab toolbox that can be installed:

> Download and install the Matlab toolbox: [PrjOptiSNW.mltbx](https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW.mltbx)

This bookdown file is a collection of mlx based vignettes for functions that are available from [PrjOptiSNW](https://github.com/FanWangEcon/PrjOptiSNW). Each Vignette file contains various examples for invoking each function.

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
