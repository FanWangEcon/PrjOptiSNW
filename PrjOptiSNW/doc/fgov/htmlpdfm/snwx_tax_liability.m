%% Taxable Income and Tax Liabilities in 2008
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/fgov/G:\repos\PrjOptiSNW\PrjOptiSNW\doc\fgov\snw_tax_liability.m 
% *snw_tax_liability*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*,> the function solves for tax liability.
% 
% We can study the effects of the 2008 Tax Rebate. The Tax rebate is a rebate 
% based on how much tax was paid, so we need to know taxable income and tax liability. 
% These differ by income, household marital status, and the count of children. 
% Given an array of pre-tax income values, we compute for from 0 to 4 kids and 
% both married and unmarried taxable-income and tax-liability at all points along 
% the income array. Deductions are from 2008 income (<https://www.irs.gov/pub/irs-prior/f1040--2008.pdf 
% 2008 IRS f1040>). Tax brackets from 2008 are here (<https://www.irs.gov/pub/irs-prior/i1040tt--2008.pdf 
% 2008 IRS 1040tt>).
%% Taxable Income and Tax Liabilities in 2008 for 4 Income Levels
% Call the function at four income levels. Solve for different kids count and 
% by marital status.

bl_visualize = true;
ar_income = [20000, 30000, 50000, 70000];
snw_tax_liability(ar_income, true);
%% Taxable Income and Tax Liabilities in 2008 for 50 Income Levels
% Call the function for incomes from 0 to 200k. Solve for different kids count 
% and by marital status.

bl_visualize = true;
ar_income = linspace(0, 2e5, 50);
snw_tax_liability(ar_income, true);