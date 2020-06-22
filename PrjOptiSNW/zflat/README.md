Main.m is the main file. Running this file runs all the files needed to compute value functions, aggregation, calibration, and planner problems

Rouwenhorst.m and mkv_recursive.m are used to discretize AR(1) productivity shocks.

VFI.m solves the household's optimization problem by means of backwards induction

utility.m computes flow utility given consumption, marital status, and number of children.

Tax.m computes income taxes.

spousal_income.m computes income from spouse in the event that the individual is married.

value_func_aux.m is called within VFI.m. It computes the value function under alternative values of assets for next period. The Matlab function fmincon is used to derive the optimal value of assets for next period.

Aggregation.m uses the transition probabilities and the policy functions derived in VFI.m to compute the stationary distribution across idiosyncratic states and to compute aggregate wealth and aggregate income. It also adjusts the average tax burden to ensure that the government balances its budget.

VFI_unemp.m computes the value function in the event of unemployment. It derives policy functions for consumption and saving in 2020, but uses the value function derived in VFI.m to compute continuation values.

V_working_proxy.m computes the value of an individual who is working in 2020 and receiving a total number of "welf_checks" welfare checks, each of which is worth "TR" dollars.

V_unemp_proxy.m computes the value of an individual who is unemployed in 2020 and receiving a total number of "welf_checks" welfare checks, each of which is worth "TR" dollars.

find_a_working.m is called within V_working_proxy.m to find the asset value on the asset grid that approximates the total value of welfare checks the individual receives.

find_a_unemp.m is called within V_unemp_proxy.m to find the asset value on the asset grid that approximates the total value of welfare checks the individual receives.

pi_unemp_calibration.m finds the factor by which unemployment probabilities in 2020 from Cajner et al. (2020) have to be adjusted to match both unemployment probabilities by wage quintiles and unemployment probabilities by age groups.

Planner.m computes the value of the planner given age, marital status, number of children, and income in 2019.