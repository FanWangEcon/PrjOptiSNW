Feasible, Optimal and Threshold allocation Line by Line Testing Files, created
structures which were used to create:
*PrjOptiAlloc/R/ffp_snw_welfarechecks_graph.R* and
*PrjOptiAlloc/R/ffp_snw_welfarechecks_input.R*

Then in the alloc_discrete_fun folder, the functionalized versions of feasible,
optimal and Threshold files solve various allocation problems and generating
graphs by calling the graph and input functions above.

Continued using the *fs_opti_feasible.Rmd* file for developing functions and
testing before functionalities are added to the PrjOptiAlloc folder functions.

THe optimal and threshold files were developed initially, and are out of date
with respect to the optimal and threshold allocation results from the
*alloc_discte_fun* folder.

Note that:

- **optimal**: allocate considering age, income, kids count and marriage status
- **feasible**: allocate considering income, kids count and marriage status
- **threshold**: allocate given kids/marriage type specific upper bounds, 1200 for
  unmarried one child household.
