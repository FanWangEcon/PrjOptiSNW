# Distribution

Six *snwx_ds_main* files. Vignettes for Simulate distribution given exogenous initial distributions of observables states, and stationary productivity shock distributions, etc.

Three types of Distributional Simulations:

1. small
2. dense
3. more-dense

Simulated either with grid_search for both VFI and distribution. Or with Bisec_vec and looped distribution.

## Distribution Speed Comparison

In Seconds (speeds include both VFi and DS?)

For bisec_vec_loop,

1. small: 2.557144
2. dense: 211.000592
3. more-dense: 2035.819094

For grid_search,

1. small: 8.496160
2. dense: 704.479941
3. more-dense: 10961.929988

## Distribtuion Outcomes Comparison

### Interest Share of Income

For bisec_vec_loop

1. small: 0.050003
2. dense: 0.119
3. more-dense: 0.12113

For grid_search:

1. small: 0.039337
2. dense: 0.093803
3. more-dense: 0.11794

### Consumption coefficient of variation

For bisec_vec_loop

1. small: 0.47406
2. dense: 0.62317
3. more-dense: 0.62567

For grid_search:

1. small: 0.49183
2. dense: 0.68246
3. more-dense: 0.63205

### Savings held by bottom 90 percentile

look at a_ss (not ap_ss), 1 minus is savings held by top 10 percentile.

For bisec_vec_loop

1. small: 0.61665
2. dense: 0.52742
3. more-dense: 0.54182

For grid_search:

1. small: 0.39457
2. dense: 0.51613
3. more-dense: 0.52491
