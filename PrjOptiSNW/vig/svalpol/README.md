# VFI Simulations

Vignettes for Solve for value and policy functions. 8 *snwx_vfi_main* files for VFI, 5 additional VFI files for unemployment *snwx_vfi_unemp*

Three types of VFI simulations:

1. small
2. dense
3. more-dense

For small, four types of simulations:

1. fmincon, 1 shock spouse, 5 shock husband
2. grid_search, 1 shock spouse, 5 shock husband
3. bisec_vec, 3 shock spouse, 5 shock husband
4. bisec_vec, 3 shock spouse, 5 shock husband

For dense and more dense, both bisec_vec and grid_search results

## VFI Speed comparison

1. dense grid_search: 516.055576
2. dense bisec_vec: 108.279731
3. more-dense grid_search: 10704.736019
4. more-dense bisec_vec: 1167.364841

## VFI Unemployment

- fmincon, grid_search and bisec_vec for small.
- grid_search and bisec_vec for dense.

Show difference between value under employment shock and no employment shock. grid_search results possibly have lower consumption under employment shock.
