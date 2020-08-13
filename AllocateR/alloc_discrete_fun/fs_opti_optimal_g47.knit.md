---
title: |
  | Optimal Allocation
  | Conditional on Age (47 Age Groups One by One), Income, Marital Status and Kids Count
  | moredense_a100zh266_e2m2 vs moredense_a65zh266zs5_e2m2_b0_calibrated
  | ybin = 2500

output:
  html_notebook:
    toc: yes
    toc_depth: 4
---



# Optimal Allocation

2019 Age, Income, Kids Count and Martial Status (non-stochastic) based optimal linear allocation (Utilitarian). Will solve for optimal allocation results given different simulation structures:

Each age one by one, to get the fully optimal age specific allocation. Do not graph.

Files:


```r
# # File Names
# st_file_type_withspouse_shock <- 'moredense_ybin25000'
# snm_simu_csv_withspouse_shock <- paste0('snwx_v_planner_',st_file_type_withspouse_shock,'.csv')

# # File Names
# st_file_type_withspouse_shock <- 'dense_ybin2500'
# snm_simu_csv_withspouse_shock <- paste0('snwx_v_planner_',st_file_type_withspouse_shock,'.csv')

# File Names
# st_file_type_withspouse_shock <- 'moredense_a100zh81zs5_e2m2'
# st_file_type_withspouse_shock <- 'moredense_a65zh133zs5_e2m2'
# st_file_type_withspouse_shock <- 'moredense_a65zh266zs5_e2m2'
# st_file_type_withspouse_shock <- 'moredense_a65zh266zs5_e2m2_b0_calibrated'
# st_file_type_withspouse_shock <- 'moredense_a65zh266zs5_e2m2_b1_calibrated'
# snm_simu_csv_withspouse_shock <- paste0('snwx_v_planner_',st_file_type_withspouse_shock,'.csv')

ls_output <- fs_opti_support()
st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
srt_simu_path <- ls_output$srt_simu_path
bl_save_img <- ls_output$bl_save_img
spt_img_save <- ls_output$spt_img_save
srt_csv_path <- ls_output$srt_csv_path
```

## Common Parameters Across


```r
# Max Phase Out given 1200*2 + 500*4 = 4400
fl_max_phaseout = 238000
it_bin_dollar_before_phaseout = 2500
# Dollar Per Check
fl_percheck_dollar = 100
# Meaning of Ymin Ymax simulated interval of 1
fl_multiple = 58056
# Number of Max Checks
it_max_checks_1st = 44
it_max_checks_2nd = 88
# Number of Tax Paying Households
fl_tax_hh = 128580000
# Number of Income Groups to Use: use 25 for 10,000 = 1
# Age Conditions
# it_max_age = 64
# it_min_age = 64
it_max_age = 64
it_min_age = 18
it_age_bins = 47
# Image Save Suffix
st_img_suf_age_ybin <- paste0(it_min_age, 't', it_max_age)
```

### Variable Names and Paths


```r
# File Path
# srt_simu_path <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Output/'

# Column Names
ar_svr_csv <- c('age', 'marital', 'kids', 'checks',	'ymin', 'mass', 'survive', 'vtilde', 'ctilde')
# Variables That Identify Individual Types
ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
ar_svr_groups_noage <- c('marital', 'kids', 'ymin_group')
ar_svr_groups_stats <- c('mass', 'survive')
# Number of Checks and Planner Value
svr_checks <- 'checks'
svr_v_value <- 'vtilde'
svr_c_value <- 'ctilde'
svr_mass <- 'mass'
```

Image and Save Control


```r
# Save Folder
# srt_csv_path <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Results/2020-08-05/csv/'
# # CSV and Image Paths
# srt_csv_path = paste0(srt_csv_path, st_file_type_withspouse_shock,'/')
# dir.create(file.path(srt_csv_path), showWarnings = FALSE, recursive = TRUE)
```

### Planner Preference


```r
ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
ar_rho <- unique(c(1,ar_rho))
# ar_rho <- c(1)
```

## Round 1 Results with Spousal Income Shocks

### Process

Call the input processing function


```r
# Call function
ls_prc_outputs_zs5_1st <- PrjOptiAlloc::ffp_snw_process_inputs(
  srt_simu_path = srt_simu_path,
  snm_simu_csv = snm_simu_csv_withspouse_shock,
  fl_max_phaseout = fl_max_phaseout,
  it_bin_dollar_before_phaseout = it_bin_dollar_before_phaseout,
  fl_percheck_dollar = fl_percheck_dollar,
  fl_multiple = fl_multiple,
  it_max_checks = it_max_checks_1st,
  fl_tax_hh = fl_tax_hh,
  it_max_age = it_max_age,
  it_min_age = it_min_age,
  it_age_bins = it_age_bins,
  ar_svr_csv = ar_svr_csv,
  ar_svr_groups = ar_svr_groups,
  ar_svr_groups_stats = ar_svr_groups_stats,
  svr_checks = svr_checks,
  svr_v_value = svr_v_value,
  svr_c_value = svr_c_value,
  svr_mass = svr_mass,
  ar_rho = ar_rho,
  bl_non_inc_adjust = TRUE,
  bl_print = FALSE,
  bl_print_verbose = FALSE)
```

```
## Warning: `cols` is now required when using unnest().
## Please use `cols = c(rev)`

## Warning: `cols` is now required when using unnest().
## Please use `cols = c(rev)`
```

```
## Adding missing grouping variables: `id_i`
```

### Outputs

REV:


```r
tb_rho_rev_c=ls_prc_outputs_zs5_1st$tb_rho_rev_c
tb_rho_rev_v=ls_prc_outputs_zs5_1st$tb_rho_rev_v
print(tb_rho_rev_c)
```

```
## # A tibble: 9 x 3
##     rho rho_val   REV
##   <int>   <dbl> <dbl>
## 1     1   1     0.358
## 2     2   0.99  0.360
## 3     3   0.963 0.365
## 4     4   0.861 0.386
## 5     5   0.482 0.487
## 6     6  -0.931 0.864
## 7     7  -6.20  0.999
## 8     8 -25.8   1.00 
## 9     9 -99     1.00
```

```r
print(tb_rho_rev_v)
```

```
## # A tibble: 9 x 3
##     rho rho_val   REV
##   <int>   <dbl> <dbl>
## 1     1   1     0.804
## 2     2   0.99  0.805
## 3     3   0.963 0.806
## 4     4   0.861 0.810
## 5     5   0.482 0.824
## 6     6  -0.931 0.876
## 7     7  -6.20  1.00 
## 8     8 -25.8   1.00 
## 9     9 -99     1.00
```

## Round 2 Results with Spousal Income Shocks


```r
bl_given_firstcheck <- TRUE
```

### Process

Call the input processing function


```r
# Call function
ls_prc_outputs_zs5_2nd <- PrjOptiAlloc::ffp_snw_process_inputs(
  srt_simu_path = srt_simu_path,
  snm_simu_csv = snm_simu_csv_withspouse_shock,
  fl_max_phaseout = fl_max_phaseout,
  it_bin_dollar_before_phaseout = it_bin_dollar_before_phaseout,
  fl_percheck_dollar = fl_percheck_dollar,
  fl_multiple = fl_multiple,
  it_max_checks = it_max_checks_2nd,
  fl_tax_hh = fl_tax_hh,
  it_max_age = it_max_age,
  it_min_age = it_min_age,
  it_age_bins = it_age_bins,
  ar_svr_csv = ar_svr_csv,
  ar_svr_groups = ar_svr_groups,
  ar_svr_groups_stats = ar_svr_groups_stats,
  svr_checks = svr_checks,
  svr_v_value = svr_v_value,
  svr_c_value = svr_c_value,
  svr_mass = svr_mass,
  ar_rho = ar_rho,
  bl_given_firstcheck = bl_given_firstcheck,
  bl_non_inc_adjust = TRUE,
  bl_print = FALSE,
  bl_print_verbose = FALSE)
```

```
## Warning: `cols` is now required when using unnest().
## Please use `cols = c(rev)`

## Warning: `cols` is now required when using unnest().
## Please use `cols = c(rev)`
```

```
## Adding missing grouping variables: `id_i`
```

### Outputs

REV:


```r
tb_rho_rev_c=ls_prc_outputs_zs5_2nd$tb_rho_rev_c
tb_rho_rev_v=ls_prc_outputs_zs5_2nd$tb_rho_rev_v
print(tb_rho_rev_c)
```

```
## # A tibble: 9 x 3
##     rho rho_val   REV
##   <int>   <dbl> <dbl>
## 1     1   1     0.337
## 2     2   0.99  0.339
## 3     3   0.963 0.343
## 4     4   0.861 0.362
## 5     5   0.482 0.455
## 6     6  -0.931 0.860
## 7     7  -6.20  1.00 
## 8     8 -25.8   1.00 
## 9     9 -99     1.00
```

```r
print(tb_rho_rev_v)
```

```
## # A tibble: 9 x 3
##     rho rho_val   REV
##   <int>   <dbl> <dbl>
## 1     1   1     0.799
## 2     2   0.99  0.801
## 3     3   0.963 0.802
## 4     4   0.861 0.806
## 5     5   0.482 0.819
## 6     6  -0.931 0.870
## 7     7  -6.20  1.00 
## 8     8 -25.8   1.00 
## 9     9 -99     1.00
```

## Print Results to CSV

Only print a limited set of results, otherwise there would be too many rows.

### Print Allocation Results

Print Allocation Results in Two Ways, Row based, with actual, first, and second round allocations as columns. And this table is outputed, at different levels of aggregations. where the allocation columns are outputed as different weighted means. Also compute fraction 0, fraction at within group max.








