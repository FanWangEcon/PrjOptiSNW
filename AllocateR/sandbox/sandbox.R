ffi_vox_checks_bidencheck <- function(ymin_group, marital, kids,
                                      ar_ycut,
                                      fl_multiple = 58056, fl_percheck_dollar = 100,
                                      it_inc_subgroups = 10) {
  # if policy change, also go modify section *Modify Maximum Allocation Point for Each Individual* in main above
  # ar_ycut_usd <- c(0, 20000, 40000, 60000, 80000, 100000, 100000000)
  # ar_ycut <- ar_ycut_usd/fl_multiple
  # fl_multiple <- 58056
  # fl_percheck_dollar <- 100
  # # Poorest married 4 kids
  # ffi_vox_checks_bidencheck(1,1,4, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_bidencheck(2,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(2,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(2,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(2,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_bidencheck(4,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(4,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(4,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(4,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_bidencheck(5,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(5,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(5,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(6,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_bidencheck(11,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(11,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(11,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_bidencheck(11,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)

  # income minimum and maximum,
  # by construction: length(ar_ycut_dollar) = max(ymin_group) + 1
  ar_ycut_dollar = ar_ycut*fl_multiple
  it_ygroups = length(ar_ycut_dollar)

  # Default no checks
  fl_check_dollor = 0

  # Only provide checks if not in the last income group, where all receive none
  if (ymin_group < it_ygroups - 1) {

    fl_slope_m0k0 = 1400/(100000-75000)
    fl_slope_m0k1 = 2800/(150000-112500)
    fl_slope_m0k2 = 4200/(150000-112500)
    fl_slope_m0k3 = 5600/(150000-112500)
    fl_slope_m0k4 = 7000/(150000-112500)

    fl_slope_m1k0 = 2800/(200000-150000)
    fl_slope_m1k1 = 4200/(200000-150000)
    fl_slope_m1k2 = 5600/(200000-150000)
    fl_slope_m1k3 = 7000/(200000-150000)
    fl_slope_m1k4 = 8400/(200000-150000)

    # start point household head
    fl_check_dollar = 1400

    # married household gets more, marital = 0, 1
    fl_check_dollar = fl_check_dollar + marital*1400

    # Households with kids: 0, 1,2,3,4
    fl_check_dollar = fl_check_dollar + kids*1400


    # lower and upper bounds on income
    fl_inc_group_lower_bound = ar_ycut_dollar[ymin_group]
    fl_inc_group_upper_bound = ar_ycut_dollar[ymin_group+1]

    # A grid of income between these two points: 10 points
    ar_inc_grid = seq(fl_inc_group_lower_bound, fl_inc_group_upper_bound,
                      length.out=it_inc_subgroups)

    # What is the tax rate at each point of these incomes given marry and kids?
    ar_check_reduce_inc_grid = matrix(data=NA, nrow=length(ar_inc_grid), ncol=1)

    # as income increases, fl_check_dollar go down
    it_ctr = 0
    for (fl_inc in ar_inc_grid) {

      it_ctr = it_ctr + 1

      # phaseout starts $75,000 for single
      # phaseout starts $112,500 for heads of household unmarried
      # phaseout starts $150,000 for heads of household married
      if (marital == 0 && kids == 0) {
        fl_check_reduce = ((max(fl_inc - 75000,0))*fl_slope_m0k0)
      } else if (marital == 0 && kids == 1) {
        fl_check_reduce = ((max(fl_inc - 112500,0))*fl_slope_m0k1)
      } else if (marital == 0 && kids == 2) {
        fl_check_reduce = ((max(fl_inc - 112500,0))*fl_slope_m0k2)
      } else if (marital == 0 && kids == 3) {
        fl_check_reduce = ((max(fl_inc - 112500,0))*fl_slope_m0k3)
      } else if (marital == 0 && kids == 4) {
        fl_check_reduce = ((max(fl_inc - 112500,0))*fl_slope_m0k4)
      } else if (marital == 1 && kids == 0) {
        fl_check_reduce = ((max(fl_inc - 150000,0))*fl_slope_m1k0)
      } else if (marital == 1 && kids == 1) {
        fl_check_reduce = ((max(fl_inc - 150000,0))*fl_slope_m1k1)
      } else if (marital == 1 && kids == 2) {
        fl_check_reduce = ((max(fl_inc - 150000,0))*fl_slope_m1k2)
      } else if (marital == 1 && kids == 3) {
        fl_check_reduce = ((max(fl_inc - 150000,0))*fl_slope_m1k3)
      } else if (marital == 1 && kids == 4) {
        fl_check_reduce = ((max(fl_inc - 150000,0))*fl_slope_m1k4)
      }

      ar_check_reduce_inc_grid[it_ctr] = max(0, fl_check_dollar - fl_check_reduce)

    }

    fl_check_dollor = mean(ar_check_reduce_inc_grid)

  }

  # Check Numbers
  fl_avg_checks = round(fl_check_dollor/fl_percheck_dollar, 0)

  return(fl_avg_checks)
}
