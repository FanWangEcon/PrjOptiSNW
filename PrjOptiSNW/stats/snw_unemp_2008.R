# 2021-11-24 07:55
#
# Find the P(U|Age,Edu). Steps:
#
# 1. Given the age death rate, what is the marginal distribution by Age, and then the joint P(Age,Edu) distribution: P(Age,Edu) = P(Age)*P(Edu), the edu distribution is independent of age, unrelated to death probability.
# 2. Get the conditional unemployment probability by Age and Edu (Vegard provided)
# 3. Use the R-Functions with rectilinear restriction to solve


# A. P(Age, Edu) Distribution -----------------

# Modify C:\Users\fan\PrjOptiSNW\PrjOptiSNW\doc\sdist\snwx_ds_bisec_vec.mlx file to run DENSE:
# We find P(Age) are:
ar_p_agegrp = c(0.1828,0.6404,0.1768);

# Additionally, we have the following statistics for share without and with college education. Permanently across population
# stat_distr_educ(1,1)=0.6970; % No college
# stat_distr_educ(1,2)=0.3030; % College
# ar_p_edugrp = c(0.6970,0.3030);
# generate a fake third educationcategory
ar_p_edugrp <- c(0.6970,0.3029*(6/10),0.3029*(4/10));

# Joint age and edu distributgion
mt_p_age_edu = (ar_p_agegrp) %*% t(ar_p_edugrp)

# B. P(U|AgeGroup) P(U|EduGroup) -----------------
# From data
ar_p_u_m_agegrp = c(0.176083333, 0.083, 0.0659166667);
ar_p_u_m_edugrp = c(0.100098317, 0.0460833333);
ar_b <- c(ar_p_u_m_agegrp, ar_p_u_m_edugrp)

# C. Use ffi_solve_3by3_rectilinear from fs_discrandvar_condi2joint -----------

# Copied from https://fanwangecon.github.io/R4Econ/statistics/discrandvar/htmlpdfr/fs_discrandvar_condi2joint.html#124_Solution_Program_for_3_by_3_Problem_with_Rectlinear_Restrictions

# Solve for unemploymnet probability:
mt_p_u_joint <- ffi_solve_3by3_rectilinear(mt_p_age_edu, ar_b, bl_verbose=TRUE)
mt_p_u_joint

# Three columns, for the three increasing age groups.
# Two rows for the two education groups
# [,1]       [,2]       [,3]
# [1,] 0.1791864 0.08610302 0.06901968
# [2,] 0.1251714 0.03208803 0.01500470
