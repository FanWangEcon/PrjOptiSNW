# Compare Three Allocations, b0xi0p25manna vs b1manna vs b1xi0tax
# Feasible allocation doubling adults and Kids First Round
# This is to generate the figure stored with path below:
# "Paper\Appendix_graphs_and_tables\Robustness_Tax_model_Actual_vs_feasible_c_allocation_round1_18-64_year-olds_by_income.png"

# Libraries
library(tidyverse)
library(REconTools)
library(scales)

# File Names, Paths Etc -------
ls_output <- fs_opti_support('b1_manna')
srt_results_root <- ls_output$srt_results_root
bl_save_img <- ls_output$bl_save_img
bl_save_img <- FALSE
srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph
# bl_save_img <- TRUE

# G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
srt_simu_root_b0_xi0p25_manna <- '2020-08-23-b0_xi0p25_manna/csv'
srt_simu_root_b1_manna <- '2020-08-23-b1_manna/csv'
srt_simu_root_b1_xi0_tax <- '2020-08-23-b1_xi0_tax/csv'
srt_main_folder_b0 <- 'b0_a64_dblbth_18t64'
srt_main_folder_b1 <- 'b1_a64_dblbth_18t64'
snm_sub_folder <- 'df_alloc_all_feasible_dblbth.csv'
st_img_name_full <- 'Robustness_Tax_model_Actual_vs_feasible_c_allocation_round1_18-64_year-olds_by_income.png'

# Load in CSV
spn_image_csv_path_b1_manna <- file.path(srt_results_root, srt_simu_root_b1_manna, srt_main_folder_b1, snm_sub_folder)
spn_image_csv_path_b0_xi0p25_manna <- file.path(srt_results_root, srt_simu_root_b0_xi0p25_manna, srt_main_folder_b0, snm_sub_folder)
spn_image_csv_path_b1_xi0_tax <- file.path(srt_results_root, srt_simu_root_b1_xi0_tax, srt_main_folder_b1, snm_sub_folder)

df_alloc_cur_b1_manna <- as_tibble(read.csv(spn_image_csv_path_b1_manna)) %>%
  mutate(optimal_c_1st_b1_manna = optimal_c_1st) %>%
    select(-X, -optimal_c_1st) %>% ungroup()
df_alloc_cur_b0_xi0p25_manna <- as_tibble(read.csv(spn_image_csv_path_b0_xi0p25_manna)) %>%
  select(id_i, optimal_c_1st) %>%
  mutate(optimal_c_1st_b0_xi0p25_manna = optimal_c_1st) %>%
  select(-optimal_c_1st) %>% ungroup()
df_alloc_cur_b0_xi0_tax <- as_tibble(read.csv(spn_image_csv_path_b1_xi0_tax)) %>%
  select(id_i, optimal_c_1st) %>%
  mutate(optimal_c_1st_b0_xi0_tax = optimal_c_1st) %>%
  select(-optimal_c_1st) %>% ungroup()

# Merging FIles Different Allocations are Additional Columns ------
df_alloc_cur <- df_alloc_cur_b1_manna %>%
  left_join(df_alloc_cur_b0_xi0p25_manna, by='id_i') %>%
  left_join(df_alloc_cur_b0_xi0_tax, by='id_i')

# Modifying dataframe ------
# Marital and Kids Level Labeling
marry_levels <- c(Single = "0", Married = "1")
kids_levels <- c("no children" = "0", "one child" = "1",
                 "two children" = "2", "three children" = "3",
                 "four children" = "4")

## Select, as factor, and recode
df_alloc_cur <- df_alloc_cur %>%
  mutate(kids = as.factor(kids),
         marital = as.factor(marital)) %>%
  mutate(kids = fct_recode(kids, !!!kids_levels),
         marital = fct_recode(marital, !!!marry_levels))

# Get value for minimum and maximum income levels for each bin
df_alloc_cur <- df_alloc_cur %>%
  filter(rho_val == 1) %>% rowwise() %>%
  mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
         y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
  ungroup() %>%
  mutate(ymin_group = as.factor(ymin_group),
         kids = as.factor(kids),
         marital = as.factor(marital))

# Get Optimal Allocation Result fro Graph
# 1st period optimal c allocation
df_alloc_cur_long <- df_alloc_cur %>%
  mutate(checks_c_optimal_b1_manna = optimal_c_1st_b1_manna,
         checks_c_optimal_b0_xi0p25_manna = optimal_c_1st_b0_xi0p25_manna,
         checks_c_optimal_b0_xi0_tax = optimal_c_1st_b0_xi0_tax) %>%
  select(-actual, -optimal_c_2nd, -optimal_v_1st, -optimal_v_2nd,
         -optimal_c_1st_b1_manna, -optimal_c_1st_b0_xi0p25_manna, -optimal_c_1st_b0_xi0_tax) %>%
  pivot_longer(cols = starts_with('checks'),
               names_to = c('allocatetype'),
               names_pattern = paste0("checks_(.*)"),
               values_to = "checks")

# Filtering data for allocation
df_alloc_cur_long <- df_alloc_cur_long %>%
  mutate(y_group_min = as.numeric(y_group_min),
         checks = as.numeric(checks))

# GGPLOT grahing ------
# GGPLOT AES and FACET_WRAP
plt_cur <- df_alloc_cur_long %>%
  ggplot(aes(x=y_group_min, y=checks*100,
             colour=allocatetype,
             shape=allocatetype)) +
  facet_wrap(~ marital + kids,
             ncol=5,
             scales = "free_x",
             labeller = label_wrap_gen(multi_line=FALSE))
plt_cur <- plt_cur + geom_line() + geom_point(size=3)

# GGPLOT legendings and coloring -------
# X and Y labels and Legends etc.
x.labels <- c('0', '50K', '100K', '150K')
x.breaks <- c(0,
              50000/58056,
              100000/58056,
              150000/58056)

# x-axis labeling
plt_cur <- plt_cur +
  scale_x_continuous(labels = x.labels, breaks = x.breaks,
                     limits = c(0, 175000/58056))
# limits = c(0, 240000/58056)

# Legend Labeling
# \u2013 allows for en-dash
ar_st_age_group_leg_labels <- c("No tax increase (b=1)", "No tax increase (b=0,\u03be=0.25)", "Tax increase (b=1,\u03be=0)")
# color #F8766D is default red
# color #F8766D is default bluish color
ar_st_colours <- c("#33cc33", hue_pal()(8)[c(2,8)])

plt_cur <- plt_cur +
  theme(
    text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.14, 0.9),
    legend.background = element_rect(fill = "white", colour = "black", linetype='solid')) +
  scale_colour_manual(values=ar_st_colours, labels=ar_st_age_group_leg_labels) +
  scale_shape_manual(values=c(17,15,18), labels=ar_st_age_group_leg_labels)

# X and Y Titling
stg_x_title_hhinc <- 'Household income (thousands of 2012 USD)'
stg_y_title_checks <- 'Stimulus check amount (2012 USD)'
plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)

# Save Image Outputs -----
# PNG and print
if (bl_save_img) {
  png(paste0(srt_paper_appendix_textgraph, st_img_name_full),
      width = 270,
      height = 216, units='mm',
      res = 300, pointsize=7)
}
print(plt_cur)
if (bl_save_img) {
  dev.off()
}
