# Feasible Allocation Graph, Compare C vs V Allocaitons.
# Feasible allocation doubling adults and Kids First Round
# This is to generate the figure stored with path below:
# "Paper\Appendix_graphs_and_tables\Alternative_Planner_Objective_Feasible_c_allocation_round1_18-64_year-olds_by_income_C_vs_V_allocation.png"

# Libraries
library(tidyverse)
library(REconTools)

# File Names, Paths Etc -------
ls_output <- fs_opti_support('b1_manna')
srt_csv_path_root <- ls_output$srt_csv_path_root
bl_save_img <- ls_output$bl_save_img
srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph
# bl_save_img <- TRUE

# G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
srt_main_folder <- 'b1_a64_dblbth_18t64'
snm_sub_folder <- 'df_alloc_all_feasible_dblbth.csv'
st_img_name_full <- 'Alternative_Planner_Objective_Feasible_c_allocation_round1_18-64_year-olds_by_income_C_vs_V_allocation.png'

# Load in CSV
spn_image_csv_path <- file.path(srt_csv_path_root, srt_main_folder, snm_sub_folder)
df_alloc_cur <- as_tibble(read.csv(spn_image_csv_path)) %>% select(-X)

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
  mutate(checks_c_optimal = optimal_c_1st,
         checks_v_optimal = optimal_v_1st) %>%
  select(-actual, -optimal_c_1st, -optimal_c_2nd, -optimal_v_1st, -optimal_v_2nd) %>%
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
ar_st_age_group_leg_labels <- c("C\u2013Optimal", "V\u2013Optimal")
# color #F8766D is default red
# color #F8766D is default bluish color
ar_st_colours <- c("#33cc33", "#1e4acf")

plt_cur <- plt_cur +
  theme(
    text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.14, 0.9),
    legend.background = element_rect(fill = "white", colour = "black", linetype='solid')) +
  scale_colour_manual(values = ar_st_colours, labels=ar_st_age_group_leg_labels) +
  scale_shape_manual(values = c(17, 15), labels=ar_st_age_group_leg_labels)

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
