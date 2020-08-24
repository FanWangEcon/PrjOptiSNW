# 20K Figure

# Libraries
library(tidyverse)
library(REconTools)

# File Names, Paths Etc -------
ls_output <- fs_opti_support()
st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
srt_simu_path <- ls_output$srt_simu_path
bl_save_img <- ls_output$bl_save_img
spt_img_save <- ls_output$spt_img_save
srt_csv_path_root <- ls_output$srt_csv_path_root
srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph


# G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
srt_main_folder <- 'b1_a64_20k_18t64'
snm_sub_folder <- 'df_alloc_all_feasible_20k.csv'

# Load in CSV
spn_image_csv_path <- file.path(srt_csv_path_root, srt_main_folder, snm_sub_folder)
df_alloc_all_optimalg4_dblbth <- as_tibble(read.csv(spn_image_csv_path)) %>% select(-X)

# Marital and Kids Level Labeling
marry_levels <- c(Single = "0", Married = "1")
kids_levels <- c("no children" = "0", "one child" = "1",
                 "two children" = "2", "three children" = "3",
                 "four children" = "4")

## Select, as factor, and recode
df_alloc_all_optimalg4_dblbth <- df_alloc_all_optimalg4_dblbth %>%
  mutate(kids = as.factor(kids),
         marital = as.factor(marital)) %>%
  mutate(kids = fct_recode(kids, !!!kids_levels),
         marital = fct_recode(marital, !!!marry_levels))


# Dataframe resetting before graphing
df_alloc_all_optimalg4_dblbth_long <- df_alloc_all_optimalg4_dblbth %>%
  filter(rho_val == 1 & mass >= 1e-5) %>% rowwise() %>%
  mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
         y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
  ungroup() %>%
  mutate(ymin_group = as.factor(ymin_group),
         kids = as.factor(kids),
         marital = as.factor(marital)) %>%
  mutate(checks_actual = actual,
         checks_optimal = optimal_c_1st) %>%
  select(-actual, -optimal_c_1st, -optimal_c_1st, -optimal_v_1st, -optimal_v_1st) %>%
  pivot_longer(cols = starts_with('checks'),
               names_to = c('allocatetype'),
               names_pattern = paste0("checks_(.*)"),
               values_to = "checks")

# Basic Graph Structure
plt_cur <- df_alloc_all_optimalg4_dblbth_long %>%
  # filter(allocatetype == 'optimal') %>%
  mutate(y_group_min = as.numeric(y_group_min),
         checks = as.numeric(checks)) %>%
  ggplot(aes(x=y_group_min, y=checks*100,
             colour=allocatetype,
             shape=allocatetype)) +
  facet_wrap(~ marital + kids, ncol=5, labeller = label_wrap_gen(multi_line=FALSE))

plt_cur <- plt_cur + geom_line() + geom_point(size=3)

# X and Y labels and Legends etc.
x.labels <- c('0', '50K', '100K', '150K', '200K')
x.breaks <- c(0,
              50000/58056,
              100000/58056,
              150000/58056,
              200000/58056)

plt_cur <- plt_cur +
  theme(
    text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.14, 0.9),
    legend.background = element_rect(fill = "white", colour = "black", linetype='solid')) +
  scale_x_continuous(labels = x.labels, breaks = x.breaks,
                     limits = c(0, 175000/58056))

stg_x_title_hhinc <- 'Household income (thousands of 2012 USD)'
stg_y_title_checks <- 'Amount of welfare checks (2012 USD)'

plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)


# PNG and print
st_img_name_full <- 'Optimalg4_20k_C_allocation_by_income_kids_marital_status_round1.png'
png(paste0(srt_paper_appendix_textgraph, st_img_name_full),
    width = 270,
    height = 216, units='mm',
    res = 300, pointsize=7)
print(plt_cur)
dev.off()

