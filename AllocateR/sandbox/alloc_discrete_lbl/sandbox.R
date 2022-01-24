# A1. Identify shifts in levels and percentages from current algorithm ---------------------
# At how many positions are V and C adjusted based on nonincreasing adjustments.
df_input_il_noninc <- df_input_il %>%
  arrange((D_il)) %>%
  group_by(id_i) %>%
  do(c_alpha_il_noninc_cumusum = ffi_alpha_non_increasing_adj_cumusum(.$c_alpha_il)) %>%
  unnest(c(c_alpha_il_noninc_cumusum)) %>%
  group_by(id_i) %>%
  mutate(D_il = row_number()) %>%
  left_join(df_input_il, by=(c('id_i'='id_i', 'D_il'='D_il'))) %>%
  mutate(c_alpha_il_gap = c_alpha_il - c_alpha_il_noninc_cumusum) %>%
  mutate(c_alpha_il_ratio = (c_alpha_il - c_alpha_il_noninc_cumusum)/c_alpha_il)

## A2. Add flatten version of array
df_input_il_noninc <- df_input_il_noninc %>%
  arrange((D_il)) %>%
  group_by(id_i) %>%
  do(c_alpha_il_noninc_flatten = ffi_alpha_non_increasing_adj_flatten(.$c_alpha_il)) %>%
  unnest(c(c_alpha_il_noninc_flatten)) %>%
  group_by(id_i) %>%
  mutate(D_il = row_number()) %>%
  left_join(df_input_il_noninc, by=(c('id_i'='id_i', 'D_il'='D_il'))) %>%
  ungroup()

# it_total_rows <- dim(df_input_il_noninc)[1]
# ar_shifts_c_alpha <- (df_input_il_noninc$c_alpha_il_noninc - df_input_il_noninc$c_alpha_il)
# view(df_input_il_noninc[(ar_shifts_c_alpha != 0),])

REconTools::ff_summ_percentiles(df_input_il_noninc)

## B. Identify cases of large changes ------------------------------
## Examine the largeest cases of modification, is the correction correct?
# 27 cases with more than 10 percent changes
ar_id_large_changes <- unique(df_input_il_noninc %>% filter(c_alpha_il_ratio > 0.025) %>% pull(id_i))
print(ar_id_large_changes)

## C. Graph Some of the largest Changes Cases
dev.off()
for (it_id in ar_id_large_changes) {

	df_id_i_1 <- df_input_il_noninc %>% filter(id_i == it_id)
	ar_c_alpha_il <- df_id_i_1$c_alpha_il
	c_alpha_il_noninc_flatten <- df_id_i_1$c_alpha_il_noninc_flatten
	c_alpha_il_noninc_cumusum <- df_id_i_1$c_alpha_il_noninc_cumusum

	mt_c_alpha_il <- cbind(ar_c_alpha_il, c_alpha_il_noninc_flatten, c_alpha_il_noninc_cumusum)
	# show
	par(new=F)
	matplot((mt_c_alpha_il), col = c("black","red","blue"), type='l')

	grid()
}
