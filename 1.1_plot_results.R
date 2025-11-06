# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 2 variables
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)

# --- Parameters

r_w_v = seq(0, 1, 0.1)

# --- Load data

df = data.frame()

for(r_w in r_w_v){
  results_path = paste0("results_csv/results_rw_", r_w, ".csv")
  temp_df = read_csv(results_path) 
  df = rbind(df, temp_df)
}

# --- Summarise data

summ_df = df %>% 
  group_by(r1, r12) %>%
  summarise(n = n(),
            mean_C = mean(C_XY),
            sd_C = sd(C_XY),
            mean_C_th = mean(E_C_XY),
            sd_C_th = sd(Var_C_XY))

# --- Plot the two curves with confidence intervals

ggplot(summ_df) +
  geom_line(aes(x = r12, y = mean_C, color = bquote("Empirical Cov")), size=1) +
  geom_ribbon(aes(x = r12, 
                  ymin = mean_C - 1.96*sd_C/sqrt(n), 
                  ymax = mean_C + 1.96*sd_C/sqrt(n)), 
              alpha = 0.2, fill="blue") +
  geom_line(aes(x = r12, y = mean_C_th, color = "Theoretical Cov"), size=1) +
  geom_ribbon(aes(x = r12, 
                  ymin = mean_C_th - 1.96*sd_C_th, 
                  ymax = mean_C_th + 1.96*sd_C_th), 
              alpha = 0.2, fill="red") +
  labs(x = bquote("Correlation between datasets ("~r[b]~")"), 
       y = "Dissimilariy Covariance", 
       title = "Dissimilarity Covariance between 2 Kernels") +
  scale_color_manual(values = c("Empirical Cov" = "blue", "Theoretical Cov" = "red")) +
  facet_wrap(~r1, labeller = labeller(r1 = function(x) paste0("r_w", "=", x) )) +
  theme(legend.title=element_blank())
ggsave("results_plot/dissimilarity_covariance_2Kernels.png", width=8, height=6)


            