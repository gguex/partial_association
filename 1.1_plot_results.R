# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 2 variables
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)

# --- Parameters


results_path = "results_csv/results_rw_0.5.csv"

# --- Plot

df = read_csv(results_path) 

summ_df = df %>% 
  group_by(r12) %>%
  summarise(n = n(),
            mean_C = mean(C_XY),
            sd_C = sd(C_XY),
            mean_C_th = mean(E_C_XY),
            sd_C_th = sd(Var_C_XY))

# Plot the two curves with confidence intervals
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
       title = bquote("Dissimilarity Covariance (with"~r[w]~"= 0.5)")) +
  scale_color_manual(values = c("Empirical C" = "blue", "Theoretical C" = "red")) +
  theme_minimal() +
  theme(legend.title=element_blank())



            