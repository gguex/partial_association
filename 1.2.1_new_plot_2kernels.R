# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 2 kernels
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)

results_path = "results_csv/nexp_2k.csv"

df = read_csv(results_path) 

df = df %>%
  mutate(z_score = (C_XY -E_C_XY)/sqrt(Var_C_XY))
  
# --- Graph 1.1

df_H0 = df %>%
  filter(sd_A==0)

mean_H0 = mean(df_H0$z_score)

df_H0 %>%
  ggplot(aes(x=z_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(aes(xintercept=mean_H0), color="red", linetype="dashed", size=1) +
  labs(x = "z-score of Dissimilarity Covariance", 
       y = "Density", 
       title = "Histogram of z-scores when sd_A = 0") +
ggsave("results_plot/nexp_H0_2k.png", width=8, height=6)

# --- Graph 1.2

df_H1 = df %>%
  group_by(sd_A) %>%
  summarise(mean_z = mean(z_score),
            q5 = quantile(z_score, 0.05),
            q95 = quantile(z_score, 0.95)) %>%
  ungroup()

df_H1 %>%
  ggplot(aes(x=sd_A, y=mean_z)) +
  geom_line(color="blue", size=1) +
  geom_ribbon(aes(ymin=q5, ymax=q95), alpha=0.8, fill="lightblue") +
  labs(x = bquote("Standard deviation of the regression matrix"),
       y = "Mean z-score of Dissimilarity Covariance", 
       title = "Mean z-score of Dissimilarity Covariance vs sd_A")
ggsave("results_plot/nexp_H1_2k.png", width=8, height=6)
