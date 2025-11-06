# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 2 kernels
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)

# --- Parameters

r_b_vec = seq(0, 1, 0.1)

# --- Load data

df = data.frame()

for(r_b in r_b_vec){
  results_path = paste0("results_csv/res_2kernels_rb_", r_b, ".csv")
  temp_df = read_csv(results_path) 
  df = rbind(df, temp_df)
}

df = df %>%
  mutate(z_score = (C_XY -E_C_XY)/sqrt(Var_C_XY))
  
# --- Graph 1.1

df_g1 = df %>%
  filter(r12==0)

z_means = df_g1 %>%
  group_by(r1) %>%
  summarise(mean_z = mean(z_score)) %>%
  ungroup()

df_g1 %>%
  ggplot(aes(x=z_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_z), color="red", linetype="dashed", size=1) +
  labs(x = "z-score of Dissimilarity Covariance", 
       y = "Density", 
       title = "Histogram of z-scores when r_b = 0") +
  facet_wrap(~r1, labeller = labeller(r1 = function(x) paste0("r_w", "=", x) ))
ggsave("results_plot/zscore_histogram_2kernels.png", width=8, height=6)

# --- Graph 1.2

df_g2 = df %>%
  #filter(semipos_def==1) %>%
  group_by(r1, r12) %>%
  summarise(mean_z = mean(z_score),
            q5 = quantile(z_score, 0.05),
            q95 = quantile(z_score, 0.95)) %>%
  ungroup()

df_g2 %>%
  ggplot(aes(x=r12, y=mean_z)) +
  geom_line(color="blue", size=1) +
  geom_ribbon(aes(ymin=q5, ymax=q95), alpha=0.8, fill="lightblue") +
  labs(x = bquote("Correlation between datasets ("~r[b]~")"),
       y = "Mean z-score of Dissimilarity Covariance", 
       title = "Mean z-score of Dissimilarity Covariance vs r_b") +
  facet_wrap(~r1, labeller = labeller(r1 = function(x) paste0("r_w", "=", x) ))
ggsave("results_plot/zscore_evo_2kernels.png", width=8, height=6)
