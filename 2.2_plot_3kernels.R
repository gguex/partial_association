# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 3 kernels
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)

# --- Parameters

r_w_vec = c(0, 0.5, 0.9)

# --- Load data

df = data.frame()

for(r_w in r_w_vec){
  results_path = paste0("results_csv/new_3k_rw_", r_w, ".csv")
  temp_df = read_csv(results_path) 
  df = rbind(df, temp_df)
}

df = df %>%
  mutate(z_score = (C_XYrZ -E_C_XYrZ)/sqrt(Var_C_XYrZ),
         z_score_n = (nC_XYrZ - E_nC_XYrZ)/sqrt(Var_nC_XYrZ))
  
  
# --- Graph 1.1

df_g1 = df %>%
  filter(r12==0)

z_means = df_g1 %>%
  group_by(r1, r13) %>%
  summarise(mean_z = mean(z_score), 
            mean_z_n = mean(z_score_n),
            semidef_pos = max(semidef_pos)) %>%
  ungroup()

df_g1 %>%
  ggplot(aes(x=z_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_z), color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=0, color="black", linetype="dotted", size=1) +
  geom_text(aes(x=mean_z_n, y=0.3, label=ifelse(semidef_pos==0, "Not semi-def pos", "")), 
            data=z_means, color="red", vjust=-1) +
  labs(x = "z-score of Dissimilarity Covariance", 
       y = "Density", 
       title = "Histogram of z-scores when r_b = 0") +
  facet_wrap(~r13+r1, labeller = labeller(
    r1 = function(x) paste0("r_w", "=", x),
    r13 = function(x) paste0("r_c", "=", x)
  ))
ggsave("results_plot/z_new_histogram_3kernels.png", width=8, height=6)

df_g1 %>%
  ggplot(aes(x=z_score_n)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_z_n), color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=0, color="black", linetype="dotted", size=1) +
  geom_text(aes(x=mean_z_n, y=0.3, label=ifelse(semidef_pos==0, "Not semi-def pos", "")), 
            data=z_means, color="red", vjust=-1) +
  labs(x = "z-score of Dissimilarity Covariance", 
       y = "Density", 
       title = "Histogram of 'naive' z-scores when r_b = 0") +
  facet_wrap(~r13+r1, labeller = labeller(
    r1 = function(x) paste0("r_w", "=", x),
    r13 = function(x) paste0("r_c", "=", x)
  ))
ggsave("results_plot/z_new_histogram_naive_3kernels.png", width=8, height=6)

# --- Graph 2.2

df_g2 = df %>%
  group_by(r1, r12, r13) %>%
  filter(semidef_pos==1) %>%
  summarise(mean_z = mean(z_score),
            q5 = quantile(z_score, 0.05),
            q95 = quantile(z_score, 0.95),
            mean_z_n = mean(z_score_n),
            q5_n = quantile(z_score_n, 0.05),
            q95_n = quantile(z_score_n, 0.95)) %>%
  ungroup()

df_g2 %>%
  ggplot() +
  geom_line(aes(x=r12, y=mean_z), color="blue", size=1) +
  geom_ribbon(aes(x=r12, ymin=q5, ymax=q95), alpha=0.8, fill="lightblue") +
  labs(x = bquote("Correlation between datasets ("~r[b]~")"),
       y = "Mean z-score of Partial Dissimilarity Covariance", 
       title = "Mean z-score of Partial Dissimilarity Covariance vs r_b") +
  facet_wrap(~r13+r1, labeller = labeller(
    r1 = function(x) paste0("r_w", "=", x),
    r13 = function(x) paste0("r_c", "=", x)
  ), scales = "free")
ggsave("results_plot/z_new_evo_3kernels.png", width=8, height=6)

# --- Graph 2.3

# df_g3 = df %>%
#   filter(r1==0.5) %>%
#   group_by(r12, r13) %>%
#   summarise(mean_C_XY_Z = mean(C_XY_Z),
#             mean_C_XYrZ = mean(C_XYrZ)) %>%
#   ungroup()

# df_g3 %>%
#   # stack plot
#   pivot_longer(cols=c("mean_C_XY_Z", "mean_C_XYrZ"), 
#                names_to="type", values_to="mean_value") %>%
#   ggplot(aes(x=r12, y=mean_value, fill=type)) +
#   geom_area(position="stack", alpha=0.8) +
#   labs(x = bquote("Correlation between datasets ("~r[b]~")"),
#        y = "Mean Dissimilarity Covariance",
#        title = "Mean Dissimilarity Covariance vs r_b when r_w = 0.5") +
#   scale_fill_discrete(name="Type", labels=c("C_XY|Z", "
# C_XYrZ")) +
#   facet_wrap(~r13, labeller = labeller(
#     r13 = function(x) paste0("r_c", "=", x),
#     scales = "free"
#   ))
# ggsave("results_plot/cov_addition_3kernels.png", width=8, height=6)



