# ----------------------------------------------------
# ----------------------------------------------------
# Plot the results for 3 kernels
# ----------------------------------------------------
# ----------------------------------------------------

library(tidyverse)
library(gridExtra)

results_path = "results_csv/nexp_3k.csv"

df = read_csv(results_path) 

# --- Load data

df = df %>%
  mutate(z_score = (C_XY -E_C_XY)/sqrt(Var_C_XY),
         zp_score = (C_XYrZ -E_C_XYrZ)/sqrt(Var_C_XYrZ),
         zpn_score = (nC_XYrZ - E_nC_XYrZ)/sqrt(Var_nC_XYrZ))
  
  
# --- Graph H0

df_H0 = df %>%
  filter(sd_A==0)

z_means = df_H0 %>%
  group_by(sd_B) %>%
  summarise(mean_z = mean(z_score),
            mean_zp = mean(zp_score), 
            mean_zpn = mean(zpn_score)) %>%
  ungroup()

g1 = df_H0 %>%
  ggplot(aes(x=z_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_z), color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=0, color="black", linetype="dotted", size=1) +
  labs(x = "z-score of Dissimilarity Covariance", 
       y = "Density") +
  facet_wrap(~sd_B, labeller = labeller(
    sd_B = function(x) paste0("sd_B", "=", x),
  ))

g2 = df_H0 %>%
  ggplot(aes(x=zp_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_zp), color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=0, color="black", linetype="dotted", size=1) +
  labs(x = "z-score of Partial Dissimilarity Covariance", 
       y = "Density") +
  facet_wrap(~sd_B, labeller = labeller(
    sd_B = function(x) paste0("sd_B", "=", x),
  ))

g3 = df_H0 %>%
  ggplot(aes(x=zpn_score)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_vline(data=z_means, 
             aes(xintercept=mean_zpn), color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=0, color="black", linetype="dotted", size=1) +
  labs(x = "z-score of Partial (naive) Dissimilarity Covariance", 
       y = "Density") +
  facet_wrap(~sd_B, labeller = labeller(
    sd_B = function(x) paste0("sd_B", "=", x),
  ))

# Plot the 3 graphs in one figure
g_all = grid.arrange(g1, g2, g3, nrow=3)
ggsave("results_plot/nexp_H0_3k.png", g_all, width=8, height=8)

# --- Graph H1

df_H1 = df %>%
  group_by(sd_B, sd_A) %>%
  summarise(mean_z = mean(z_score),
            q5 = quantile(z_score, 0.05),
            q95 = quantile(z_score, 0.95),
            mean_zp = mean(zp_score),
            q5_p = quantile(zp_score, 0.05),
            q95_p = quantile(zp_score, 0.95),
            mean_zpn = mean(zpn_score),
            q5_zpn = quantile(zpn_score, 0.05),
            q95_zpn = quantile(zpn_score, 0.95)) %>%
  ungroup()

df_H1 %>%
  ggplot() +
  geom_line(aes(x=sd_A, y=q5), color="Dissimilarity Covariance", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=q95), color="Dissimilarity Covariance", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=mean_z), color="Dissimilarity Covariance", size=1) +
  geom_line(aes(x=sd_A, y=q5_p), color="red", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=q95_p), color="red", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=mean_zp), color="red", size=1) +
  geom_line(aes(x=sd_A, y=q5_pn), color="green", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=q95_pn), color="green", alpha=0.4, lty=3, size=1) +
  geom_line(aes(x=sd_A, y=mean_zpn), color="green", size=1) +
  labs(x = "Standard deviation of the regression matrix",
       y = "Mean z-score of Dissimilarity Covariance") +
  scale_color_manual(name="Method", 
                     values=c("blue", "red", "green"), 
                     labels=c("Dissimilarity Covariance", 
                              "Partial Dissimilarity Covariance", 
                              "Partial (naive) Dissimilarity Covariance")) +
  facet_wrap(~sd_B, labeller = labeller(
    sd_B = function(x) paste0("sd_B", "=", x)
  ), scales = "free")
