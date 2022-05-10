library(tidyverse)
library(mvnfast)

### Not independent case
func2 = function(num.samples = 100, p = 0.4, rho = 0, alt.hyp = 3) {
  # Sigma has correlation structure rho and 1 on diagonal
  # num.samples = 100
  # p = 0.4
  # rho = 0.5
  # alt.hyp = 3
  Sigma = matrix(nrow = num.samples, ncol = num.samples, data = rho)
  diag(Sigma) = 1
  # For known p, simulate the true hypothesis state
  theta = rbinom(num.samples, 1, p)
  # Generate positively correlated data
  # MVN = MASS::mvrnorm(n = 1, 
  #                     mu =  ifelse(theta, alt.hyp, 0), 
  #                     Sigma = Sigma)
  MVN = matrix(numeric(num.samples), nrow = 1, ncol = num.samples)
  mvnfast::rmvn(n = 1, mu = ifelse(theta, alt.hyp, 0), sigma = Sigma, A = MVN)
  # Get p-values for BH and adjust them
  MVN = t(MVN)
  MVN_p = pnorm(MVN, mean = 0, sd = 1, lower.tail = F)
  BH_p =  p.adjust(MVN_p, method = "BH")
  adaptBH_p = p.adjust(MVN_p*((num.samples - sum(MVN_p <= 0.5) + 1)/(num.samples*(1-0.5))), method = "BH")
  bonf_p = p.adjust(MVN_p, method = "bonf")
  # Copy/paste code from before to get Lfdr
  df = tibble(theta, 
              MVN_p,
              BH_p,
              adaptBH_p,
              bonf_p,
              f = (1-p)*dnorm(MVN, 0, 1) + p*dnorm(MVN, alt.hyp, 1),
              Lfdr = (1-p)*dnorm(MVN, 0, 1)/f) %>%
    arrange(Lfdr) %>%
    mutate(q = cumsum(Lfdr)/(1:length(Lfdr)),
           `p<a` = MVN_p < 0.05,
           `q<a` = q < 0.05,
           `BH<a` = BH_p < 0.05,
           `adaptBH<a` = adaptBH_p < 0.05,
           `bonf<a` = bonf_p < 0.05,
           # Calculate errors for BH and Lfdr methods
           Lfdr_Type1error = `q<a`&!theta,
           Lfdr_Type2error = !`q<a`&theta,
           BH_Type1error = `BH<a` & !theta,
           BH_Type2error = theta & !`BH<a`,
           adaptBH_Type1error = `adaptBH<a` & !theta,
           adaptBH_Type2error = theta & !`adaptBH<a`,
           Bonf_Type1error = `bonf<a`&!theta,
           Bonf_Type2error = !`bonf<a`&theta,
           None_Type1error = `p<a`&!theta,
           None_Type2error = !`p<a`&theta) %>% 
    summarise(theta = sum(theta),
              Lfdr_disc = sum(`q<a`),
              BH_disc = sum(`BH<a`),
              adaptBH_disc = sum(`adaptBH<a`),
              Bonf_disc = sum(`bonf<a`),
              None_disc = sum(`p<a`),
              Lfdr_T1 = sum(Lfdr_Type1error),
              Lfdr_T2 = sum(Lfdr_Type2error),
              BH_T1 = sum(BH_Type1error),
              BH_T2 = sum(BH_Type2error),
              adaptBH_T1 = sum(adaptBH_Type1error),
              adaptBH_T2 = sum(adaptBH_Type2error),
              Bonf_T1 = sum(Bonf_Type1error),
              Bonf_T2 = sum(Bonf_Type2error),
              None_T1 = sum(None_Type1error),
              None_T2 = sum(None_Type2error),
              p, rho, num.samples, alt.hyp)
  return(df)
}

params = expand.grid(rho = rep(seq(0, 0.8, by = 0.2), times = 10),
                     num.samples = c(30, 100, seq(500, 4000, by = 500)),
                     p = c(seq(0.1, 0.9, by = 0.1)), 
                     alt_hyp = 2:5)  

parallel::mcmapply(FUN = func2, 
                   mc.cores = parallel::detectCores(), 
                   p = params$p, 
                   num.samples = params$num.samples, 
                   rho = params$rho, 
                   alt.hyp = params$alt_hyp) -> c
sim_data = t(c) %>% 
  as_tibble(column_name = rownames(c)) %>% 
  mutate_all(as.double) %>%
  mutate(BH_FDR = BH_T1/pmax(BH_disc, 1),
         BH_FNR = ifelse(theta == 0, 0, BH_T2/theta),
         Lfdr_FDR = Lfdr_T1/pmax(Lfdr_disc, 1),
         Lfdr_FNR = ifelse(theta == 0, 0, Lfdr_T2/theta),
         adaptBH_FNR = ifelse(theta == 0, 0, adaptBH_T2/theta),
         adaptBH_FDR = adaptBH_T1/pmax(adaptBH_disc, 1), 
         Bonf_FNR = ifelse(theta == 0, 0, Bonf_T2/theta),
         Bonf_FDR = Bonf_T1/pmax(Bonf_disc, 1),
         None_FNR = ifelse(theta == 0, 0, None_T2/theta),
         None_FDR = None_T1/pmax(None_disc, 1),
         .keep = "unused") %>%
  pivot_longer(cols = contains("_"),
               names_to = c("Method", "Type"),
               names_sep = "_", 
               values_to = "Prob") # %>% 

write_csv(sim_data, "sim_data_final3.csv")

sim_data = read_csv("sim_data_final.csv") %>% 
  bind_rows(read_csv("sim_data_final2.csv"), 
            read_csv("sim_data_final3.csv"))
###################
# FDR, independent plots
fdr_a1 = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          num.samples == 1000,
                                          rho == 0,
                                          alt.hyp == 3),
               aes(x = p, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  # geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Proportion of alt hypotheses" ~(pi)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw()

fdr_b1 = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          p == 0.2,
                                          num.samples != 30,
                                          rho == 0,
                                          alt.hyp == 3),
               aes(x = num.samples, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Number of tests" ~(m)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw()


fdr_c1 = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          p == 0.2,
                                          num.samples == 1000,
                                          #rho == 0,
                                          alt.hyp == 3),
               aes(x = rho, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Amount of correlation" ~(rho)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw()

fdr_d1 = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          p == 0.2,
                                          num.samples == 1000,
                                          rho == 0,
                                          #alt.hyp == 3
                                          ),
aes(x = alt.hyp, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Alternative hypothesis mean" ~(mu[1])),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw()

####################
# FNR, independent plots
fnr_a1 = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          num.samples == 1000,
                                          rho == 0,
                                          alt.hyp == 3),
               aes(x = p, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  # geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Proportion of alt hypotheses" ~(pi)),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme_bw()

fnr_b1 = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          p == 0.2,
                                          num.samples != 30,
                                          rho == 0,
                                          alt.hyp == 3),
               aes(x = num.samples, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Number of tests" ~(m)),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()


fnr_c1 = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          p == 0.2,
                                          num.samples == 1000,
                                          #rho == 0,
                                          alt.hyp == 3),
               aes(x = rho, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Amount of correlation" ~(rho)),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

fnr_d1 = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          p == 0.2,
                                          num.samples == 1000,
                                          rho == 0,
                                          #alt.hyp == 3
                                          ),
aes(x = alt.hyp, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Alternative hypothesis mean" ~(mu[1])),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()


###################
# FDR, not independent plots
fdr_a = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          num.samples == 1000,
                                          rho == 0.4,
                                          alt.hyp == 3),
       aes(x = p, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  # geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Proportion of alt hypotheses" ~(pi)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw()

fdr_b = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples != 30,
                                          rho == 0.4,
                                          alt.hyp == 3),
               aes(x = num.samples, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Number of tests" ~(m)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw()


fdr_c = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples == 1000,
                                          #rho == 0,
                                          alt.hyp == 3),
               aes(x = rho, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Amount of correlation" ~(rho)),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw()

fdr_d = ggplot(data = sim_data %>% filter(Type == "FDR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples == 1000,
                                          rho == 0.4,
                                          #alt.hyp == 3
                                          ),
               aes(x = alt.hyp, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Alternative hypothesis mean" ~(mu[1])),
       y = "FDR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw()

####################
# FNR, not independent plots
fnr_a = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          num.samples == 1000,
                                          rho == 0.4,
                                          alt.hyp == 3),
               aes(x = p, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  # geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Proportion of alt hypotheses" ~(pi)),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(0, 1)) +
  theme_bw()

fnr_b = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples != 30,
                                          rho == 0.4,
                                          alt.hyp == 3),
               aes(x = num.samples, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Number of tests" ~(m)), 
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_bw()


fnr_c = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples == 1000,
                                          #rho == 0,
                                          alt.hyp == 3),
               aes(x = rho, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Amount of correlation" ~(rho)),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_bw()

fnr_d = ggplot(data = sim_data %>% filter(Type == "FNR",
                                          Method %in% c("adaptBH", "Lfdr", "BH"),
                                          p == 0.2,
                                          num.samples == 1000,
                                          rho == 0.4,
                                          #alt.hyp == 3
                                          ),
aes(x = alt.hyp, y = Prob, color = Method, shape = Method)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #geom_point(size = 0.3, alpha = 0.3) +
  # geom_hline(yintercept = 0.05) +
  # facet_grid(rho ~ num.samples) +
  labs(x = expression("Alternative hypothesis mean" ~(mu[1])),
       y = "FNR") +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Number of tests",
  #                                        breaks = NULL, labels = NULL)) +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Amount of correlation",
  #                                        breaks = NULL, labels = NULL)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

###################
# cowplots
A = cowplot::plot_grid(fdr_a1, fdr_b1, fdr_c1, fdr_d1, labels = "A")
B = cowplot::plot_grid(fnr_a1, fnr_b1, fnr_c1, fnr_d1, labels = "B")
C = cowplot::plot_grid(fdr_a, fdr_b, fdr_c, fdr_d, labels = "C")
D = cowplot::plot_grid(fnr_a, fnr_b, fnr_c, fnr_d, labels = "D")
cowplot::save_plot("IndepFDR_plot.png", A)
cowplot::save_plot("IndepFNR_plot.png", B)
cowplot::save_plot("FDR_plot.png", C)
cowplot::save_plot("FNR_plot.png", D)
