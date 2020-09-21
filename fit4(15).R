
rm(list=ls())  # clear memory
library(rstan)
library(tidyverse)
library(bayesplot)
library(deSolve)
library(readxl)

WT_viral_load_data <- read_excel("PR8-viral-load.xls", sheet = 1) %>% rename(Days = Days, ViralLoad = `Viral load`)
KO_viral_load_data <- read_excel("PR8-viral-load.xls", sheet = 2) %>% rename(Days = Days, ViralLoad = `Viral load`)

WT_macrophage_data <- read_excel("Macrophage.xls", sheet = 1) %>% rename(Days = Days, Macrophage = 'Macro level')
KO_macrophage_data <- read_excel("Macrophage.xls", sheet = 2) %>% rename(Days = Days, Macrophage = 'Macro level')


time_data <- c(1,2,3,5,7)
V_data_WT <- WT_viral_load_data$ViralLoad
V_data_KO <- KO_viral_load_data$ViralLoad


time_data_macrophage <- c(1,3,5,7)
Macrophage_data_WT <- WT_macrophage_data$Macrophage
Macrophage_data_KO <- KO_macrophage_data$Macrophage




data_combined_muc1_V_W <-  list(N_T_WT = length(time_data),
                                N_V_WT = length(V_data_WT),
                                time_data_WT = time_data,
                                log_viral_load_data_WT = log(V_data_WT), # above V_WT
                                N_T_KO = length(time_data),
                                N_V_KO = length(V_data_KO),
                                time_data_KO = time_data,
                                log_viral_load_data_KO = log(V_data_KO), # above V_KO
                                N_T_Macrophage_WT = length(time_data_macrophage),
                                N_Macrophage_WT = length(Macrophage_data_WT),
                                time_data_Macrophage_WT = time_data_macrophage,
                                Macrophage_data_WT = log(Macrophage_data_WT), # above Macrophage_WT
                                N_T_Macrophage_KO = length(time_data_macrophage),
                                N_Macrophage_KO = length(Macrophage_data_KO),
                                time_data_Macrophage_KO = time_data_macrophage,
                                Macrophage_data_KO = log(Macrophage_data_KO), # above Macrophage_KO
                                t0 = 0,
                                T0 = 4e+8,
                                I0 = 0,
                                MP0 = 0)

init1 <- list(
  theta = c(0.6, 0.8, 2e-5, 1.3, 5e-3, 1.3, 1e-7, 1e-6, 3e+2, 1, 1e-2, 120, 0.7, 0.8, 5),
  sigma = c(0.5, 0.5))

#init2 <- list(
#  theta = c(0.8, 1, 1.5e-5, 5, 2.2e-3, 3, 1e-5, 4e+3, 4e-3, 0.8, 1, 1e-5, 2, 1.5e-1),
#  sigma = c(0.5,0.5))

options(mc.cores=parallel::detectCores()) # to utilise all cores available in your computer

Estimate_Full4 <- stan("Fit4.stan",
                       data = data_combined_muc1_V_W,
                       pars = c("theta","sigma"),
                       seed = 20200921,  # set random seed for reproducibility
                       iter = 2000,
                       chains = 1,
                       init = list(init1),
                       warmup = 1000,
                       control = list(adapt_delta = 0.99, max_treedepth = 15))


print(Estimate_Full2, pars = c("theta","sigma"))
stan_dens(Estimate_Full2, pars = c("theta"), separate_chains = TRUE,nrow = 6)


# extract posterior samples for selected parameters
posterior_samples_all = rstan::extract(Estimate_Full2, pars = c("theta"), inc_warmup = TRUE, permuted = FALSE)
posterior_samples_merged_after_burnin = rstan::extract(Estimate_Full2, pars = c("theta"))

# show markov chains
color_scheme_set("brewer-Spectral")
mcmc_trace(posterior_samples_all, n_warmup = 1000,
           facet_args = list(nrow = 6, labeller = label_parsed))

# show all marginal posterior distributions
posterior_sample_table = data.frame(epsilon1 = posterior_samples_merged_after_burnin$theta[,1],
                                    k1 = posterior_samples_merged_after_burnin$theta[,2],
                                    beta = posterior_samples_merged_after_burnin$theta[,3],
                                    delta_I = posterior_samples_merged_after_burnin$theta[,4],
                                    p = posterior_samples_merged_after_burnin$theta[,5],
                                    c = posterior_samples_merged_after_burnin$theta[,6],
                                    kappa_M = posterior_samples_merged_after_burnin$theta[,7],
                                    s = posterior_samples_merged_after_burnin$theta[,8],
                                    delta_M = posterior_samples_merged_after_burnin$theta[,9],
                                    epsilon2 = posterior_samples_merged_after_burnin$theta[,10],
                                    k2 = posterior_samples_merged_after_burnin$theta[,11],
                                    q_V = posterior_samples_merged_after_burnin$theta[,12],
                                    eta = posterior_samples_merged_after_burnin$theta[,13],
                                    V0 = posterior_samples_merged_after_burnin$theta[,14])


ggplot(posterior_sample_table, aes(x=epsilon1)) + geom_histogram(breaks=seq(0,1,1/100))+theme_bw()
ggplot(posterior_sample_table, aes(x=epsilon2)) + geom_histogram(breaks=seq(0,1,1/100))+theme_bw()


ggplot(posterior_sample_table, aes(x=k1)) + geom_histogram(breaks=seq(0,2.5,1/100))+theme_bw()
ggplot(posterior_sample_table, aes(x=k2)) + geom_histogram(breaks=seq(0,2.5,1/100))+theme_bw()

ggplot(posterior_sample_table, aes(x=beta)) + geom_histogram(breaks=seq(0,2e-4,2e-4/100))+theme_bw()
ggplot(posterior_sample_table, aes(x=delta_I)) + geom_histogram(breaks=seq(0,2,2/100))+theme_bw()
ggplot(posterior_sample_table, aes(x=q_V)) + geom_histogram(breaks=seq(0,1e-4,1e-4/30))+theme_bw()



# show parameter correlations
pairs(posterior_samples_merged_after_burnin$theta[,10:12], col = c("blue", "red"))#[,1:7], labels = c("epsilon1","k1","beta","delta_I", "p", "c","kappa_M"), col = c( "blue", "red"))


# posterior predictive check
Within_host_model = function(t, y, param){
  dydt1  = -(1 - param[1] * t^10/ (t^10 + param[2]^10)) * param[3] * y[1] * y[3]
  dydt2 = (1 - param[1] * t^10/ (t^10 + param[2]^10)) * param[3] * y[1] * y[3] - param[4] * y[2]
  dydt3 = param[5] * y[2] - param[6] * y[3] - param[7] * y[3] * y[5]
  dydt4 = param[8] - param[9] * y[4] - (1 - param[10] * t^10/(t^10 + param[11]^10)) * param[12] * y[3] * y[4] + param[13] * y[5];
  dydt5 = (1 - param[10] * t^10/(t^10 + param[11]^10)) * param[12] * y[3] * y[4] - param[13] * y[5] - param[9] * y[5];
  
  list(c(dydt1, dydt2, dydt3, dydt4, dydt5))
  
}


t_ppc = seq(0, 10, 0.1)
V_ppc_WT = matrix(, nrow = 1000, ncol = length(t_ppc))
M_ppc_WT = matrix(, nrow = 1000, ncol = length(t_ppc))

V_ppc_KO = matrix(, nrow = 1000, ncol = length(t_ppc))
M_ppc_KO = matrix(, nrow = 1000, ncol = length(t_ppc))

lower_95PI_WT = t_ppc
median_95PI_WT = t_ppc
upper_95PI_WT = t_ppc

lower_95PI_KO = t_ppc
median_95PI_KO = t_ppc
upper_95PI_KO = t_ppc


lower_95PI_WT_M = t_ppc
median_95PI_WT_M = t_ppc
upper_95PI_WT_M = t_ppc

lower_95PI_KO_M = t_ppc
median_95PI_KO_M = t_ppc
upper_95PI_KO_M = t_ppc



for (i in 1:1000){
  y_init = c(data_combined_muc1_V_W$T0, data_combined_muc1_V_W$I0, posterior_sample_table$V0[i], data_combined_muc1_V_W$MR0, data_combined_muc1_V_W$MA0)
  param_fit = c(posterior_sample_table$epsilon1[i],
                posterior_sample_table$k1[i],
                posterior_sample_table$beta[i], 
                posterior_sample_table$delta_I[i], 
                posterior_sample_table$p[i], 
                posterior_sample_table$c[i],
                posterior_sample_table$kappa_M[i],
                posterior_sample_table$s[i],
                posterior_sample_table$delta_M[i],
                posterior_sample_table$epsilon2[i],
                posterior_sample_table$k2[i],
                posterior_sample_table$q_V[i],
                posterior_sample_table$eta[i])
  
  param_fit_hat = c(posterior_sample_table$epsilon1[i] * 0,
                    posterior_sample_table$k1[i],
                    posterior_sample_table$beta[i], 
                    posterior_sample_table$delta_I[i], 
                    posterior_sample_table$p[i], 
                    posterior_sample_table$c[i],
                    posterior_sample_table$kappa_M[i],
                    posterior_sample_table$s[i],
                    posterior_sample_table$delta_M[i],
                    posterior_sample_table$epsilon2[i] * 0,
                    posterior_sample_table$k2[i],
                    posterior_sample_table$q_V[i],
                    posterior_sample_table$eta[i])
  
  
  model_output_WT = ode(times = t_ppc, y = y_init, func = Within_host_model, parms = param_fit, method = "bdf")
  model_output_KO = ode(times = t_ppc, y = y_init, func = Within_host_model, parms = param_fit_hat, method = "bdf")
  
  V_ppc_WT[i,] = model_output_WT[,4]
  M_ppc_WT[i,] = model_output_WT[,6]
  
  V_ppc_KO[i,] = model_output_KO[,4]
  M_ppc_KO[i,] = model_output_KO[,6]
}






for (i in 1:length(t_ppc)){
  temp_WT = unname(quantile(V_ppc_WT[,i], probs = c(0.025, 0.5, 0.975)))
  temp_KO = unname(quantile(V_ppc_KO[,i], probs = c(0.025, 0.5, 0.975)))
  
  
  temp_WT_M = unname(quantile(M_ppc_WT[,i], probs = c(0.025, 0.5, 0.975)))
  temp_KO_M = unname(quantile(M_ppc_KO[,i], probs = c(0.025, 0.5, 0.975)))
  
  lower_95PI_WT[i] = temp_WT[1]
  median_95PI_WT[i] = temp_WT[2]
  upper_95PI_WT[i] = temp_WT[3]
  
  lower_95PI_KO[i] = temp_KO[1]
  median_95PI_KO[i] = temp_KO[2]
  upper_95PI_KO[i] = temp_KO[3]
  
  
  lower_95PI_WT_M[i] = temp_WT_M[1]
  median_95PI_WT_M[i] = temp_WT_M[2]
  upper_95PI_WT_M[i] = temp_WT_M[3]
  
  lower_95PI_KO_M[i] = temp_KO_M[1]
  median_95PI_KO_M[i] = temp_KO_M[2]
  upper_95PI_KO_M[i] = temp_KO_M[3]
  
}

data_plot = data.frame(time = rep(c(1,2,3,5,7), each = 5),
                       V_WT = log10(exp(data_combined_muc1_V_W$log_viral_load_data_WT)),
                       V_KO = log10(exp(data_combined_muc1_V_W$log_viral_load_data_KO)))

fit_plot = data.frame(time = t_ppc,
                      lower_95PI_WT = lower_95PI_WT,
                      V_WT = median_95PI_WT,
                      upper_95PI_WT = upper_95PI_WT,
                      lower_95PI_KO = lower_95PI_KO,
                      V_KO = median_95PI_KO,
                      upper_95PI_KO = upper_95PI_KO)

ggplot(fit_plot, aes(time))+
  geom_point(data = data_plot, aes(time, V_WT), size = 3) +
  geom_line(data = fit_plot, aes(time, log10(V_WT))) +
  geom_ribbon(aes(ymin = log10(lower_95PI_WT), ymax = log10(upper_95PI_WT)), alpha = 0.5, na.rm = TRUE) +
  theme_bw() # V_WT

ggplot(fit_plot, aes(time))+
  geom_point(data = data_plot, aes(time, V_KO), size = 3) +
  geom_line(data = fit_plot, aes(time, log10(V_KO))) +
  geom_ribbon(aes(ymin = log10(lower_95PI_KO), ymax = log10(upper_95PI_KO)), alpha = 0.5, na.rm = TRUE) +
  theme_bw() # V_KO



data_plot_M_WT = data.frame(time = rep(c(1,3,5,7), c(5,5,5,4)),
                            M_WT = log10(exp(data_combined_muc1_V_W$Macrophage_data_WT)))


fit_plot_M_WT = data.frame(time = t_ppc,
                           lower_95PI_WT = lower_95PI_WT_M,
                           M_WT = median_95PI_WT_M,
                           upper_95PI_WT = upper_95PI_WT_M)

ggplot(fit_plot_M_WT, aes(time))+
  geom_point(data = data_plot_M_WT, aes(time, M_WT), size = 3) +
  geom_line(data = fit_plot_M_WT, aes(time, log10(M_WT))) +
  geom_ribbon(aes(ymin = log10(lower_95PI_WT), ymax = log10(upper_95PI_WT)), alpha = 0.5, na.rm = TRUE) +
  theme_bw() # M_WT





data_plot_M_KO = data.frame(time = rep(c(1,3,5,7), c(5,4,5,6)),
                            M_KO = log10(exp(data_combined_muc1_V_W$Macrophage_data_KO)))


fit_plot_M_KO = data.frame(time = t_ppc,
                           lower_95PI_KO = lower_95PI_KO_M,
                           M_KO = median_95PI_KO_M,
                           upper_95PI_KO = upper_95PI_KO_M)

ggplot(fit_plot_M_KO, aes(time))+
  geom_point(data = data_plot_M_KO, aes(time, M_KO), size = 3) +
  geom_line(data = fit_plot_M_KO, aes(time, log10(M_KO))) +
  geom_ribbon(aes(ymin = log10(lower_95PI_KO), ymax = log10(upper_95PI_KO)), alpha = 0.5, na.rm = TRUE) +
  theme_bw() # M_KO


# =============== Save the fitting results and generate outputs in CSV file  ===================
write.csv(posterior_samples_merged_after_burnin, file="Muc1_fit_full.csv")
save.image(file = "fitting_results_Mu1")
