
#rm(list=ls())  # clear memory
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
                                T0 = 1e+7,
                                I0 = 0,
                                V0 = 30)



init1 <- list(
  theta = c(1, 1e-5, 2, 10, 
            2.3, 1e-7, 4e+2,1,
            2, 2e-2),
  sigma = c(0.5, 0.5))

init2 <- list(
  theta = c(0.5, 8e-6, 0.5, 3, 
            1, 1e-7, 4e+2,2,
            1.5, 2e-2),
  sigma = c(0.5, 0.5))

## =========== Viral parameters only ==============


init1 <- list(
  theta = c(0.5, 0.5, 8e-6, 0.5, 3, 3, 400),
  sigma = c(0.5, 0.5))


init2 <- list(
  theta = c(0.2, 1, 4e-6, 4, 10, 20),
  sigma = c(0.5, 1))


init3 <- list(
  theta = c(0.7, 0.1, 1e-6, 2, 50, 40),
  sigma = c(1, 0.5))


init4 <- list(
  theta = c(0.1, 1.5, 8e-5, 3, 20, 50),
  sigma = c(1, 1))

init5 <- list(
  log10_theta = c(log10(1.5), log10(8e-6), log10(3), log10(20), log10(20),
            log10(1e-5), log10(400), log10(1.5), log10(1e-2)),
  sigma = c(1, 1))

options(mc.cores=parallel::detectCores()) # to utilise all cores available in your computer

Estimate_Full_06_OCT <- stan("PriorStan2.stan",
                             data = data_combined_muc1_V_W,
                             pars = c("theta_KO","sigma"),
                             seed = 20201910,  # set random seed for reproducibility
                             iter = 2000,
                             chains = 1,
                             init = list(init5),
                             warmup = 1000,
                             control = list(adapt_delta = 0.99, max_treedepth = 15))


#print(Estimate_Full_06_OCT, pars = c("theta_WT"))
print(Estimate_Full_06_OCT, pars = c("theta_KO"))
stan_dens(Estimate_Full_06_OCT, pars = c("theta_KO"), separate_chains = TRUE,nrow = 6)



# ASSESSING AND FIXING DIVERGENCES AND TREEDEPTH PROBLEMS

mack_diagnostics <- rstan::get_sampler_params(Estimate_Full_06_OCT) %>% 
  set_names(1) %>% 
  map_df(as_data_frame, .id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= 1000)

mack_diagnostics %>% 
  group_by(warmup, chain) %>% 
  summarise(percent_divergent = mean(divergent__ > 0)) %>% 
  ggplot()+ 
  geom_col(aes(chain, percent_divergent, fill = warmup), position = 'dodge', color = 'black') + 
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs") + 
  theme_bw() # plot divergent rate for each chain 

mack_diagnostics %>% 
  ggplot(aes(iteration, treedepth__, color = chain)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 15), color = 'red') + 
  theme_bw() # plot treedepth for each chain 


mack_diagnostics %>% 
  ggplot(aes(iteration, stepsize__, color = chain)) + 
  geom_line() + 
  lims(y = c(0,.1)) + 
  theme_bw() # plot stepsize for each chain 


# PARAMETER DIAGNOSTICS
para_summary <- summary(Estimate_Full_06_OCT)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

# n_eff 
para_summary %>% 
  ggplot(aes(n_eff)) + 
  geom_histogram(binwidth = 1000/50) + 
  geom_vline(aes(xintercept = 1000), color = 'red')

# R_hat

para_summary %>% 
  ggplot(aes(Rhat)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = 1.01), color = 'red')

# parameter 

para_summary %>% 
  filter(variable %in% c('theta[1]', 'theta[2]', 'theta[3]','theta[4]', 'theta[5]', 'epsilon')) %>% 
  ggplot() + 
  geom_linerange(aes(variable, ymin = `2.5%`,ymax = `97.5%`)) + 
  geom_crossbar(aes(variable, mean, ymin = `25%`, ymax = `75%`), fill= 'grey') + 
  facet_wrap(~variable, scales = 'free') + 
  theme_bw()


# plot sampling parameters 
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

params_cp <- as.data.frame(rstan::extract(Estimate_Full_06_OCT, permuted=FALSE))

names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:1000

plot(params_cp$iter, params_cp$theta.10, col=c_dark, pch=16, cex=0.8,
     xlab="Iteration", ylab="epsilon2 (Chain1)", ylim=c(0,1))
plot(params_cp$iter, params_cp$`chain:2.theta.10`, col=c_dark, pch=16, cex=0.8,
     xlab="Iteration", ylab="epsilon2 (Chain2)", ylim=c(0,1))

running_means <- sapply(params_cp$iter, function(n) mean(params_cp$theta.10[1:n]))
plot(params_cp$iter, running_means, col=c_dark, pch=16, cex=0.8, ylim=c(0, 1),
     xlab="Iteration", ylab="epsilon2")



divergent <- get_sampler_params(Estimate_Full_06_OCT, inc_warmup=FALSE)[[1]][,'divergent__']
params_cp$divergent <- divergent

div_params_cp <- params_cp[params_cp$divergent == 1,]
nondiv_params_cp <- params_cp[params_cp$divergent == 0,]

plot(div_params_cp$theta.9, div_params_cp$theta.10,
     col="green", pch=16, cex=0.8, xlab="phi", ylab="epsilon2",
     xlim=c(0,50), ylim=c(0,1))
points(nondiv_params_cp$theta.9, nondiv_params_cp$theta.10,
       col=c_dark, pch=16, cex=0.8)


# show mcmc_parcoord
draws <- as.array(Estimate_Full_06_OCT,pars = c("theta_KO"))
np <- nuts_params(Estimate_Full_06_OCT)
str(np)

color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.9)
mcmc_parcoord(
  draws,
  transform = function(x) {(x - mean(x)) / sd(x)},
  size = 0.25,
  alpha = 0.8,
  np = np,
  np_style = div_style
)
d <- mcmc_parcoord_data(draws, np = np)
head(d)
tail(d)


# extract posterior samples for selected parameters
posterior_samples_all = rstan::extract(Estimate_Full_06_OCT, pars = c("theta_KO",'sigma'), inc_warmup = TRUE, permuted = FALSE)
posterior_samples_merged_after_burnin = rstan::extract(Estimate_Full_06_OCT, pars = c( "theta_KO",'sigma'))

color_scheme_set("brewer-Spectral")
mcmc_trace(posterior_samples_all, n_warmup = 1000,
           facet_args = list(nrow = 4, labeller = label_parsed))

# show all marginal posterior distributions
posterior_sample_table = data.frame(epsilon1 = posterior_samples_merged_after_burnin$theta_KO[,1],
                                    #k1 = posterior_samples_merged_after_burnin$theta[,2],
                                    beta = posterior_samples_merged_after_burnin$theta_KO[,2],
                                    delta_I = posterior_samples_merged_after_burnin$theta_KO[,3],
                                    p = posterior_samples_merged_after_burnin$theta_KO[,4],
                                    delta_V = posterior_samples_merged_after_burnin$theta_KO[,5],
                                    kappa_M = posterior_samples_merged_after_burnin$theta_KO[,6],
                                    s = posterior_samples_merged_after_burnin$theta_KO[,7],
                                    #phi = posterior_samples_merged_after_burnin$theta[,8],
                                    epsilon2 = posterior_samples_merged_after_burnin$theta_KO[,8],
                                    #k2 = posterior_samples_merged_after_burnin$theta[,9],
                                    delta_M = posterior_samples_merged_after_burnin$theta[,9])



# ================== plot posterior and prior for each parameter =================

ggplot(posterior_sample_table, aes(x = log10(epsilon1))) + 
  geom_histogram(breaks=seq(-2,2,2/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dunif(x, -1, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(epsilon1, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(epsilon1))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(epsilon1, 0.975))),
             linetype="dashed",color = 'red') + 
  #geom_vline(aes(xintercept=log10(1)),
  #           linetype="dashed", color = 'black') + 
  lims(x = c(-1,1)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # epsilon1 


beta_WT <-  posterior_sample_table$beta
beta_KO <-  posterior_sample_table$beta * posterior_sample_table$epsilon1
beta.df <- data.frame(beta = log10(c(beta_WT, beta_KO)), type = as.factor(rep(c('WT', 'KO'), each = length(beta_WT))))
library(plyr)
mu <- ddply(beta.df, "type", summarise, grp.mean=mean(beta))
head(mu)
p <- ggplot(beta.df, aes(x=beta, color=type)) +
  geom_histogram(aes(y = ..density..),fill = "white", position="dodge", alpha = 1, binwidth = 1/30)+
  theme(legend.position="top")+
  #scale_color_grey()+scale_fill_grey() + 
  theme_classic() 

p +  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),
                linetype="dashed") +
  stat_function(fun = function(x) {dnorm(x, -6, 1)}, aes(color = 'prior')) +
  lims(x = c(-7,-4)) + xlab('log10(beta)')          # beta



ggplot(posterior_sample_table, aes(x = log10(delta_I))) + 
  geom_histogram(breaks=seq(-3,3,0.05),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, 0, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(delta_I, 0.025))),
               linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(delta_I))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(delta_I, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-1, 3)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # delta_I





ggplot(posterior_sample_table, aes(x = log10(p))) + 
  geom_histogram(breaks=seq(-6,2,8/100),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, -2, 2)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(p, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(p))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(p, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-6, 3)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # p




ggplot(posterior_sample_table, aes(x = log10(delta_V))) + 
  geom_histogram(breaks=seq(-2,2,4/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, 0, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(delta_V, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(delta_V))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(delta_V, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-3,2)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # delta_V


ggplot(posterior_sample_table, aes(x = log10(kappa_M))) + 
  geom_histogram(breaks=seq(-8,-4,4/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, -6, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(kappa_M, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(kappa_M))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(kappa_M, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-8,-4)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="white", prior="purple")) # kappa_M

V.decay <- posterior_sample_table$kappa_M * posterior_sample_table$s/posterior_sample_table$delta_M + posterior_sample_table$delta_V
V.decay.df <- data.frame(V.decay = V.decay)
ggplot(V.decay.df, aes(x = log10(V.decay))) + 
  geom_histogram(breaks=seq(-2,1,3/50), aes(y = ..density..,color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, 0, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(V.decay, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(V.decay))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(V.decay, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-3,1)) + 
  theme_classic() +
  scale_colour_manual(name="Distribution",
                     values=c(posterior="black", prior="purple")) # V.decay



ggplot(posterior_sample_table, aes(x = log10(s))) + 
  geom_histogram(breaks=seq(1,3.5,2.5/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, 2, 0.5)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(s, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(s))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(s, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(0,4)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # s


ggplot(posterior_sample_table, aes(x = log10(epsilon2))) + 
  geom_histogram(breaks=seq(-2,2,2/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dunif(x, -1, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(epsilon2, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(epsilon2))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(epsilon2, 0.975))),
             linetype="dashed",color = 'red') + 
  lims(x = c(-1,1)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # epsilon2




ggplot(posterior_sample_table, aes(x = log10(delta_M))) + 
  geom_histogram(breaks=seq(-4,0,4/50),aes(y = ..density.., color = 'posterior')) + 
  stat_function(fun = function(x) {dnorm(x, -2, 1)}, aes(color = 'prior')) +
  geom_vline(aes(xintercept=log10(quantile(delta_M, 0.025))),
             linetype="dashed",color = 'red') + 
  geom_vline(aes(xintercept=log10(mean(delta_M))),
             linetype="solid",color = 'red') + 
  geom_vline(aes(xintercept=log10(quantile(delta_M, 0.975))),
             linetype="dashed",color = 'red') +  
  lims(x = c(-4,0)) + 
  theme_classic() + 
  scale_colour_manual(name="Distribution",
                      values=c(posterior="black", prior="purple")) # delta_M





# show parameter correlations
pairs(posterior_samples_merged_after_burnin$theta_KO[,1:9], col = c("blue", "red"), labels = c('epsilon1','beta','delta_I','p','delta_V','kappa_M','s','epsilon2','delta_M'))#[,1:7], labels = c("epsilon1","k1","beta","delta_I", "p", "c","kappa_M"), col = c( "blue", "red"))
pairs(posterior_samples_merged_after_burnin$theta[,9:12], col = c("blue", "red"), labels = c('phi','epsilon2','k2','delta_M'))#[,1:7], labels = c("epsilon1","k1","beta","delta_I", "p", "c","kappa_M"), col = c( "blue", "red"))



# posterior predictive check
Within_host_model = function(t, y, theta){
  
  dydt1 = -theta[1] * theta[2] * y[1] * y[3];
  dydt2 =  theta[1] * theta[2] * y[1] * y[3] - theta[3] * y[2];
  dydt3 =  theta[4] * y[2] - theta[5] * y[3] - theta[6] * y[4] * y[3];
  dydt4 =  theta[7] + theta[8] * y[2] - theta[9] * y[4];
  
  list(c(dydt1, dydt2, dydt3, dydt4))
  
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

R0_WT = c()
R0_KO = c()



for (i in 1:1000){
  y_init = c(data_combined_muc1_V_W$T0, data_combined_muc1_V_W$I0, data_combined_muc1_V_W$V0, posterior_sample_table$s[i]/posterior_sample_table$delta_M[i])
  param_fit = c(1,
                #posterior_sample_table$k1[i],
                posterior_sample_table$beta[i], 
                posterior_sample_table$delta_I[i], 
                posterior_sample_table$p[i], 
                posterior_sample_table$delta_V[i],
                posterior_sample_table$kappa_M[i],
                posterior_sample_table$s[i],
                1,
                #posterior_sample_table$k2[i],
                posterior_sample_table$delta_M[i])
  
  param_fit_hat = c(posterior_sample_table$epsilon1[i],
                #posterior_sample_table$k1[i],
                posterior_sample_table$beta[i], 
                posterior_sample_table$delta_I[i], 
                posterior_sample_table$p[i], 
                posterior_sample_table$delta_V[i],
                posterior_sample_table$kappa_M[i],
                posterior_sample_table$s[i],
                #posterior_sample_table$phi[i],
                posterior_sample_table$epsilon2[i],
                #posterior_sample_table$k2[i],
                posterior_sample_table$delta_M[i])
  
  model_output_WT = ode(times = t_ppc, y = y_init, func = Within_host_model, parms = param_fit, method = "bdf")
  model_output_KO = ode(times = t_ppc, y = y_init, func = Within_host_model, parms = param_fit_hat, method = "bdf")
  
  V_ppc_WT[i,] = model_output_WT[,4]
  M_ppc_WT[i,] = model_output_WT[,5] 
  
  V_ppc_KO[i,] = model_output_KO[,4]
  M_ppc_KO[i,] = model_output_KO[,5] 
  
  R0_WT[i] <-  posterior_sample_table$p[i] *  posterior_sample_table$beta[i] * 1e+7 / (posterior_sample_table$delta_I[i]*(posterior_sample_table$delta_V[i] + posterior_sample_table$kappa_M[i] * posterior_sample_table$s[i] / posterior_sample_table$delta_M[i]))
  R0_KO[i] <-  posterior_sample_table$epsilon1[i] * posterior_sample_table$p[i] *  posterior_sample_table$beta[i] * 1e+7 / (posterior_sample_table$delta_I[i]*(posterior_sample_table$delta_V[i] + posterior_sample_table$kappa_M[i] * posterior_sample_table$s[i] / posterior_sample_table$delta_M[i]))
  
}

# ======== plot R0 ============


R0.df <- data.frame(R0 = c(R0_WT, R0_KO), type = as.factor(rep(c('WT', 'KO'), each = length(R0_KO))))
mu_R0 <- ddply(R0.df, "type", summarise, grp.mean=mean(R0))
head(mu_R0)
p <- ggplot(R0.df, aes(x=R0, color=type)) +
  geom_histogram(aes(y = ..density..),fill = "white", position="dodge", alpha = 1, binwidth = 1)+
  theme(legend.position="top")+
  theme_classic() 

p +  geom_vline(data=mu_R0, aes(xintercept=grp.mean, color=type),
                linetype="dashed") +
  lims(x = c(1,50)) + xlab('R0')         #R0

pair.comp.df <- data.frame(R0 = R0_KO, 
                           epsilon1 = posterior_samples_merged_after_burnin$theta_KO[,1], 
                           epsilon2 = posterior_samples_merged_after_burnin$theta_KO[,8],
                           M_number = posterior_samples_merged_after_burnin$theta_KO[,7]/posterior_samples_merged_after_burnin$theta_KO[,9])
pairs(pair.comp.df, col = c("blue", "red"), labels = c('R0','epsilon1','epsilon2','M_number'))



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
  ylim(log10(1e+3), log10(1e+6))+
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
  ylim(log10(1e+3), log10(1e+6))+
  theme_bw() # M_KO


# =============== Save the fitting results and generate outputs in CSV file  ===================
write.csv(posterior_samples_merged_after_burnin, file="Muc1_fit_OCT_21.csv")
save.image(file = "fitting_results_Mu1")












