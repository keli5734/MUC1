functions{
  real[]  Muc1_func(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    
    real dydt[4];
    
    dydt[1] = -(1 - theta[1] * t^10/(t^10 + theta[2]^10)) * theta[3] * y[1] * y[3];
    dydt[2] =  (1 - theta[1] * t^10/(t^10 + theta[2]^10)) * theta[3] * y[1] * y[3] - theta[4] * y[2];
    dydt[3] = theta[5] * y[2] - theta[6] * y[3] - theta[7] * y[3] * (y[4] + (theta[8]*y[3]*3e+2)/((theta[9] + 1e-2 + theta[8]*y[3])*1e-2));
    dydt[4] = theta[10] * (1 - theta[11] * t^10/(t^10 + theta[12]^10)) * (1 - y[4]/1.8e+5) * ((theta[8]*y[3]*3e+2)/((theta[9] + 1e-2 + theta[8]*y[3])*1e-2))  - 0.5 * y[4];
    
    return dydt;
  }
}

data{
  int<lower = 0> N_T_WT; // length of time data WT
  int<lower = 0> N_V_WT; // length of viral load WT
  real time_data_WT[N_T_WT]; // time series data  WT
  real log_viral_load_data_WT[N_V_WT]; // viral load data WT
  
  int<lower = 0> N_T_KO; // length of time data KO
  int<lower = 0> N_V_KO; // length of viral load KO
  real time_data_KO[N_T_KO]; // time series data KO
  real log_viral_load_data_KO[N_V_KO]; // viral laod data KO
  
  int<lower = 0> N_T_Macrophage_WT;
  int<lower = 0> N_Macrophage_WT;
  real time_data_Macrophage_WT[N_T_Macrophage_WT];
  real Macrophage_data_WT[N_Macrophage_WT];
  
  int<lower = 0> N_T_Macrophage_KO;
  int<lower = 0> N_Macrophage_KO;
  real time_data_Macrophage_KO[N_T_Macrophage_KO];
  real Macrophage_data_KO[N_Macrophage_KO];
  

  real t0; // initial time point 
  real T0; // initial target cells
  real I0; // initial infected cells 
  real MP0; // initial M_P
}

transformed data{
  real x_r[0];
  int x_i[0];
}

parameters{
  real<lower=0> theta[13]; // theta[1:15] = c(epsilon1, k1, beta, delta_I, p, delta_V, kappa_M, q_V, eta, phi, epsilon2, k2, V0)
  real<lower=0> sigma[2];
}

transformed parameters{
  real y_hat_WT[N_T_WT,4];
  real Y0[4];
  vector[N_V_WT] V_pred_WT; 
  
  real y_hat_KO[N_T_KO,4];
  vector[N_V_KO] V_pred_KO;
  real theta_hat[13];
  
  real y_hat_WT_Macrophage[N_T_Macrophage_WT,4];
  vector[N_Macrophage_WT] Macrophage_pred_WT;
  
  real y_hat_KO_Macrophage[N_T_Macrophage_KO,4];
  vector[N_Macrophage_KO] Macrophage_pred_KO;

  
  
  Y0[1] = T0;
  Y0[2] = I0;
  Y0[3] = theta[13];
  Y0[4] = MP0;
  
  theta_hat[1] = 0;
  
  for (j in 2:10){
    theta_hat[j] = theta[j];
  }
  theta_hat[11] = 0;
  theta_hat[12] = theta[12];
  theta_hat[13] = theta[13];
  

  
  y_hat_WT = integrate_ode_bdf(Muc1_func, Y0, t0, time_data_WT, theta, x_r, x_i);
  for(i in 1:N_T_WT){
    V_pred_WT[(i-1)*5+1:i*5] = rep_vector(y_hat_WT[i,3],5);
  }
  
  y_hat_KO = integrate_ode_bdf(Muc1_func, Y0, t0, time_data_KO, theta_hat, x_r, x_i);
  for(i in 1:N_T_KO){
    V_pred_KO[(i-1)*5+1:i*5] = rep_vector(y_hat_KO[i,3],5);
  }
  
  
  
  y_hat_WT_Macrophage = integrate_ode_bdf(Muc1_func, Y0, t0, time_data_Macrophage_WT, theta, x_r, x_i);
  Macrophage_pred_WT[1:5] = rep_vector(y_hat_WT_Macrophage[1,4]+3e+4,5);
  Macrophage_pred_WT[6:10] = rep_vector(y_hat_WT_Macrophage[2,4]+3e+4,5);
  Macrophage_pred_WT[11:15] = rep_vector(y_hat_WT_Macrophage[3,4]+3e+4,5);
  Macrophage_pred_WT[16:19] = rep_vector(y_hat_WT_Macrophage[4,4]+3e+4,4);
  
  
  y_hat_KO_Macrophage = integrate_ode_bdf(Muc1_func, Y0, t0, time_data_Macrophage_KO, theta_hat, x_r, x_i);
  Macrophage_pred_KO[1:5] = rep_vector(y_hat_KO_Macrophage[1,4]+3e+4,5);
  Macrophage_pred_KO[6:9] = rep_vector(y_hat_KO_Macrophage[2,4]+3e+4,4);
  Macrophage_pred_KO[10:14] = rep_vector(y_hat_KO_Macrophage[3,4]+3e+4,5);
  Macrophage_pred_KO[15:20] = rep_vector(y_hat_KO_Macrophage[4,4]+3e+4,6);

}

model{
// priors
theta[1] ~ normal(0.5,0.3); // epsilon1
theta[2] ~ normal(0.5,1); // k1
theta[3] ~ lognormal(log(2.5e-5),1); // beta
theta[4] ~ lognormal(log(1.3),1); // delta_I
theta[5] ~ lognormal(log(5e-3),1); // p
theta[6] ~ lognormal(log(1.3),1); // delta_V
theta[7] ~ lognormal(log(1e-7),1); // kappa_M
theta[8] ~ lognormal(log(1e-6),1); // q_V
//theta[9] ~ lognormal(log(3e+2),1); // s
theta[9] ~ lognormal(log(1),1); // eta
//theta[11] ~ lognormal(log(1e-2),1); // delta_M
theta[10] ~ lognormal(log(120), 1); // phi 
theta[11] ~ normal(0.5,0.3); // epsilon2
theta[12] ~ normal(0.5,1); // k2
theta[13] ~ lognormal(log(3.5e-1),1); // V0
sigma[1] ~ normal(0,1);
sigma[2] ~ normal(0,1);


//for (i in 1: N_V_WT){
//  target += normal_lpdf(log(log_viral_load_data_WT[i]) | log(V_pred_WT[i]), sigma[1]);
//}

//for (i in 1: N_V_KO){
//  target += normal_lpdf(log(log_viral_load_data_KO[i]) | log(V_pred_KO[i]), sigma[1]);
//}

//for (i in 1: N_TNF_WT){
//  target += normal_lpdf(TNF_data_WT[i] | TNF_pred_WT[i], sigma[2]);
//}

//for (i in 1: N_TNF_KO){
//  target += normal_lpdf(TNF_data_KO[i] | TNF_pred_KO[i], sigma[2]);
//}

log_viral_load_data_WT ~ normal(log(V_pred_WT), sigma[1]); // --> measured 
log_viral_load_data_KO ~ normal(log(V_pred_KO), sigma[1]); 

Macrophage_data_WT ~ normal(log(Macrophage_pred_WT), sigma[2]);
Macrophage_data_KO ~ normal(log(Macrophage_pred_KO), sigma[2]);


}




