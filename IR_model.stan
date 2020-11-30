functions{
  real[]  Muc1_func(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    
    
    //int n_TIV_innate;
    //int n_E;
    //int n_p;
    //int tau_E;
    //int tau_p;
    //int n_total;
    real dydt[20];
    
    //n_TIV_innate = 4;
    //n_E = 200;
    //n_p = 200;
    //tau_E = 6;
    //tau_p = 4;
    //n_total = n_TIV_innate + 1 + 1 + (n_E-1) + 1 + 1 + 1 + (n_p-1) + 1 + 1 + 1;
    
    
    dydt[1] =  0.8 * y[1] * (1 - (y[1] + y[2])/ 1e+7) - (1 - theta[1]* y[2]/(y[2] + theta[10]) ) * theta[2] * y[1] * y[3];
    dydt[2] =  (1 - theta[1]* y[2]/(y[2] + theta[10]) ) * theta[2] * y[1] * y[3] - theta[3] * y[2] - 2e-5 * y[11] * y[2];
    dydt[3] =  theta[4] * y[2] - theta[5] * y[3] - theta[6] * y[4] * y[3] - 0.8 * y[19] * y[3] - 0.4 * y[20] * y[3];
    dydt[4] =  330 + (1 - theta[7]* y[2]/(y[2] + theta[10]) ) * theta[9] * y[2] - theta[8] * y[4];
    dydt[5] = -2*y[3]/(y[3] + 1e+4) * y[5];
    dydt[6] = 2*y[3]/(y[3] + 1e+4) * y[5] - 5 * 1.0 / 10 * y[6];
    
    for (i in 7 : 10){
      dydt[i] = 5 * 1.0 / 10 * (y[i-1] - y[i]);
    }
    dydt[11] = 1400 * 5 * 1.0 / 10 * y[10] - 0.57 * y[11];
    dydt[12] = -0.06 * y[3]/(y[3] + 1e+4) * y[12];
    dydt[13] = 0.06 * y[3]/(y[3] + 1e+4) * y[12] - 5 * 1.0 / 8 * y[13];
    
    for (j in 14 : 17){
      dydt[j] = 5 * 1.0 / 8 * (y[j-1] - y[j]);
    }
    dydt[18] = 8 * 5 * 1.0 / 8 * y[17] - 0.5 * y[18];
    dydt[19] = 12 * y[18] - 2 * y[19]; 
    dydt[20] = 4 * y[18] - 0.015 * y[20];
     
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
  real V0; // initial viral inoculation
  real E_naive0; // 
  real E0[6]; //
  real B_naive0; 
  real B0[6]; 
  real A0;
  real AL0;
}

transformed data{
  real x_r[0];
  int x_i[0];
}

parameters{
  
  real log10_theta[5]; // use log10() normal priors
  real<lower=0> theta[5]; // use lognormal priros delta_I, delta_V, delta_M and uniform distributions for epsilon1 and epsilon2
  real<lower=0> sigma[2];
}

transformed parameters{
  real y_hat_WT[N_T_WT,20];
  vector[N_V_WT] V_pred_WT; 
  real Y0_WT[20];
  real theta_WT[10];
  
  real y_hat_KO[N_T_KO,20];
  vector[N_V_KO] V_pred_KO;
  real Y0_KO[20];
  real theta_KO[10];
  
  real y_hat_WT_Macrophage[N_T_Macrophage_WT,20];
  vector[N_Macrophage_WT] Macrophage_pred_WT;
  
  real y_hat_KO_Macrophage[N_T_Macrophage_KO,20];
  vector[N_Macrophage_KO] Macrophage_pred_KO;

  
  theta_WT[1] = theta[1]; // epsilon1 
  theta_WT[2] = 10^log10_theta[1]; // beta 
  theta_WT[3] = theta[2]; // delta_I
  theta_WT[4] = 10^log10_theta[2]; // p 
  theta_WT[5] = theta[3]; // delta_V
  theta_WT[6] = 10^log10_theta[3]; // kappa_M;
  theta_WT[7] = theta[4]; // epsilon2
  theta_WT[8] = theta[5]; // delta_M
  theta_WT[9] = 10^log10_theta[4]; // phi 
  theta_WT[10] = 10^log10_theta[5]; // I_50

  theta_KO[1] = 0;
  theta_KO[2] = 10^log10_theta[1];
  theta_KO[3] = theta[2];
  theta_KO[4] = 10^log10_theta[2];
  theta_KO[5] = theta[3];
  theta_KO[6] = 10^log10_theta[3];
  theta_KO[7] = 0;
  theta_KO[8] = theta[5];
  theta_KO[9] = 10^log10_theta[4];
  theta_KO[10] = 10^log10_theta[5];

  
  
  Y0_WT[1] = T0;
  Y0_WT[2] = I0;
  Y0_WT[3] = V0;
  Y0_WT[4] = 330/theta_WT[8];
  Y0_WT[5] = E_naive0;
  for (i in 6:11){
    Y0_WT[i] = 0;
  }
  Y0_WT[12] = B_naive0;
  for(i in 13:20){
    Y0_WT[i] = 0;
  }

  Y0_KO[1] = T0;
  Y0_KO[2] = I0;
  Y0_KO[3] = V0;
  Y0_KO[4] = 330/theta_KO[8];
  Y0_KO[5] = E_naive0;
  for (i in 6:11){
    Y0_KO[i] = 0;
  }
  Y0_KO[12] = B_naive0;
  for(i in 13:20){
    Y0_KO[i] = 0;
  }

  
  y_hat_WT = integrate_ode_bdf(Muc1_func, Y0_WT, t0, time_data_WT, theta_WT, x_r, x_i);
  for(i in 1:N_T_WT){
    V_pred_WT[(i-1)*5+1:i*5] = rep_vector(y_hat_WT[i,3],5);
  }
  
  y_hat_KO = integrate_ode_bdf(Muc1_func, Y0_KO, t0, time_data_KO, theta_KO, x_r, x_i);
  for(i in 1:N_T_KO){
    V_pred_KO[(i-1)*5+1:i*5] = rep_vector(y_hat_KO[i,3],5);
  }
  
  
  
  y_hat_WT_Macrophage = integrate_ode_bdf(Muc1_func, Y0_WT, t0, time_data_Macrophage_WT, theta_WT, x_r, x_i);
  Macrophage_pred_WT[1:5] = rep_vector(y_hat_WT_Macrophage[1,4],5);
  Macrophage_pred_WT[6:10] = rep_vector(y_hat_WT_Macrophage[2,4],5);
  Macrophage_pred_WT[11:15] = rep_vector(y_hat_WT_Macrophage[3,4],5);
  Macrophage_pred_WT[16:19] = rep_vector(y_hat_WT_Macrophage[4,4],4);
  
  
  y_hat_KO_Macrophage = integrate_ode_bdf(Muc1_func, Y0_KO, t0, time_data_Macrophage_KO, theta_KO, x_r, x_i);
  Macrophage_pred_KO[1:5] = rep_vector(y_hat_KO_Macrophage[1,4],5);
  Macrophage_pred_KO[6:9] = rep_vector(y_hat_KO_Macrophage[2,4],4);
  Macrophage_pred_KO[10:14] = rep_vector(y_hat_KO_Macrophage[3,4],5);
  Macrophage_pred_KO[15:20] = rep_vector(y_hat_KO_Macrophage[4,4],6);

}

model{
// priors
//log10_theta[1] ~ normal(0,2); // epsilon1
//log10_theta[2] ~ normal(-6,4); // beta
//theta[1] ~ lognormal(log(0.89), 1); // delta_I
//log10_theta[3] ~ normal(-2,4); // p
//theta[2] ~ lognormal(log(28.4), 1); // delta_V
//log10_theta[4] ~ normal(-6,4); // kappa_M
//log10_theta[5] ~ normal(0,2); // epsilon2
//theta[3] ~ lognormal(log(4.2e-3), 1); // delta_M
//log10_theta[6] ~ normal(0,5); // converting factor of I. 


theta[1] ~ uniform(0,1); // epsilon1 
log10_theta[1] ~ normal(-6,4); // beta
theta[2] ~ lognormal(log(0.89), 1); // delta_I
log10_theta[2] ~ normal(-2,4); //p 
theta[3] ~ lognormal(log(28.4), 1); // delta_V
log10_theta[3] ~ normal(-6,4); //kappa_M
theta[4] ~ uniform(0,1); // epsilon2 
theta[5] ~ lognormal(log(4.2e-3), 1); // delta_M
log10_theta[4] ~ normal(0,4); // phi
log10_theta[5] ~ normal(0,5); // I_50


sigma[1] ~ normal(0,1);
sigma[2] ~ normal(0,1);



//log_viral_load_data_WT ~ normal(log(V_pred_WT), sigma[1]); // --> measured 
//log_viral_load_data_KO ~ normal(log(V_pred_KO), sigma[1]); 

//Macrophage_data_WT ~ normal(log(Macrophage_pred_WT), sigma[2]);
//Macrophage_data_KO ~ normal(log(Macrophage_pred_KO), sigma[2]);


for (i in 1: N_V_WT){
  target += normal_lpdf(log_viral_load_data_WT[i] | log(V_pred_WT[i]), sigma[1]);
  }
for (i in 1: N_V_KO){
  target += normal_lpdf(log_viral_load_data_KO[i] | log(V_pred_KO[i]), sigma[1]);
  }

for (i in 1: N_Macrophage_WT){
  target += normal_lpdf(Macrophage_data_WT[i] | log(Macrophage_pred_WT[i]), sigma[2]);
  }
for (i in 1: N_Macrophage_KO){
  target += normal_lpdf(Macrophage_data_KO[i] | log(Macrophage_pred_KO[i]), sigma[2]);
  }

}




