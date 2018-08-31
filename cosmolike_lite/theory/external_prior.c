double log_L_BAO();
double log_L_SN();
double log_L_Planck();
double log_L_Planck_BAO_SN();
double log_L_ia();
double log_L_clusterMobs();
double log_L_wlphotoz();
double log_L_clphotoz();
double log_L_shear_calib();
double log_like_f_red();
double do_matrix_mult_invcov(int n_param, double invcov[n_param][n_param], double param_diff[n_param]); //CH
double log_L_Planck15_BAO_w0wa(); //CH
double log_L_Planck15_BAO_H070p6_JLA_w0wa(); //CH
double log_L_Planck18_BAO_Riess18_Pantheon_w0wa(); //CH
double log_L_Planck18_BAO_w0wa(); //CH
double log_L_Planck18_w0(); //CH

double log_L_SRD_SN_Y1_RENEE();
double log_L_SRD_SN_Y10_RENEE();
double log_L_SRD_SL_Y1_TOM();
double log_L_SRD_SL_Y10_TOM();
double log_L_SRD_LSST_Y1_LSS();

void set_ia_priors();
void set_lin_bias_priors();
void set_HOD_redmagic_priors();
void set_baryon_priors();

void set_shear_priors_stage3();
void set_shear_priors_stage4();
void relax_shear_priors();

void set_wlphotoz_priors_stage3();
void set_wlphotoz_priors_stage4();
void set_clphotoz_priors_redmagic();
void set_clphotoz_priors_benchmark();
void set_clphotoz_priors_cmass();
void set_clphotoz_priors_source();

void set_clusterMobs_priors();



void set_HOD_redmagic_priors()
{
  flat_prior.HOD_rm[0][0] = 10.0;flat_prior.HOD_rm[0][1] = 15.0;
  flat_prior.HOD_rm[1][0] = 0.1;flat_prior.HOD_rm[1][1] = 1.0;
  flat_prior.HOD_rm[2][0] = 10.0;flat_prior.HOD_rm[2][1] = 15.0;
  flat_prior.HOD_rm[3][0] = 10.0;flat_prior.HOD_rm[3][1] = 15.0;
  flat_prior.HOD_rm[4][0] = 0.5;flat_prior.HOD_rm[4][1] = 1.5;
  flat_prior.fc_rm[0] = 0.1; flat_prior.fc_rm[1] = 1.0;
  flat_prior.cg_rm[0] = 0.2;flat_prior.cg_rm[1] = 5.;
}

void set_lin_bias_priors()
{
   printf("todo\n");
}

void set_baryon_priors()
{
   printf("todo\n");
}

void set_ia_priors()
{
  prior.LF_alpha[0]=0.;
  prior.LF_P[0]=0.;
  prior.LF_Q[0]=0.;
  prior.LF_red_alpha[0]=0.;
  prior.LF_red_P[0]=0.;
  prior.LF_red_Q[0]=0.;

  prior.LF_alpha[1]=0.05;
  prior.LF_P[1]=0.5;
  prior.LF_Q[1]=0.5;
  prior.LF_red_alpha[1]=0.1;
  prior.LF_red_P[1]=0.5;
  prior.LF_red_Q[1]=0.5;

}

void set_clusterMobs_priors()
{
  if (strcmp(survey.name,"LSST_Y10")==0 || strcmp(survey.name,"LSST_Y1")==0 ){
    // flat prior
    prior.cluster_Mobs_lgN0[0]=0.5;
    prior.cluster_Mobs_lgN0[1]=5.0;

    prior.cluster_Mobs_alpha[0]=0.0;
    prior.cluster_Mobs_alpha[1]=2.0;
    
    prior.cluster_Mobs_beta[0]=-1.5;
    prior.cluster_Mobs_beta[1]=1.5;
    
    prior.cluster_Mobs_sigma0[0]=0.0;
    prior.cluster_Mobs_sigma0[1]=1.5;
    
    prior.cluster_Mobs_sigma_qm[0] = -1.5;
    prior.cluster_Mobs_sigma_qm[1] = 1.5;

    prior.cluster_Mobs_sigma_qz[0] = -1.5;
    prior.cluster_Mobs_sigma_qz[1] = 1.5;

  } else {
    prior.cluster_Mobs_lgM0[0]=1.72+log(1.e+14*0.7);
    prior.cluster_Mobs_alpha[0]=1.08;
    prior.cluster_Mobs_beta[0]=0.0;
    prior.cluster_Mobs_sigma[0]=0.25;
    
    prior.cluster_Mobs_lgM0[1]=1.0;
    prior.cluster_Mobs_alpha[1]=0.2;
    prior.cluster_Mobs_beta[1]=0.5;
    prior.cluster_Mobs_sigma[1]=0.2;

    prior.cluster_completeness[0] = 0.9;
    prior.cluster_completeness[1] = 0.05;
  }
  
  like.clusterMobs=1;
}

void set_shear_priors_stage3()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 3\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.05;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}

void set_shear_priors_stage4()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.004;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}

void relax_shear_priors()
{
   int i;
   printf("Relaxing shear bias prior\n");
   for (i=0;i<tomo.shear_Nbin; i++){
      prior.shear_calibration_m[i][0] = 0.0;
      prior.shear_calibration_m[i][1] = 1.e3;
      printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
   }
   like.shearcalib=1;
}

void set_wlphotoz_priors_stage3()
{
  int i;  
  printf("Setting Gaussian shear photo-z Priors stage 3\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.01;
    prior.sigma_zphot_shear[i][1]=0.2*nuisance.sigma_zphot_shear[i];
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}

void set_wlphotoz_priors_SV()
{
  int i;  
  printf("Setting Gaussian shear photo-z Priors SV\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=0.0;
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.05;
    printf("Mean (of bias)=%le, Sigma (of bias)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
  }
  like.wlphotoz=1;
}

void set_wlphotoz_priors_stage4()
{
  int i;  
  printf("Setting Gaussian shear photo-z Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.002;
    prior.sigma_zphot_shear[i][1]=0.002; // http://lsst.org/files/docs/Phot-z-plan.pdf
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]);
  }
  like.wlphotoz=1;
}

void set_clphotoz_priors_redmagic()
{
  int i;
  printf("Setting Gaussian clustering photo-z Priors redmagic\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    //rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.0004;
    prior.sigma_zphot_clustering[i][1]=0.0006;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_priors_source()
{
  int i;
  printf("Setting Gaussian clustering photo-z Priors redmagic\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    //rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = prior.bias_zphot_shear[i][1];
    prior.sigma_zphot_clustering[i][1]= prior.sigma_zphot_shear[i][1];
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_priors_LSST_gold()
{
  int i;
  printf("Setting Gaussian clustering photo-z Priors redmagic\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    //rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.001;
    prior.sigma_zphot_clustering[i][1]=0.002;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_priors_LSST_SRD()
{
  int i;
  printf("Setting Gaussian clustering photo-z Priors redmagic\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    //rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.001;
    prior.sigma_zphot_clustering[i][1]=0.002;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_priors_benchmark()
{
  int i;
  printf("Setting Gaussian clustering photo-z Priors benchmark\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
  // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
  //rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.001;
    prior.sigma_zphot_clustering[i][1]=0.002; 
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}


void set_clphotoz_priors_cmass()
{
   printf("Assuming perfect spec-z for CMASS\n");
   like.clphotoz=0;
}



double log_L_BAO()
{
  double log_L = 0.;
  int i;
  for (i = 0; i < BAO.N; i++){
    //printf("%d %e %e %e\n", i,BAO.z[i], BAO.data[i],dist_BAO(BAO.z[i]));
    log_L -= pow((dist_BAO(BAO.z[i])-BAO.data[i])/(BAO.data[i]*BAO.sigma[i]),2.0);
  }
  return 0.5*log_L;
}

double log_L_Planck()
{
  double log_L = 0.;
  log_L-=(cosmology.Omega_m - 0.3156)*(cosmology.Omega_m - 0.3156)/(0.0091*0.0091);
  log_L-=(cosmology.sigma_8 - 0.831)*(cosmology.sigma_8 - 0.831)/(0.013*0.013);
  log_L-=(cosmology.omb - 0.0491685)*(cosmology.omb - 0.0491685)/(0.00035106195*0.00035106195);
  log_L-=(cosmology.n_spec - 0.9645)*(cosmology.n_spec - 0.9645)/(0.0049*0.0049);
  log_L-=(cosmology.h0 - 0.6727)*(cosmology.h0 - 0.6727)/(0.0066*0.0066);
  return 0.5*log_L;
}

double log_L_SN()
{
  double log_L = 0.;
  log_L-=(cosmology.Omega_m - 0.3156)*(cosmology.Omega_m - 0.3156)/(0.08856208*0.08856208);
  log_L-=(cosmology.h0 - 0.6727)*(cosmology.h0 - 0.6727)/(0.2492126*0.2492126);
  log_L-=(cosmology.w0+1.0)*(cosmology.w0+1.0)/(0.12428845*0.12428845);
  log_L-=(cosmology.wa - 0.0)*(cosmology.wa - 0.0)/(0.979977*0.979977);
  return 0.5*log_L;
}

double log_L_Planck_BAO_SN() // using the Planck 15 fid values and error bars as described in Aubourg et al 2014
{
  double log_L = 0.;
  double omegab=cosmology.omb*cosmology.h0*cosmology.h0;
  double omegab_fid=0.0491685*0.6727*0.6727;
  double w_pivot=cosmology.w0+(1.0-(1.0/(1.0+0.266)))*cosmology.wa; // see table 4 in aubourg et al 2014

  log_L-=(cosmology.Omega_m - 0.3156)*(cosmology.Omega_m - 0.3156)/(0.011*0.011);
  log_L-=(cosmology.h0 - 0.6727)*(cosmology.h0 - 0.6727)/(0.011*0.011);
  log_L-=(omegab - omegab_fid)*(omegab - omegab_fid)/(0.0003*0.0003);
  log_L-=(w_pivot +1.0)*(w_pivot +1.0)/(0.11*0.11);
  log_L-=(cosmology.wa - 0.0)*(cosmology.wa - 0.0)/(0.4*0.4);
 
 return 0.5*log_L;
 }




//CH
double do_matrix_mult_invcov(int n_param, double invcov[n_param][n_param], double param_diff[n_param]) 
{
  int c,r;
  double sum,prod[n_param];
 
  for (r = 0; r < n_param; r++) {
    sum = 0.0;  
    for (c = 0; c < n_param; c++) {
      sum = sum + invcov[r][c]*param_diff[c];
    }
    prod[r] = sum;
  }

  sum = 0.0;
  for (c = 0; c < n_param; c++) {
    sum = sum + param_diff[c]*prod[c];
  }

  return sum;
}

//CH
double log_L_Planck15_BAO_w0wa()
{
  double log_L = 0.;
  int n_param = 7;
  double param_fid[n_param], param_diff[n_param];
  double table[n_param][n_param]; 
  int c, r;

  table[0][0] = 3.44277e+05;
  table[0][1] = -3.28153e+03;
  table[0][2] = 4.22375e+04;
  table[0][3] = 5.74005e+04;
  table[0][4] = 1.40240e+04;
  table[0][5] = -2.22404e+06;
  table[0][6] = 2.33225e+05;
  table[1][1] = 6.75674e+03;
  table[1][2] = -5.86507e+03;
  table[1][3] = 1.83006e+03;
  table[1][4] = 4.61213e+02;
  table[1][5] = -6.38218e+03;
  table[1][6] = -3.54112e+03;
  table[2][2] = 1.12446e+05;
  table[2][3] = -1.11375e+04;
  table[2][4] = -2.73720e+03;
  table[2][5] = -6.67431e+04;
  table[2][6] = -5.45340e+03;
  table[3][3] = 1.65304e+04;
  table[3][4] = 3.98154e+03;
  table[3][5] = -5.04165e+05;
  table[3][6] = 4.44333e+04;
  table[4][4] = 9.64647e+02;
  table[4][5] = -1.22127e+05;
  table[4][6] = 1.06349e+04;
  table[5][5] = 2.06119e+07;
  table[5][6] = -1.01467e+06;
  table[6][6] = 2.70588e+05;

  for (c = 0; c < n_param; c++) {
    for (r = c + 1; r < n_param; r++) {
      table[r][c] = table[c][r];
    }
  }

  param_fid[0] = 3.5098925e-01;
  param_fid[1] = 8.0467525e-01;
  param_fid[2] = 9.6406051e-01;
  param_fid[3] = -5.0551771e-01;
  param_fid[4] = -1.46884482e+00;
  param_fid[5] = 5.462452e-02;
  param_fid[6] = 6.3983876e-01;

  param_diff[0] = cosmology.Omega_m-param_fid[0]; 
  param_diff[1] = cosmology.sigma_8-param_fid[1]; 
  param_diff[2] = cosmology.n_spec-param_fid[2]; 
  param_diff[3] = cosmology.w0-param_fid[3]; 
  param_diff[4] = cosmology.wa-param_fid[4];
  param_diff[5] = cosmology.omb-param_fid[5];
  param_diff[6] = cosmology.h0-param_fid[6];
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

 double log_L_ia()
{
  if (like.IA ==3){return 0.;}
  double log_L = 0.;
  log_L -=  pow(( nuisance.LF_alpha - prior.LF_alpha[0])/ prior.LF_alpha[1],2.0);
  log_L -=  pow(( nuisance.LF_P - prior.LF_P[0])/ prior.LF_P[1],2.0);
  log_L -=  pow(( nuisance.LF_Q - prior.LF_Q[0])/ prior.LF_Q[1],2.0);
  log_L -=  pow(( nuisance.LF_red_alpha - prior.LF_red_alpha[0])/ prior.LF_red_alpha[1],2.0);
  log_L -=  pow(( nuisance.LF_red_P - prior.LF_red_P[0])/ prior.LF_red_P[1],2.0);
  log_L -=  pow(( nuisance.LF_red_Q - prior.LF_red_Q[0])/ prior.LF_red_Q[1],2.0);
  return 0.5*log_L;
}

//CH 
double log_L_Planck15_BAO_H070p6_JLA_w0wa()
{
  double log_L = 0.;
  int n_param = 7;
  double param_fid[n_param], param_diff[n_param];
  double table[n_param][n_param]; 
  int c, r;
  table[0][0] = 1.25403e+06;
  table[0][1] = -1.68487e+04;
  table[0][2] = 1.07247e+05;
  table[0][3] = 1.46946e+05;
  table[0][4] = 3.72662e+04;
  table[0][5] = -3.15490e+06;
  table[0][6] = 1.22441e+06;
  table[1][1] = 6.40426e+03;
  table[1][2] = -7.19010e+03;
  table[1][3] = 2.79656e+02;
  table[1][4] = 8.82725e+01;
  table[1][5] = 1.14697e+04;
  table[1][6] = -1.85914e+04;
  table[2][2] = 1.14261e+05;
  table[2][3] = -4.71337e+03;
  table[2][4] = -1.19919e+03;
  table[2][5] = -1.15936e+05;
  table[2][6] = 7.10418e+04;
  table[3][3] = 2.42113e+04;
  table[3][4] = 6.07750e+03;
  table[3][5] = -5.72657e+05;
  table[3][6] = 1.38794e+05;
  table[4][4] = 1.53328e+03;
  table[4][5] = -1.44471e+05;
  table[4][6] = 3.49984e+04;
  table[5][5] = 2.60247e+07;
  table[5][6] = -1.26930e+06;
  table[6][6] = 1.44367e+06;

  for (c = 0; c < n_param; c++) {
    for (r = c + 1; r < n_param; r++) {
      table[r][c] = table[c][r];
    }
  }

  param_fid[0] = 3.08672e-01;
  param_fid[1] = 8.41330e-01;
  param_fid[2] = 9.63707e-01;
  param_fid[3] = -9.35997e-01;
  param_fid[4] = -3.82183e-01;
  param_fid[5] = 4.79953e-02;
  param_fid[6] = 6.80887e-01;
    
  param_diff[0] = cosmology.Omega_m-param_fid[0]; 
  param_diff[1] = cosmology.sigma_8-param_fid[1]; 
  param_diff[2] = cosmology.n_spec-param_fid[2]; 
  param_diff[3] = cosmology.w0-param_fid[3]; 
  param_diff[4] = cosmology.wa-param_fid[4];
  param_diff[5] = cosmology.omb-param_fid[5];
  param_diff[6] = cosmology.h0-param_fid[6];
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

//CH 2018/07/17
double log_L_Planck18_BAO_Riess18_Pantheon_w0wa()
{
  double log_L = 0.;
  int n_param = 7;
  double param_fid[n_param], param_diff[n_param];
  double table[n_param][n_param]; 
  int c, r;
  table[0][0] = 4.30741e+06;
  table[0][1] = -6.29653e+05;
  table[0][2] = 2.30239e+05;
  table[0][3] = 2.40526e+05;
  table[0][4] = 5.75393e+04;
  table[0][5] = -1.09309e+07;
  table[0][6] = 3.63645e+06;
  table[1][1] = 1.84962e+05;
  table[1][2] = -4.30980e+04;
  table[1][3] = 9.10127e+03;
  table[1][4] = 3.25867e+03;
  table[1][5] = 1.38575e+06;
  table[1][6] = -4.94260e+05;
  table[2][2] = 1.26977e+05;
  table[2][3] = -7.56824e+03;
  table[2][4] = -2.11814e+03;
  table[2][5] = -4.02713e+05;
  table[2][6] = 1.58351e+05;
  table[3][3] = 4.14100e+04;
  table[3][4] = 1.04261e+04;
  table[3][5] = -9.36525e+05;
  table[3][6] = 2.14609e+05;
  table[4][4] = 2.64176e+03;
  table[4][5] = -2.28031e+05;
  table[4][6] = 5.15943e+04;
  table[5][5] = 4.94180e+07;
  table[5][6] = -7.29198e+06;
  table[6][6] = 3.34690e+06;
  param_fid[0] = 2.97717e-01;
  param_fid[1] = 8.33421e-01;
  param_fid[2] = 9.64631e-01;
  param_fid[3] = -9.99383e-01;
  param_fid[4] = -2.88583e-01;
  param_fid[5] = 4.64837e-02;
  param_fid[6] = 6.93855e-01;

  for (c = 0; c < n_param; c++) {
    for (r = c + 1; r < n_param; r++) {
      table[r][c] = table[c][r];
    }
  }

  param_diff[0] = cosmology.Omega_m-param_fid[0]; 
  param_diff[1] = cosmology.sigma_8-param_fid[1]; 
  param_diff[2] = cosmology.n_spec-param_fid[2]; 
  param_diff[3] = cosmology.w0-param_fid[3]; 
  param_diff[4] = cosmology.wa-param_fid[4];
  param_diff[5] = cosmology.omb-param_fid[5];
  param_diff[6] = cosmology.h0-param_fid[6];
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

//CH 2018/07/17
double log_L_Planck18_BAO_w0wa()
{
  double log_L = 0.;
  int n_param = 7;
  double param_fid[n_param], param_diff[n_param];
  double table[n_param][n_param]; 
  int c, r;
  table[0][0] = 6.78890e+05;
  table[0][1] = -1.57721e+05;
  table[0][2] = 5.15938e+04;
  table[0][3] = 4.25312e+04;
  table[0][4] = 9.74958e+03;
  table[0][5] = -3.57564e+06;
  table[0][6] = 4.44898e+05;
  table[1][1] = 1.12864e+05;
  table[1][2] = -1.19553e+04;
  table[1][3] = 2.29815e+04;
  table[1][4] = 6.07805e+03;
  table[1][5] = 6.20572e+05;
  table[1][6] = -8.00370e+04;
  table[2][2] = 1.12060e+05;
  table[2][3] = -1.19553e+04;
  table[2][4] = -2.96474e+03;
  table[2][5] = -1.18012e+05;
  table[2][6] = -1.70556e+03;
  table[3][3] = 2.30273e+04;
  table[3][4] = 5.65183e+03;
  table[3][5] = -4.46117e+05;
  table[3][6] = 4.30261e+04;
  table[4][4] = 1.39546e+03;
  table[4][5] = -1.05193e+05;
  table[4][6] = 9.95900e+03;
  table[5][5] = 2.75140e+07;
  table[5][6] = -1.65229e+06;
  table[6][6] = 4.48558e+05;
  param_fid[0] = 3.41667e-01;
  param_fid[1] = 7.95474e-01;
  param_fid[2] = 9.64545e-01;
  param_fid[3] = -5.77797e-01;
  param_fid[4] = -1.31630e+00;
  param_fid[5] = 5.33105e-02;
  param_fid[6] = 6.48927e-01;

  for (c = 0; c < n_param; c++) {
    for (r = c + 1; r < n_param; r++) {
      table[r][c] = table[c][r];
    }
  }

  param_diff[0] = cosmology.Omega_m-param_fid[0]; 
  param_diff[1] = cosmology.sigma_8-param_fid[1]; 
  param_diff[2] = cosmology.n_spec-param_fid[2]; 
  param_diff[3] = cosmology.w0-param_fid[3]; 
  param_diff[4] = cosmology.wa-param_fid[4];
  param_diff[5] = cosmology.omb-param_fid[5];
  param_diff[6] = cosmology.h0-param_fid[6];
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

//CH 2018/07/17
double log_L_Planck18_w0()
{
  double log_L = 0.;
  int n_param = 6;
  double param_fid[n_param], param_diff[n_param];
  double table[n_param][n_param]; 
  int c, r;
  //var_derived: ['omegam', 'sigma8', 'ns', 'w', 'omegab', 'h0']
  table[0][0] = 2.96336e+06;
  table[0][1] = -3.29741e+05;
  table[0][2] = 1.05809e+05;
  table[0][3] = 2.71328e+05;
  table[0][4] = -9.32882e+06;
  table[0][5] = 2.53296e+06;
  table[1][1] = 1.70713e+05;
  table[1][2] = -3.67360e+04;
  table[1][3] = 4.38186e+04;
  table[1][4] = 6.95686e+05;
  table[1][5] = -2.09522e+05;
  table[2][2] = 1.20693e+05;
  table[2][3] = -2.22763e+04;
  table[2][4] = -7.31831e+04;
  table[2][5] = 4.37545e+04;
  table[3][3] = 7.39254e+04;
  table[3][4] = -1.27156e+06;
  table[3][5] = 2.69293e+05;
  table[4][4] = 4.79954e+07;
  table[4][5] = -6.50670e+06;
  table[5][5] = 2.42902e+06;
  param_fid[0] = 3.03785e-01;
  param_fid[1] = 8.22263e-01;
  param_fid[2] = 9.65449e-01;
  param_fid[3] = -1.04080e+00;
  param_fid[4] = 4.75908e-02;
  param_fid[5] = 6.86273e-01;

  for (c = 0; c < n_param; c++) {
    for (r = c + 1; r < n_param; r++) {
      table[r][c] = table[c][r];
    }
  }

  param_diff[0] = cosmology.Omega_m-param_fid[0]; 
  param_diff[1] = cosmology.sigma_8-param_fid[1]; 
  param_diff[2] = cosmology.n_spec-param_fid[2]; 
  param_diff[3] = cosmology.w0-param_fid[3]; 
  param_diff[4] = cosmology.omb-param_fid[4];
  param_diff[5] = cosmology.h0-param_fid[5];
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

double log_L_wlphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.bias_zphot_shear[i][1] == 0.){
      printf("external_prior.c: called log_L_wlphotoz while prior.bias_zphot_shear[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -=  pow((nuisance.bias_zphot_shear[i] - prior.bias_zphot_shear[i][0])/ prior.bias_zphot_shear[i][1],2.0);
  }
  if(redshift.shear_photoz!=4) {
    if (prior.sigma_zphot_shear[0][1] == 0.){
      printf("external_prior.c: called log_L_wlphotoz while prior.sigma_zphot_shear[0][1] not set.\nEXIT\n");
      exit(1);
    }
    log_L -=  pow((nuisance.sigma_zphot_shear[0] - prior.sigma_zphot_shear[0][0])/ prior.sigma_zphot_shear[0][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_clphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.clustering_Nbin; i++){
      if (prior.bias_zphot_clustering[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.bias_zphot_clustering[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }

   log_L -=  pow((nuisance.bias_zphot_clustering[i] - prior.bias_zphot_clustering[i][0])/ prior.bias_zphot_clustering[i][1],2.0);
  }
  if(redshift.clustering_photoz!=4) {
    if (prior.sigma_zphot_clustering[0][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.sigma_zphot_clustering[0][1] not set.\nEXIT\n");
      exit(1);
    }
    log_L -=  pow((nuisance.sigma_zphot_clustering[0] - prior.sigma_zphot_clustering[0][0])/ prior.sigma_zphot_clustering[0][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_shear_calib()
{
  int i;
  double log_L = 0.; 
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.shear_calibration_m[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.shear_calibration_m[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -=  pow((nuisance.shear_calibration_m[i] - prior.shear_calibration_m[i][0])/ prior.shear_calibration_m[i][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_clusterMobs()
{
  double log_L = 0.; 
  if (strcmp(Cluster.model, "default")==0){
    int i;
    log_L -=  pow((nuisance.cluster_Mobs_lgM0 - prior.cluster_Mobs_lgM0[0])/ prior.cluster_Mobs_lgM0[1],2.0);
    log_L -=  pow((nuisance.cluster_Mobs_alpha - prior.cluster_Mobs_alpha[0])/ prior.cluster_Mobs_alpha[1],2.0);
    log_L -=  pow((nuisance.cluster_Mobs_beta - prior.cluster_Mobs_beta[0])/ prior.cluster_Mobs_beta[1],2.0);
    log_L -=  pow((nuisance.cluster_Mobs_sigma - prior.cluster_Mobs_sigma[0])/ prior.cluster_Mobs_sigma[1],2.0);
    for (i=0;i<4; i++){
      log_L -=  pow((nuisance.cluster_completeness[i] - prior.cluster_completeness[0])/ prior.cluster_completeness[1],2.0);
    }
    log_L = 0.5*log_L;
  }
  return log_L;
}


double log_like_f_red(){ //use f_red in each tomography bin to constrain LF parameters, assuming 10% uncertainty in f_red measurement; evaluated at bin center to avoid cosmology dependence
  static double **table;
  int i;
  double f_red,n_all,chi_sqr = 0.;
  if (table ==0){
    double a;
    table = create_double_matrix(0, tomo.shear_Nbin-1, 0, 2);
    for (i = 0; i < tomo.shear_Nbin; i++){
      a = 1./(1.+0.5*(tomo.shear_zmin[i]+tomo.shear_zmax[i]));
      table[i][0] = a;
      table[i][1] = f_red_LF(survey.m_lim,a);
      table[i][2] = n_all_LF(survey.m_lim,a);
    }
  }
  for (i = 0; i < tomo.shear_Nbin; i++){
    f_red = f_red_LF(survey.m_lim, table[i][0]);
    n_all = n_all_LF(survey.m_lim, table[i][0]);
    chi_sqr += pow((f_red-table[i][1])/(0.1*table[i][1]),2.0); //assume 10% uncertainty in f_red determination
    chi_sqr += pow((n_all-table[i][2])/(0.05*table[i][2]),2.0); //assume 5% uncertainty in n_gal determination
  }
  return -0.5*chi_sqr;
}


