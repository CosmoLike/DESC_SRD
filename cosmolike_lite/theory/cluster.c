/* cluster number counts & cluster lensing routines*/

double lgM_obs(double N200, double a); //mass observable relation
double lgMmin(double a, int nN); //lower integration boundary for mass integrals over richness bin nN at scale factor a - change to log(limits.M_min) for general scatter distribution
double lgMmax(double a, int nN); //upper integration boundary for mass integrals over richness bin nN at scale factor a - change to log(limits.M_max) for general scatter distribution
double Mobs_x(double Nobs, double lgM, double a); // x in the covolution for scaling relation P(x)
double int_n_Mobs(double lgM, void* params);//number density integrand for cluster halo model calculations
double n_N200(int nz, int nN); //projected cluster density [1/radian]
double N_N200(int nz, int nN); // cluster numbercounts in redshift bin nz, richness bin nN
double dN200_dchi (double a, int nN);// d N_N200/dchi
double b_cluster(int nz, int nM); // cluster bias in redshift bin nz, richness bin nN
double W_cluster(double a, double nz,double nN); //projection weight for cluster power spectra - not including cluster bias!

double C_cl_g_tomo(double l, int nzc1, int nN1, int nzl1);//angular cluster x galaxy clustering power spectrum  - only needed in covariance computation
double C_ccl_tomo(double l, int nz1, int nN1, int nz2, int nN2);//angular cluster clustering power spectrum - only needed in covariance computation
double C_cgl_tomo_nointerp(double l, int nz, int nN, int zs);//angular cluster lensing power spectrum



/*********** N200-M relation routines ***********/
double lgM_obs(double N200, double a){
  return nuisance.cluster_Mobs_lgM0 + log(N200/nuisance.cluster_Mobs_N_pivot)*nuisance.cluster_Mobs_alpha + log(1./a)*nuisance.cluster_Mobs_beta;
}
double scatter_lgM_obs(double N200, double a){
  return nuisance.cluster_Mobs_sigma;
}

/*********** M->N200 relation routines with the scatter as a function os M and z ***********/
double lgN200_model(double M, double a){
  return nuisance.cluster_Mobs_lgN0 + log(M/(3.e+14))*nuisance.cluster_Mobs_alpha + log(1./a)*nuisance.cluster_Mobs_beta;
}
double scatter_lgN200_model_mz(double M, double a){
  return nuisance.cluster_Mobs_sigma0 + log(M/(3.e+14))*nuisance.cluster_Mobs_sigma_qm + log(1./a)*nuisance.cluster_Mobs_sigma_qz;
}

/********** systematics models *************/ //simple place-holder routines for now, to be 
double N200_completeness(int nz, int nN){
  return nuisance.cluster_completeness[nz];
}
double P_cm_offcentering(double k, double M, double a){
  double fc = nuisance.cluster_centering_f0 + nuisance.cluster_centering_alpha*log(M/nuisance.cluster_centering_M_pivot);
  if (nuisance.cluster_centering_f0 > 0.99 || fc >0.99){return 1.0;}
  return fc +(1.-fc)*exp(-0.5*pow(nuisance.cluster_centering_sigma*k/a,2.0));
}
/*********** cluster counts *******************/
double lgMmin(double a,int nN){ // speed up of mass integrals for Gaussian scatter - change to log(limits.M_min) for general scatter distribution
  if (strcmp(Cluster.model, "default")==0){
    return fmax(log(limits.M_min),lgM_obs(Cluster.N_min[nN],a)-5.*scatter_lgM_obs(Cluster.N_min[nN],a));
  } else if(strcmp(Cluster.model, "Murata_etal_2018")==0){ // foward model with mass- and redshift-dependent MOR scatter
    return log(1.e12);
  } else {
    fprintf(stderr, "Error: Cluster.model should be 'default' or 'Rykoff_2012'");
    exit(1);
  }
}
double lgMmax(double a,int nN){// speed up of mass integrals for Gaussian scatter - change to log(limits.M_max) for general scatter distribution
  if (strcmp(Cluster.model, "default")==0){
    return fmin(log(limits.M_max),lgM_obs(Cluster.N_max[nN],a)+5.*scatter_lgM_obs(Cluster.N_max[nN],a));
  } else if(strcmp(Cluster.model, "Murata_etal_2018")==0){
    return log(limits.M_max);
  } else {
    fprintf(stderr, "Error: Cluster.model should be 'default' or 'Rykoff_2012'");
    exit(1);
  }
}

double Mobs_x(double Nobs, double lgM, double a){
  if(strcmp(Cluster.model, "default")==0){ // Eq. B4 in Rykoff et al. (2012) http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1104.2089
    return (lgM_obs(Nobs,a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Nobs,a));
  } else if (strcmp(Cluster.model, "Murata_etal_2018")==0){ // Eq. 16 in Murata et al. (2018) http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1707.01907
    return (log(Nobs)-lgN200_model(exp(lgM),a))/(sqrt(2.0)*scatter_lgN200_model_mz(exp(lgM), a));
  } else {
    fprintf(stderr, "Error: Cluster.model should be 'default' or 'Rykoff_2012'");
    exit(1);
  }
}

double int_n_Mobs(double lgM, void* params){
  double *array = (double *) params;
  double x1,x2;
  x1 = Mobs_x(Cluster.N_min[(int)array[1]], lgM, array[0]);
  x2 = Mobs_x(Cluster.N_max[(int)array[1]], lgM, array[0]);
  return exp(lgM)*massfunc(exp(lgM),array[0])*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
}

double int_b_Mobs(double lgM, void* params){
  double *array = (double *) params;
  double x1,x2;
  x1 = Mobs_x(Cluster.N_min[(int)array[1]], lgM, array[0]);
  x2 = Mobs_x(Cluster.N_max[(int)array[1]], lgM, array[0]);
  return exp(lgM)*massfunc(exp(lgM),array[0])*B1(exp(lgM),array[0])*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
}

double int_M_Mobs(double lgM, void* params){
  double *array = (double *) params;
  double x1,x2;
  x1 = Mobs_x(Cluster.N_min[(int)array[1]], lgM, array[0]);
  x2 = Mobs_x(Cluster.N_max[(int)array[1]], lgM, array[0]);
  return exp(lgM)*exp(lgM)*massfunc(exp(lgM),array[0])*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
}

double int_da_n_Mobs(double a, void* params){
  double *array = (double *) params;
  array[0] = a;
  return pow(f_K(chi(a)),2.0)*dchi_da(a)*int_gsl_integrate_high_precision(int_n_Mobs, (void*)array,lgMmin(a,(int) array[1]),lgMmax(a,(int) array[1]),NULL,1000);
}

double n_N200 (int nz, int nN){
  static cosmopara C;
  static nuisancepara N;
  static double **table;
  if (recompute_clusters(C,N)){
    if (table==0) {table = create_double_matrix(0, tomo.cluster_Nbin-1, 0, Cluster.N200_Nbin-1);}
    double array[2];
    int z, n;
    for (z = 0; z < tomo.cluster_Nbin; z++){
      for (n = 0; n< Cluster.N200_Nbin; n++){
        array[1] = (double) n;
	table[z][n] = int_gsl_integrate_medium_precision(int_da_n_Mobs, (void*)array, 1./(tomo.cluster_zmax[z]+1), 1./(tomo.cluster_zmin[z]+1),NULL,1000);
      }
    }
    update_cosmopara(&C); update_nuisance(&N);
  }
  return table[nz][nN];
}

double dN200_dchi (double a, int nN){
  double array[2]= {a, (double) nN};
  return 4.*M_PI*survey.area/41253.0*pow(f_K(chi(a)),2.0)*int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array,lgMmin(a,(int) array[1]),lgMmax(a,(int) array[1]),NULL,1000);
}

double N_N200(int nz, int nN){
  return 4.*M_PI*survey.area/41253.0*n_N200(nz,nN)*N200_completeness(nz, nN);
}

/******** cluster n(z) and bias **********/

double b_cluster (int nz, int nN){ //mean bias of clusters in redshift bin nz, richness bin nN, evaluated at center of redshift bin
  static cosmopara C;
  static nuisancepara N;
  static double **table;
  if (recompute_clusters(C,N)){
    if (table==0) {table = create_double_matrix(0, tomo.cluster_Nbin-1, 0, Cluster.N200_Nbin-1);}
    double array[2],a;
    int z, n;
    for (z = 0; z < tomo.cluster_Nbin; z++){
      array[0] = 1./(1.+0.5*(tomo.cluster_zmin[z]+ tomo.cluster_zmax[z]));
      a= array[0];
      for (n = 0; n< Cluster.N200_Nbin; n++){
        array[1] = (double) n;
        table[z][n] = int_gsl_integrate_medium_precision(int_b_Mobs, (void*)array, lgMmin(a,n),lgMmax(a,n),NULL,1000)/int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array, lgMmin(a,n),lgMmax(a,n),NULL,1000);
      }
    }
    update_cosmopara(&C); update_nuisance(&N);
  }
  return table[nz][nN];
}
double dV_cluster(double z, void* params){ //placeholder routine, factor in selection function for actual analaysis
  return pow(f_K(chi(1./(1.+z))),2.0)/(hoverh0(1./(1+z)));
}
double norm_z_cluster(int nz){
  static cosmopara C;
  static double **table;
  if(recompute_expansion(C)){
    if (table==0) {table = create_double_matrix(0, tomo.cluster_Nbin-1, 0, 1);}
    int n;
    double array[2];
    for (n = 0; n < tomo.cluster_Nbin; n++){
      table[n][0] = int_gsl_integrate_high_precision(dV_cluster, (void*) array, tomo.cluster_zmin[n],tomo.cluster_zmax[n], NULL,1000);
    }
    update_cosmopara(&C);
  }
  return table[nz][0];
}
double zdistr_cluster(double z, int nz, int nN){ //simplfied selection function, disregards evolution of N-M relation + mass function within redshift bin
  if (z >tomo.cluster_zmax[nz]|| z< tomo.cluster_zmin[nz]){return 0.;}
  double a[2];
  return dV_cluster(z, &a)/norm_z_cluster(nz);
}


/************ cluster lensing routines ************/
double W_cluster(double a, double nz,double nN){
  return zdistr_cluster(1./a-1.,(int)nz,(int) nN)*hoverh0(a);
}

double int_P_cm_1h(double lgM, void *params){
  double x1,x2,m,a,u_m;
  double *array = (double *) params;
  m = exp(lgM); a = array[0];
  x1 = Mobs_x(Cluster.N_min[(int)array[1]], lgM, a);
  x2 = Mobs_x(Cluster.N_max[(int)array[1]], lgM, a);
  u_m = m/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(m,a),array[2],m,a); //density profile
  return m*massfunc(m,a)*u_m*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))*P_cm_offcentering(array[2],m,a);
}

double P_cm_1h(double k, double a,int nz, int nN){ 
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.,amin =0., amax = 0.;
  static int N_a = 20, N_k_nlin = 30;
  static int NN =-1;
  
  static double **table_P_cm=0;
  double klog,aa;
  int i,j;
  double array[3];
  if (nN != NN){
    NN = nN; 
    if (table_P_cm == 0){
      table_P_cm = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    }
    amax = 1./(1+tomo.cluster_zmin[0]);
    amin = 1./(1+tomo.cluster_zmax[tomo.cluster_Nbin-1]);
    da = (amax - amin)/(N_a-1.);
    logkmin = log(0.1*cosmology.coverH0);
    logkmax = log(100.*cosmology.coverH0);
    dk = (logkmax - logkmin)/(N_k_nlin-1.);
    aa= amin;
    array[1] = (double) nN;
    for (i=0; i<N_a; i++, aa +=da) {
      array[0] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        array[2]= exp(klog);
	table_P_cm[i][j] = log(int_gsl_integrate_medium_precision(int_P_cm_1h, (void*)array, lgMmin(aa,nN), lgMmax(aa,nN),NULL,1000)/int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array, lgMmin(aa,nN), lgMmax(aa,nN),NULL,1000));
      }
    }
  }
  if (k < 0.1*cosmology.coverH0){return 0.;}
  klog = log(k);
  aa = fmin(a,amax-1.1*da);
  return exp(interpol2d(table_P_cm, N_a, amin, amax, da, aa, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0));
}
double P_cm(double k,double a,int nz,int nN){
  return P_cm_1h(k,a,nz,nN)+b_cluster(nz,nN)*Pdelta(k,a);
}
double int_for_C_cgl_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  
  ell       = ar[3]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_cluster(a,ar[0],ar[1])*W_kappa(a,fK,ar[2])*dchi_da(a)/fK/fK;
  if (res !=0){res= res*P_cm(k,a, (int) ar[0], (int) ar[1]);}
  return res;
}
double int_for_C_ccl_tomo(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k;
  
  ell       = ar[3]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= b_cluster((int)(ar[0]),(int)(ar[1]))*b_cluster((int)(ar[0]),(int)(ar[2]))*W_cluster(a,ar[0],ar[1])*W_cluster(a,ar[0],ar[2])*dchi_da(a)/fK/fK;
  res= res*p_lin(k,a);
  return res;
}

double int_for_C_cl_g_tomo(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k;
  
  ell       = ar[3]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= b_cluster((int)(ar[0]),(int)(ar[1]))*W_cluster(a,ar[0],ar[1])*W_gal(a,ar[2])*dchi_da(a)/fK/fK;
  res= res*p_lin(k,a);
  return res;
}

double C_cl_g_tomo(double l, int nzc1, int nN1, int nzl1){
  double array[4] ={(double)nzc1, (double)nN1,  (double)nzl1, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_g_tomo, (void*) array, 1./(1+fmin(tomo.cluster_zmax[nzc1],tomo.clustering_zmax[nzl1])),1./(1+fmax(tomo.cluster_zmin[nzc1],tomo.clustering_zmin[nzl1])), NULL,1000);
}
double C_ccl_tomo(double l, int nz1, int nN1, int nz2, int nN2){
  if (nz1 != nz2){return 0.;}
  double array[4] ={(double)nz1, (double)nN1,  (double)nN2, l};
  return int_gsl_integrate_low_precision(int_for_C_ccl_tomo, (void*) array, 1./(1+tomo.cluster_zmax[nz1]),1./(1+tomo.cluster_zmin[nz1]), NULL,1000);
}
double C_cgl_tomo_nointerp(double l, int nz, int nN, int zs){
  double array[4] = {(double)nz,(double) nN, (double) zs,l};
  return int_gsl_integrate_low_precision(int_for_C_cgl_tomo, (void*) array, 1./(1+tomo.cluster_zmax[nz]),1./(1+tomo.cluster_zmin[nz]), NULL,1000);
}
