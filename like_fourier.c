#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>


#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "init_SRD.c"

double C_shear_tomo_sys(double ell,int z1,int z2);
double C_cgl_tomo_sys(double ell_Cluster,int zl,int nN, int zs);
double C_gl_tomo_sys(double ell,int zl,int zs);
void set_data_shear(int Ncl, double *ell, double *data, int start);
void set_data_ggl(int Ncl, double *ell, double *data, int start);
void set_data_clustering(int Ncl, double *ell, double *data, int start);
void set_data_cluster_N(double *data, int start);
void set_data_cgl(double *ell_Cluster, double *data, int start);
void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q, double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope);
double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q, double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope);
double write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in);
double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in);
int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_ell(void);

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}
int get_N_ell(void){
  return like.Ncl;
}


double C_shear_tomo_sys(double ell, int z1, int z2)
{
  double C;
  // C= C_shear_tomo_nointerp(ell,z1,z2);
  // if(like.IA==1) C+=C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  
  if(like.IA!=1) C= C_shear_tomo_nointerp(ell,z1,z2);
  if(like.IA==1) C= C_shear_shear_IA(ell,z1,z2);
  //if(like.IA==1) C = C_shear_tomo_nointerp(ell,z1,z2)+C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  if(like.IA==2) C += C_II_lin_nointerp(ell,z1,z2)+C_GI_lin_nointerp(ell,z1,z2);  
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
  //printf("%le %d %d %le\n",ell,z1,z2,C_shear_tomo_nointerp(ell,z1,z2)+C_II_JB_nointerp(ell,z1,z2)+C_GI_JB_nointerp(ell,z1,z2));
return C;
}

double C_gl_tomo_sys(double ell,int zl,int zs)
{
  double C;
  // C=C_gl_tomo_nointerp(ell,zl,zs); 
  // if(like.IA==1) C += C_gI_nointerp(ell,zl,zs);
  
  if(like.IA!=1) C=C_gl_tomo_nointerp(ell,zl,zs);
  if(like.IA==1) C = C_ggl_IA(ell,zl,zs);
  if(like.IA==2) C += C_gI_lin_nointerp(ell,zl,zs);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}

double C_cgl_tomo_sys(double ell_Cluster, int zl,int nN, int zs)
{
  double C;
  C=C_cgl_tomo_nointerp(ell_Cluster,zl,nN,zs);
  //if(like.IA!=0) C += 
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}      

void set_data_shear(int Ncl, double *ell, double *data, int start)
{
  int i,z1,z2,nz;
  double a;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < Ncl; i++){
      if (ell[i]<like.lmax_shear){ data[Ncl*nz+i] = C_shear_tomo_sys(ell[i],z1,z2);}
      else {data[Ncl*nz+i] = 0.;}
    }
  }
}

void set_data_ggl(int Ncl, double *ell, double *data, int start)
{
  int i, zl,zs,nz;  
  for (nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
    zl = ZL(nz); zs = ZS(nz);
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],zl)){
        data[start+(Ncl*nz)+i] = C_gl_tomo_sys(ell[i],zl,zs);
      }
      else{
        data[start+(Ncl*nz)+i] = 0.;
      }
    } 
  }
}

void set_data_clustering(int Ncl, double *ell, double *data, int start){
  int i, nz;
  for (nz = 0; nz < tomo.clustering_Npowerspectra; nz++){
    //printf("%d %e %e\n",nz, gbias.b[nz][1],pf_photoz(gbias.b[nz][1],nz));
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],nz)){data[start+(Ncl*nz)+i] = C_cl_tomo_nointerp(ell[i],nz,nz);}
      else{data[start+(Ncl*nz)+i] = 0.;}
      //printf("%d %d %le %le\n",nz,nz,ell[i],data[Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + nz)+i]);
    }
  }
}

void set_data_cluster_N(double *data, int start){
  int nN, nz;
  for (nz = 0; nz < tomo.cluster_Nbin; nz++){
    for (nN = 0; nN < Cluster.N200_Nbin; nN++){
      data[start+Cluster.N200_Nbin*nz+nN] = N_N200(nz, nN);
    }
  }
}


void set_data_cgl(double *ell_Cluster, double *data, int start)
{
  int zl,zs,nN,nz,i,j;
  for(nN = 0; nN < Cluster.N200_Nbin; nN++){
    for (nz = 0; nz < tomo.cgl_Npowerspectra; nz++){
      zl = ZC(nz); zs = ZSC(nz);
      for (i = 0; i < Cluster.lbin; i++){
        j = start;
        j += (nz*Cluster.N200_Nbin+nN)*Cluster.lbin +i;
        data[j] = C_cgl_tomo_sys(ell_Cluster[i],zl,nN,zs);
      }
    }
  }
}


int set_cosmology_params(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=S8;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  cosmology.h0=H0;
  cosmology.MGSigma=MGSigma;
  cosmology.MGmu=MGmu;

  if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
  if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  return 1;
}

void set_nuisance_shear_calib(double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10)
{
  nuisance.shear_calibration_m[0] = M1;
  nuisance.shear_calibration_m[1] = M2;
  nuisance.shear_calibration_m[2] = M3;
  nuisance.shear_calibration_m[3] = M4;
  nuisance.shear_calibration_m[4] = M5;
  nuisance.shear_calibration_m[5] = M6;
  nuisance.shear_calibration_m[6] = M7;
  nuisance.shear_calibration_m[7] = M8;
  nuisance.shear_calibration_m[8] = M9;
  nuisance.shear_calibration_m[9] = M10;
}

int set_nuisance_shear_photoz(double SP1,double SP2,double SP3,double SP4,double SP5,double SP6,double SP7,double SP8,double SP9,double SP10,double SPS1)
{
  int i;
  nuisance.bias_zphot_shear[0]=SP1;
  nuisance.bias_zphot_shear[1]=SP2;
  nuisance.bias_zphot_shear[2]=SP3;
  nuisance.bias_zphot_shear[3]=SP4;
  nuisance.bias_zphot_shear[4]=SP5;
  nuisance.bias_zphot_shear[5]=SP6;
  nuisance.bias_zphot_shear[6]=SP7;
  nuisance.bias_zphot_shear[7]=SP8;
  nuisance.bias_zphot_shear[8]=SP9;
  nuisance.bias_zphot_shear[9]=SP10;
  
  for (i=0;i<tomo.shear_Nbin; i++){ 
    nuisance.sigma_zphot_shear[i]=SPS1;
    if (nuisance.sigma_zphot_shear[i]<0.001) return 0;
  }
  return 1;
}

int set_nuisance_clustering_photoz(double CP1,double CP2,double CP3,double CP4,double CP5,double CP6,double CP7,double CP8,double CP9,double CP10,double CPS1)
{
  int i;
  nuisance.bias_zphot_clustering[0]=CP1;
  nuisance.bias_zphot_clustering[1]=CP2;
  nuisance.bias_zphot_clustering[2]=CP3;
  nuisance.bias_zphot_clustering[3]=CP4;
  nuisance.bias_zphot_clustering[4]=CP5;
  nuisance.bias_zphot_clustering[5]=CP6;
  nuisance.bias_zphot_clustering[6]=CP7;
  nuisance.bias_zphot_clustering[7]=CP8;
  nuisance.bias_zphot_clustering[8]=CP9;
  nuisance.bias_zphot_clustering[9]=CP10;
  
  for (i=0;i<tomo.clustering_Nbin; i++){ 
    nuisance.sigma_zphot_clustering[i]=CPS1;
    if (nuisance.sigma_zphot_clustering[i]<0.001) return 0;
  }
  return 1;
}

int set_nuisance_ia(double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q)
{

  nuisance.A_ia=A_ia;  
  nuisance.beta_ia=beta_ia;
  nuisance.eta_ia=eta_ia;
  nuisance.eta_ia_highz=eta_ia_highz;
  nuisance.LF_alpha=LF_alpha;
  nuisance.LF_P=LF_P;
  nuisance.LF_Q=LF_Q;
  nuisance.LF_red_alpha=LF_red_alpha;
  nuisance.LF_red_P=LF_red_P;
  nuisance.LF_red_Q=LF_red_Q;

  if (nuisance.A_ia < 0.0 || nuisance.A_ia > 10.0) return 0;
  if (nuisance.beta_ia < -1.0 || nuisance.beta_ia > 3.0) return 0;
  if (nuisance.eta_ia < -3.0 || nuisance.eta_ia> 3.0) return 0;
  if (nuisance.eta_ia_highz < -1.0 || nuisance.eta_ia_highz> 1.0) return 0;
  // printf("%le %le %le %le %le %le %le %le %le %le\n",nuisance.A_ia,nuisance.beta_ia,nuisance.eta_ia,nuisance.eta_ia_highz,nuisance.LF_alpha,nuisance.LF_P, nuisance.LF_Q,nuisance.LF_red_alpha,nuisance.LF_red_P,nuisance.LF_red_Q);
  // if(like.IA!=0){
  //  if (check_LF()) return 0;
  // }
return 1;
}

int set_nuisance_cluster_Mobs(double cluster_Mobs_lgN0,  double cluster_Mobs_alpha, double cluster_Mobs_beta, double cluster_Mobs_sigma0, double cluster_Mobs_sigma_qm, double cluster_Mobs_sigma_qz)
{
  //  nuisance.cluster_Mobs_lgM0 = mass_obs_norm;  //fiducial : 1.72+log(1.e+14*0.7); could use e.g. sigma = 0.2 Gaussian prior
  //  nuisance.cluster_Mobs_alpha = mass_obs_slope; //fiducial: 1.08; e.g. sigma = 0.1 Gaussian prior
  //  nuisance.cluster_Mobs_beta = mass_z_slope; //fiducial: 0.0; e.g. sigma = 0.1 Gaussian prior
  //  nuisance.cluster_Mobs_sigma = mass_obs_scatter; //fiducial 0.25; e.g. sigma = 0.05 Gaussian prior

  // fiducial values and priors from Murata et al. (2018) except for redshift-related parameters
  nuisance.cluster_Mobs_lgN0 = cluster_Mobs_lgN0; //fiducial: 3.207, flat prior [0.5, 5.0]
  nuisance.cluster_Mobs_alpha = cluster_Mobs_alpha; //fiducial: 0.993, flat prior [0.0, 2.0]
  nuisance.cluster_Mobs_beta = cluster_Mobs_beta; //fiducial: 0.0, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma0 = cluster_Mobs_sigma0; //fiducial: 0.456, flat prior [0.0, 1.5]
  nuisance.cluster_Mobs_sigma_qm = cluster_Mobs_sigma_qm; //fiducial: -0.169, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma_qz = cluster_Mobs_sigma_qz; //fiducial: 0.0, flat prior [-1.5, 1.5]

  if (nuisance.cluster_Mobs_lgN0 < 0.5 || nuisance.cluster_Mobs_lgN0 > 5.0) return 0;
  if (nuisance.cluster_Mobs_alpha < 0.0 || nuisance.cluster_Mobs_alpha > 2.0) return 0;
  if (nuisance.cluster_Mobs_beta < -1.5 || nuisance.cluster_Mobs_beta > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma0 < 0.0|| nuisance.cluster_Mobs_sigma0 > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma_qm < -1.5 && nuisance.cluster_Mobs_sigma_qm > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma_qz < -1.5 && nuisance.cluster_Mobs_sigma_qz > 1.5)return 0;

return 1;
}

int set_nuisance_gbias(double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8,double B9, double B10)
{
  int i;
  gbias.b[0] = B1;
  gbias.b[1] = B2;
  gbias.b[2] = B3;
  gbias.b[3] = B4;
  gbias.b[4] = B5;
  gbias.b[5] = B6;
  gbias.b[6] = B7;
  gbias.b[7] = B8;
  gbias.b[8] = B9;
  gbias.b[9] = B10;
  if(like.bias==1){
    for (i = 0; i < 10; i++){
      if (gbias.b[i] < 0.8 || gbias.b[i] > 3.0) return 0;
    }
  }
  return 1;
} 



double log_L_SRD_SN_Y1_RENEE() // computed from /Users/teifler/Dropbox/cosmolike_store/LSSTawakens/ReneeSNChains/lsst_y1_jan29_nostarts_varyM_oldomb_rand.txt
{
  double log_L = 0.;
  double param_diff[4];
  //double inv[4][4]={{3.79965832e+03,-1.67422234e+01, 1.49107053e+03, 2.03605752e+02},{ -1.67422234e+01, 4.73633438e+02,-2.97547647e+00,-2.21008677e-01},{1.49107053e+03,-2.97547647e+00, 7.06767724e+02, 7.78147349e+01},{2.03605752e+02,-2.21008677e-01, 7.78147349e+01, 1.20866157e+01}};

  double inv[4][4]={{2.54176679e+03,2.52990975e+01,8.90935426e+02,1.36402504e+02},{2.52990975e+01,4.81470268e+02,1.01904094e+01,6.95186662e-01},{8.90935426e+02,1.01904094e+01,4.02864434e+02,4.60275519e+01},{1.36402504e+02,6.95186662e-01,4.60275519e+01,8.41721305e+00}};
  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.h0-0.6727; 
  param_diff[2] = cosmology.w0+1.0;  
  param_diff[3] = cosmology.wa;
 
  log_L = -0.5*do_matrix_mult_invcov(4,inv,param_diff);
  return log_L;
 }

double log_L_SRD_SN_Y10_RENEE() // computed from /Users/teifler/Dropbox/cosmolike_store/LSSTawakens/ReneeSNChains/lsst_y10_jan29_nostarts_varyM_oldomb_rand.txt
{
  double log_L = 0.;
  double param_diff[4];
  //double inv[4][4]={{3.11715166e+04,-0-6.05075241e+01, 1.86310837e+04, 1.66544964e+03},{-6.05075241e+01, 5.05572825e+02, -4.08007880e-01,-4.03507564e+00}, {1.86310837e+04,-4.08007880e-01, 1.25651434e+04, 8.67021009e+02},{1.66544964e+03,-4.03507564e+00, 8.67021009e+02,1.06630810e+02}};
  double inv[4][4]={{1.16609820e+04, 5.07477631e+01, 5.54091287e+03, 6.60416966e+02},{5.07477631e+01, 4.82820782e+02,-1.02706789e+01, 3.04465946e+00},{5.54091287e+03,-1.02706789e+01, 3.29658857e+03, 2.95282164e+02},{6.60416966e+02, 3.04465946e+00, 2.95282164e+02, 3.99928278e+01}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.h0-0.6727; 
  param_diff[2] = cosmology.w0+1.0; 
  param_diff[3] = cosmology.wa; 

  log_L = -0.5*do_matrix_mult_invcov(4,inv, param_diff);
  return log_L;
 }




double log_L_SRD_SL_Y1_TOM() // computed from /Users/teifler/Dropbox/cosmolike_store/LSSTawakens/TomSLChains/SL_LSSTY1.txt
{
  double log_L = 0.;
  double param_diff[4];
  double inv[4][4]={{9.57953440e+01,1.90114212e+01,1.13561665e+00,2.33866782e+02},{1.90114212e+01,1.31631828e+01,-1.88534652e-02,1.45671419e+02},{1.13561665e+00,-1.88534652e-02,6.29296427e-01,-6.10433226e+00},{2.33866782e+02,1.45671419e+02,-6.10433226e+00,2.71769694e+03}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.w0+1.0; 
  param_diff[2] = cosmology.wa; 
  param_diff[3] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(4,inv, param_diff);
  return log_L;
 }

double log_L_SRD_SL_Y10_TOM() // computed from /Users/teifler/Dropbox/cosmolike_store/LSSTawakens/TomSLChains/SL_LSSTY10.txt
{
  double log_L = 0.;
  double param_diff[4];
  double inv[4][4]={{2.44407663e+02,1.31525702e+02,1.62494864e+01,1.03889766e+03},{1.31525702e+02,1.90748959e+02,8.03881315e+00,1.99392208e+03},{1.62494864e+01,8.03881315e+00,2.23917868e+00,-6.64709076e+00},{1.03889766e+03,1.99392208e+03,-6.64709076e+00,3.31389425e+04}};


  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.w0+1.0; 
  param_diff[2] = cosmology.wa; 
  param_diff[3] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(4,inv, param_diff);
  return log_L;
 }


double log_L_SRD_LSST_Y1_clusters() // // including MOR for Y1, see fisher.py
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{3.858101e+04, 3.588867e+04, 2.231543e+03, -9.737292e+02, -2.076847e+02, -1.567836e+03, 5.165868e+02},{3.588867e+04, 4.046241e+04, 2.204910e+03, -4.306084e+02, -1.921012e+02, -1.153062e+03, 4.146891e+02},{2.231543e+03, 2.204910e+03, 4.640472e+02, -3.849041e+01, -1.164098e+01, -2.997845e+02, 8.606974e+01},{-9.737292e+02, -4.306084e+02, -3.849041e+01, 1.057829e+02, 1.234460e+01, 3.159550e+01, -1.021792e+01},{-2.076847e+02, -1.921012e+02, -1.164098e+01, 1.234460e+01, 3.965104e+00, 6.461167e+00, -2.324549e+00},{-1.567836e+03, -1.153062e+03, -2.997845e+02, 3.159550e+01, 6.461167e+00, 9.845550e+04, -2.104350e+02},{5.165868e+02, 4.146891e+02, 8.606974e+01, -1.021792e+01, -2.324549e+00, -2.104350e+02, 2.121782e+02}};


  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }

double log_L_SRD_LSST_Y10_clusters() // including MOR for Y10, see fisher.py
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{8.011883e+04, 7.501415e+04, 4.972257e+03, -1.009514e+03, -2.332928e+02, -4.035026e+03, 1.280726e+03},{7.501415e+04, 9.212235e+04, 5.029562e+03, 1.256229e+03, -8.904458e+00, -3.884022e+03, 1.236633e+03},{4.972257e+03, 5.029562e+03, 8.996244e+02, -1.335280e+01, -1.271912e+01, -9.151661e+02, 2.624954e+02},{-1.009514e+03, 1.256229e+03, -1.335280e+01, 3.943470e+02, 5.659304e+01, -3.776163e+01, 7.168167e+00},{-2.332928e+02, -8.904458e+00, -1.271912e+01, 5.659304e+01, 1.167193e+01, 4.161086e+00, -1.872393e+00},{-4.035026e+03, -3.884022e+03, -9.151661e+02, -3.776163e+01, 4.161086e+00, 1.000908e+05, -6.424232e+02},{1.280726e+03, 1.236633e+03, 2.624954e+02, 7.168167e+00, -1.872393e+00, -6.424232e+02, 3.273443e+02}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }

double log_L_SRD_stage3_exceptw0wa() 
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{8.52388933e+03,-2.38254882e-03,4.32405053e+00,0.00000000e+00,0.00000000e+00,1.22554260e+02,5.48470403e+01},{-2.38254882e-03,3.34465215e+01,3.86664237e-01,0.00000000e+00,0.00000000e+00,-8.71498515e+01,-7.09315084e+00},{4.32405053e+00,3.86664237e-01,2.47843770e+02,0.00000000e+00,0.00000000e+00,-1.21079947e+02,4.07696833e+00},{0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},{0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},{1.22554260e+02,-8.71498515e+01,-1.21079947e+02,-0.00000000e+00,-0.00000000e+00,7.17126110e+04,3.07891510e+03},{5.48470403e+01,-7.09315084e+00,4.07696833e+00,0.00000000e+00,0.00000000e+00,3.07891510e+03,5.31816639e+02}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }

double log_L_SRD_stage3() 
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{8.52889541e+03,9.34266060e-03,2.21513407e+00,1.85008170e+01,6.12185662e+00,1.17720010e+02,5.30670581e+01},{9.34266060e-03,3.34475078e+01,3.79135616e-01,-8.02192934e-02,5.88736315e-02,-8.69202865e+01,-7.08029577e+00},{2.21513407e+00,3.79135616e-01,2.48739187e+02,-7.46022981e+00,-2.69922217e+00,-1.19693897e+02,4.78085262e+00},{1.85008170e+01,-8.02192934e-02,-7.46022981e+00,8.42937127e+01,1.68856256e+01,-4.89063661e+01,-8.77194357e+00},{6.12185662e+00,5.88736315e-02,-2.69922217e+00,1.68856256e+01,9.55489400e+00,5.27704214e+00,-1.38597499e+00},{1.17720010e+02,-8.69202865e+01,-1.19693897e+02,-4.89063661e+01,5.27704214e+00,7.17777988e+04,3.08491105e+03},{5.30670581e+01,-7.08029577e+00,4.78085262e+00,-8.77194357e+00,-1.38597499e+00,3.08491105e+03,5.32751808e+02}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }


double log_L_SRD_LSST_Y1_3x2pt() 
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{6.461019e+05, 4.105098e+05, 9.063469e+04, -2.968463e+04, -6.656297e+03, -3.266829e+05, 7.947019e+04},{4.105098e+05, 3.102396e+05, 4.061214e+04, -2.051833e+04, -5.069455e+03, -1.156574e+05, 2.911602e+04},{9.063469e+04, 4.061214e+04, 3.652491e+04, -1.033492e+03, -1.571074e+02, -9.436089e+04, 2.396154e+04},{-2.968463e+04, -2.051833e+04, -1.033492e+03, 2.081364e+03, 4.721662e+02, 7.414854e+03, -1.864708e+03},{-6.656297e+03, -5.069455e+03, -1.571074e+02, 4.721662e+02, 1.200287e+02, 1.340576e+03, -3.293614e+02},{-3.266829e+05, -1.156574e+05, -9.436089e+04, 7.414854e+03, 1.340576e+03, 5.388793e+05, -9.884641e+04},{7.947019e+04, 2.911602e+04, 2.396154e+04, -1.864708e+03, -3.293614e+02, -9.884641e+04, 2.373400e+04}};


  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }


double log_L_SRD_LSST_Y10_3x2pt() 
{
  double log_L = 0.;
  double param_diff[7];
  double inv[7][7]={{1.732125e+06, 1.363116e+06, 2.420885e+05, -6.384395e+04, -1.991510e+04, -7.567872e+05, 1.858361e+05},{1.363116e+06, 1.196537e+06, 1.469326e+05, -5.148233e+04, -1.796618e+04, -4.526364e+05, 1.126267e+05},{2.420885e+05, 1.469326e+05, 1.011204e+05, -1.370858e+03, -4.391526e+02, -2.001808e+05, 5.228884e+04},{-6.384395e+04, -5.148233e+04, -1.370858e+03, 4.153843e+03, 1.137490e+03, 1.286287e+04, -3.253747e+03},{-1.991510e+04, -1.796618e+04, -4.391526e+02, 1.137490e+03, 3.720224e+02, 4.198170e+03, -1.020253e+03},{-7.567872e+05, -4.526364e+05, -2.001808e+05, 1.286287e+04, 4.198170e+03, 8.507285e+05, -1.685201e+05},{1.858361e+05, 1.126267e+05, 5.228884e+04, -3.253747e+03, -1.020253e+03, -1.685201e+05, 4.064231e+04}};

  param_diff[0] = cosmology.Omega_m-0.3156; 
  param_diff[1] = cosmology.sigma_8-0.831;
  param_diff[2] = cosmology.n_spec-0.9645; 
  param_diff[3] = cosmology.w0+1.0; 
  param_diff[4] = cosmology.wa; 
  param_diff[5] = cosmology.omb-0.0491685; 
  param_diff[6] = cosmology.h0-0.6727; 
  
  log_L = -0.5*do_matrix_mult_invcov(7,inv, param_diff);
  return log_L;
 }

double log_L_DESI() 
{
  double log_L = 0.;
  double param_diff[2];
  double inv[2][2]={{1411.40888852,449.46652105},{449.46652105,   291.94320759}};
  
  param_diff[0] = cosmology.w0+1.0; 
  param_diff[1] = cosmology.wa; 
  
  log_L = -0.5*do_matrix_mult_invcov(2,inv, param_diff);
  return log_L;
 }


double scatter(double N200, double z)
{
  double lnM, M, a;
  lnM = 1./nuisance.cluster_Mobs_alpha*(log(N200)-nuisance.cluster_Mobs_lgN0-nuisance.cluster_Mobs_beta*log(1.+z));
  M = exp(lnM)*3.e14;
  a = 1./(1.+z);
  return scatter_lgN200_model_mz(M, a);
}

double scatter_fid(double N200, double z)
{
  double cluster_Mobs_lgN0 = 3.207;
  double cluster_Mobs_alpha = 0.993;
  double cluster_Mobs_beta = 0.0;
  double cluster_Mobs_sigma0 = 0.456;
  double cluster_Mobs_sigma_qm = 0.0; // SRD suggestion from Eduardo
//  double cluster_Mobs_sigma_qm = -0.169;
  double cluster_Mobs_sigma_qz = 0.0;

  double lnM, M, a;
  lnM = 1./cluster_Mobs_alpha*(log(N200)-cluster_Mobs_lgN0-cluster_Mobs_beta*log(1.+z));
  M = exp(lnM)*3.e14;
  a = 1./(1.+z);
  return cluster_Mobs_sigma0 + log(M/(3.e+14))*cluster_Mobs_sigma_qm + log(1./a)*cluster_Mobs_sigma_qz;
}

double log_L_SRD_LSST_CL_EDUARDO() // values for scatter function are "arbitrary" but based on some experience (Hironao/Eduarod)
{
  double log_L = 0.;
  double s_100_0p2 = scatter(100., 0.2);
  double s_30_0p2 = scatter(30., 0.2);
  double s_100_0p8 = scatter(100., 0.8); // remove if redshift dependent scatter is excluded in analysis
  double f_100_0p2 = scatter_fid(100., 0.2);
  double f_30_0p2 = scatter_fid(30., 0.2);
  double f_100_0p8 = scatter_fid(100., 0.8); // remove if redshift dependent scatter is excluded in analysis
  double sig = 0.1; // 0.05 weighting the giants 

  log_L = -0.5*(exp(pow((s_100_0p2-f_100_0p2)/sig, 2.0))+ exp(pow((s_30_0p2-f_30_0p2)/sig, 2.0)) + exp(pow((s_100_0p8-f_100_0p8)/sig, 2.0)));

  return log_L;
}


double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q, double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope)
{
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0, log_L_GRS=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
    ell_Cluster= create_double_vector(0, Cluster.lbin-1);
    darg=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
    for (l=0;l<Cluster.lbin;l++){
      ell_Cluster[l]=exp(log(Cluster.l_min)+(l+0.5)*darg);
    }
  }
  if (set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e8;
  }
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  if (set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e8; 
  }
  if (set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10)==0){
    printf("Bias out of bounds\n");
    return -1.0e8;
  }
  if (set_nuisance_cluster_Mobs(mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm, mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope)==0){
    printf("Mobs out of bounds\n");
    return -1.0e8;
  }
       
  //printf("like %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  // if(like.Aubourg_Planck_BAO_SN==1) log_L_prior+=log_L_Planck_BAO_SN();
  // if(like.SN==1) log_L_prior+=log_L_SN();
  // if(like.BAO==1) log_L_prior+=log_L_BAO();
  // if(like.Planck==1) log_L_prior+=log_L_Planck();
  // if(like.IA!=0) log_L_prior+=log_L_ia();
  // if(like.IA!=0) log_L_prior+=log_like_f_red();
  // if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  // if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  // if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();
  // if(like.clusterMobs==1) log_L_prior+=log_L_clusterMobs();
  
  // if(like.SRD_SN_Y10==1) log_L_prior+=log_L_SRD_SN_Y10_RENEE();
  // if(like.SRD_SN_Y1==1) log_L_prior+=log_L_SRD_SN_Y1_RENEE();
  // if(like.SRD_SL_Y10==1) log_L_prior+=log_L_SRD_SL_Y10_TOM();
  // if(like.SRD_SL_Y1==1) log_L_prior+=log_L_SRD_SL_Y1_TOM();
  
  // log_L_prior+=log_L_DESI();

   log_L_prior+=log_L_SRD_SN_Y10_RENEE();
   log_L_prior+=log_L_SRD_SL_Y10_TOM();
  log_L_prior+=log_L_SRD_LSST_Y10_3x2pt(); 
  log_L_prior+=log_L_SRD_LSST_Y10_clusters(); 
  log_L_prior+=log_L_SRD_stage3(); 
  // log_L_prior+=log_L_SRD_SN_Y1_RENEE();
  // log_L_prior+=log_L_SRD_SL_Y1_TOM();
  // log_L_prior+=log_L_SRD_LSST_Y1_3x2pt();
  // log_L_prior+=log_L_SRD_LSST_Y1_clusters();
  // log_L_SRD_LSST_CL_EDUARDO();
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  // int start=0;  
  
  // if(like.shear_shear==1) {
  //   set_data_shear(like.Ncl, ell, pred, start);
  //   start=start+like.Ncl*tomo.shear_Npowerspectra;
  // }
  // if(like.shear_pos==1){
  //   set_data_ggl(like.Ncl, ell, pred, start);
  //   start=start+like.Ncl*tomo.ggl_Npowerspectra;
  // } 
  // if(like.pos_pos==1){
  //   set_data_clustering(like.Ncl,ell,pred, start);
  //   start=start+like.Ncl*tomo.clustering_Npowerspectra;
  // }
  // if(like.clusterN==1){ 
  //   set_data_cluster_N(pred,start);
  //   start=start+tomo.cluster_Nbin*Cluster.N200_Nbin;
  // }
  // if(like.clusterWL==1){
  //   set_data_cgl(ell_Cluster,pred, start);
  // }
  chisqr=0.0;
  // for (i=0; i<like.Ndata; i++){
  //   for (j=0; j<like.Ndata; j++){
  //     a=(pred[i]-data_read(1,i))*invcov_read(1,i,j)*(pred[j]-data_read(1,j));
  //     chisqr=chisqr+a;
  //   }
  //  //  if (fabs(data_read(1,i)/pred[i]-1.0) >1.e-4){
  //  // printf("%d %le %le %le\n",i,data_read(1,i),pred[i],data_read(1,i)/pred[i]);
  //  //  }
  // }
  // if (chisqr<0.0){
  //   printf("errror: chisqr < 0\n");
  // }
  // if (chisqr<-1.0) exit(EXIT_FAILURE);
  
  return -0.5*chisqr+log_L_prior;
}


void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q, double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope)
{
    
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
    ell_Cluster= create_double_vector(0, Cluster.lbin-1);
    darg=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
    for (l=0;l<Cluster.lbin;l++){
      ell_Cluster[l]=exp(log(Cluster.l_min)+(l+0.5)*darg);    
    }
  }
// for (l=0;l<like.Ncl;l++){
//   printf("%d %le\n",i,ell[l]);
// }
  set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu);
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1);
  set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1);
  set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q);
  set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10);
  set_nuisance_cluster_Mobs(mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm, mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope);
  
  int start=0;  
  if(like.shear_shear==1) {
    set_data_shear(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    //printf("ggl\n");
    set_data_ggl(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    //printf("clustering\n");
    set_data_clustering(like.Ncl,ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }

  if(like.clusterN==1){ 
    set_data_cluster_N(pred,start);
    start= start+tomo.cluster_Nbin*Cluster.N200_Nbin;
  }
  if(like.clusterWL==1){
    set_data_cgl(ell_Cluster,pred, start);
  }
  FILE *F;
  char filename[300];
  if (strstr(details,"FM") != NULL){
    sprintf(filename,"%s",details);
  }
  else {sprintf(filename,"datav_MCMC/%s_%s_%s",survey.name,like.probes,details);}
  F=fopen(filename,"w");
  for (i=0;i<like.Ndata; i++){  
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);
}

double write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in)
{
  compute_data_vector(details, ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,
    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,
    in.lf[0], in.lf[1], in.lf[2], in.lf[3], in.lf[4], in.lf[5], 
    in.m_lambda[0], in.m_lambda[1], in.m_lambda[2], in.m_lambda[3],
    in.m_lambda[4], in.m_lambda[5]);
  return 0;
}

double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in)
{
  double like = log_multi_like(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,

    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,
    in.lf[0], in.lf[1], in.lf[2], in.lf[3], in.lf[4], in.lf[5], 
    in.m_lambda[0], in.m_lambda[1], in.m_lambda[2], in.m_lambda[3],
    in.m_lambda[4], in.m_lambda[5]);
  
  return like;
}

 int main(void)
{
  clock_t begin, end;
  double time_spent;

/* here, do your time-consuming job */

  // init_cosmo();
  // //init_fisher_precision();
  // init_binning_fourier(20,20.0,15000.0,3000.0,21.0,5,5);
  // init_survey("LSST_Y10");
  // init_galaxies("zdistris/zdistri_model_z0=1.100000e-01_beta=6.800000e-01_Y10_source","zdistris/zdistri_model_z0=2.800000e-01_beta=9.000000e-01_Y10_lens", "gaussian", "gaussian", "SRD");  
  // init_clusters();
  // init_IA("NLA_HF", "GAMA");
  // init_priors("none","none","none","none");
  // //init_probes("pos_pos");
  // init_probes("pos_pos");
  // // init_probes("all_2pt");
  // init_Pdelta("halofit",0.8,0.35);
  // init_Pdelta("linear",0.8,0.35);
  // compute_data_vector("fid_IA",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.376695e+00,1.451179e+00,1.528404e+00,1.607983e+00,1.689579e+00,1.772899e+00,1.857700e+00,1.943754e+00,2.030887e+00,2.118943e+00,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.03,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  // compute_data_vector("fid_IA_sc0006_zdep",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.376695e+00,1.451179e+00,1.528404e+00,1.607983e+00,1.689579e+00,1.772899e+00,1.857700e+00,1.943754e+00,2.030887e+00,2.118943e+00,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.03,-0.00491861688284,-0.00320013979348,-0.00181278284029,8.59238903942e-05,0.006,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  // compute_data_vector("fid_IA_sigph2x",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.376695e+00,1.451179e+00,1.528404e+00,1.607983e+00,1.689579e+00,1.772899e+00,1.857700e+00,1.943754e+00,2.030887e+00,2.118943e+00,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.06,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  // compute_data_vector("fid_IA_sph005_lph003",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.376695e+00,1.451179e+00,1.528404e+00,1.607983e+00,1.689579e+00,1.772899e+00,1.857700e+00,1.943754e+00,2.030887e+00,2.118943e+00,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  
  
  
  init_cosmo();
  //init_fisher_precision();
  init_binning_fourier(20,20.0,15000.0,3000.0,21.0,5,5);
  init_survey("LSST_Y1");
  init_galaxies("zdistris/zdistri_model_z0=1.300000e-01_beta=7.800000e-01_Y1_source","zdistris/zdistri_model_z0=2.600000e-01_beta=9.400000e-01_Y1_lens", "gaussian", "gaussian", "SRD");
  init_clusters();
  init_IA("NLA_HF", "GAMA");
  init_priors("none","none","none","none");
  //init_probes("pos_pos");
  init_probes("pos_pos");
  // init_probes("all_2pt");
  init_Pdelta("halofit",0.8,0.35);
  // init_Pdelta("linear",0.8,0.35);
  compute_data_vector("fid_IA",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.562362e+00,1.732963e+00,1.913252e+00,2.100644e+00,2.293210e+00,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.03,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  compute_data_vector("fid_IA_sc0006_zdep",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.562362e+00,1.732963e+00,1.913252e+00,2.100644e+00,2.293210e+00,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.03,-0.00508007596207,-0.00365009514776,-0.00255711630406,-0.00109525588403,0.006,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  compute_data_vector("fid_IA_sigph2x",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.562362e+00,1.732963e+00,1.913252e+00,2.100644e+00,2.293210e+00,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.06,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  compute_data_vector("fid_IA_sph005_lph003",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.562362e+00,1.732963e+00,1.913252e+00,2.100644e+00,2.293210e+00,1.0,1.0,1.0,1.0,1.0,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  
  // // //init_data_inv("cov/cov_LSST_2.600000e+01_1.800000e+04_Rmin10_Ncl25_Ntomo10_2pt_clusterN_clusterWL_inv","datav/LSST_all_2pt_clusterN_clusterWL_fid");
  

  // begin = clock();
  //log_multi_like(0.3156 ,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  // end = clock();
  // time_spent = (double)(end - begin) / CLOCKS_PER_SEC;      
  // printf("timespent %le\n",time_spent);
  
  return 0;
}


