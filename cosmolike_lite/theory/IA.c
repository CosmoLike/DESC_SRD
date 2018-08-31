#include "IA.h"
void set_LF_GAMA(void);
void set_LF_DEEP2(void);
int check_LF (void);
double M_abs(double mag, double a); //absolute magnitude corresponding to aparent magnitude at a; h = 1 units, incl k-corrections
double n_all_LF(double mag, double a);
double f_red_LF(double mag, double a);
double A_LF(double mag, double a);

double W_source(double a, double nz);

double A_IA_Joachimi(double a);
double P_II (double k, double a);
double P_dI (double k, double a);
double P_II_JB(double k, double a);
double P_dI_JB(double k, double a);
/*** fast routines for total shear-shear + ggl signal for like.IA = 1,3****/
double C_shear_shear_IA(double s, int ni, int nj);
double C_ggl_IA(double s, int ni, int nj);
/*** tabulated versions of these for angular correlation functions ****/
double C_shear_shear_IA_tab(double l, int ni, int nj);
double C_ggl_IA_tab(double l, int ni, int nj);

double C_gI_nointerp(double s, int ni, int nj); //gI term using NL-LA P_deltaI power spectrum
double C_gI_JB_nointerp(double s, int ni, int nj); //gI term using Blazek+ P_deltaI power spectrum
double C_gI_lin_nointerp(double s, int ni, int nj);

double C_II_nointerp(double s, int ni, int nj);
double C_II_JB_nointerp(double s, int ni, int nj);
double C_II_lin_nointerp(double s, int ni, int nj);

double C_GI_nointerp(double s, int ni, int nj);
double C_GI_JB_nointerp(double s, int ni, int nj);
double C_GI_lin_nointerp(double s, int ni, int nj); 


/******** select LF ************/
void set_LF_GAMA(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_GAMA[0][i];
    LF_coefficients[1][i] =LF_coefficients_GAMA[1][i];
  }
}
void set_LF_DEEP2(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_DEEP2[0][i];
    LF_coefficients[1][i] =LF_coefficients_DEEP2[1][i];
  }
}
/************** normalization rouintes *********************/
double M_abs(double mag, double a){ //in h = 1 units, incl. Poggianti 1997 k+e-corrections
  static double *table;
  static double dz;
  double ke;
  if (table ==0){
    int i;
    //read in + tabulate k+e corrections for early types, restframe r band
    //interpolated from http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A%2BAS/122/399
    table = create_double_vector(0,30);
    dz = 0.1;
    /*FILE *ein;
    double d1,d2,d3;
    ein = fopen(survey.Kcorrect_File,"r");
    EXIT_MISSING_FILE(ein, "LoverL0",survey.Kcorrect_File)*/
    for (i = 0; i< 31; i++){
   /*   fscanf(ein,"%le %le %le\n",&d1,&d2,&d3);
       table[i] = d2+d3; */
      table[i] = KE[i];
    }
   /* fclose(ein);*/
  }
  double z=1./a-1.0;
  if(z >= 3.0) {z=2.99;}  // no acceptable k-korrection exists for k>3, also no meaningful IA model
  ke = interpol(table, 31, 0., 3.0, 0.1, z, 1.0, 1.0);
  
  return mag - 5.0*log10(f_K(chi(a))/a*cosmology.coverH0) -25.0 -ke;
}
double n_all_LF(double mag, double a){
   double Q, alpha, Mstar, Phistar,P,Mlim, LF_all[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
  
  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];
  
  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units
  
   
  return LF_all[0]*gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double f_red_LF(double mag, double a){
  double Q, alpha, Mstar, Phistar,P,Q_red, alpha_red, Mstar_red, Phistar_red,P_red,Mlim, LF_all[3], LF_red[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //red galaxies
  alpha_red = LF_coefficients[0][1] + nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;
  P_red = LF_coefficients[0][4]+nuisance.LF_red_P;
  Phistar_red = LF_coefficients[0][0];
  
  LF_red[0] = Phistar_red*pow(10.0,0.4*P_red*(1./a-1));
  LF_red[1] = Mstar_red-Q_red*(1./a-1. - 0.1);
  LF_red[2]= alpha_red;
  
  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];
  
  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units
  
  return LF_red[0]/LF_all[0]*gsl_sf_gamma_inc(LF_red[2]+1,pow(10.0,-0.4*(Mlim-LF_red[1])))/gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double A_LF(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
    
    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //red galaxies
  alpha_red = LF_coefficients[0][1]+ nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;
  
    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;
  
    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


double A_LF_all(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
    
    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //all galaxies
  alpha_red = LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar_red = LF_coefficients[1][2]; //in h = 1 units
  Q_red = LF_coefficients[1][3]+nuisance.LF_Q;
  
    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;
  
    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


int check_LF(void){ //return 1 if combination of all + red galaxy LF parameters is unphysical, i.e. if f_red > 1 for some z < redshift.shear_zdistrpar_zmax
  double a=1./(1+redshift.shear_zdistrpar_zmax)+0.005;
  while (a < 1.){
    if( M_abs(survey.m_lim,a) < LF_coefficients[1][2]-(LF_coefficients[1][3]+nuisance.LF_Q)*(1./a-1. - 0.1) || M_abs(survey.m_lim,a) < LF_coefficients[0][2]-(LF_coefficients[0][3]+nuisance.LF_red_Q)*(1./a-1. - 0.1)){return 1;}
    if (f_red_LF(survey.m_lim,a) > 1.0){return 1;}
    a+=0.01;
  }
  return 0;
}

/*=========================================================*/
/*=============  Intrinsic Alignment models  ==============*/
/*=========================================================*/

double W_source(double a, double nz){
  return zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a);
}


double A_IA_Joachimi(double a){
  double z, A_red, highz = 0.75;
  z = 1./a-1;
  A_red = nuisance.A_ia*A_LF(survey.m_lim,a)*f_red_LF(survey.m_lim,a); //A_0*<(L/L_0)^beta>*f_red
  if (a < 1./(1.+highz)){ //z > highz, factor in uncertainty in extrapolation of redshift scaling
    return A_red*pow((1.+z)/nuisance.oneplusz0_ia,nuisance.eta_ia)*pow((1.+z)/(1.+highz),nuisance.eta_ia_highz);
  }
  else{//standard redshift scaling
    return A_red*pow((1.+z)/nuisance.oneplusz0_ia,nuisance.eta_ia);
  }
}

double P_II (double k, double a){
  return pow(A_IA_Joachimi(a),2.0)*pow(cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a),2.0)*Pdelta(k,a); //Joachimi+ 11 Eq.(23) + Eq. (B.6)
}
double P_dI (double k, double a){
  return A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*Pdelta(k,a); // Joachimi+ 11 Eq.(6)
}

double P_II_lin (double k, double a){
  return pow(A_IA_Joachimi(a),2.0)*pow(cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a),2.0)*p_lin(k,a); 
}

double P_dI_lin (double k, double a){
  return A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*p_lin(k,a); 
}

double P_II_JB(double k, double a){
  static cosmopara C;
   
  static double **table_P =0;
  static double da =0., amin = 0., amax =1.;
  static double dlgk = 0., lgkmin = 0., lgkmax = 0.;
  static int N_k = 500;
  
  if (recompute_cosmo3D(C)){
    int i,j;
    double kk,aa,Da4,b1,Az,sig2,ksm,sk,A0,a0,Da0;
    ksm = 0.1/cosmology.coverH0;
    amin = 1./(0.99+redshift.shear_zdistrpar_zmax);
    amax= 1./(1.001+redshift.shear_zdistrpar_zmin);
    da = (amax - amin)/(Ntable.N_a-1.);

    lgkmin = log(IA_JB_terms[0][0]*cosmology.coverH0);
    lgkmax = log(IA_JB_terms[N_k-1][0]*cosmology.coverH0);
    dlgk = (lgkmax - lgkmin)/(N_k-1.);
    a0 = 1./1.3;
    Da0 = growfac(a0)/growfac(0.9999);
    A0 = A_IA_Joachimi(a0)/(1.0+3.1*Da0*Da0)*cosmology.Omega_m*nuisance.c1rhocrit_ia;
    if (table_P==0){
      table_P = create_double_matrix(0, Ntable.N_a-1, 0, N_k-1);
    }
    aa= amin;
    for (i=0; i<Ntable.N_a; i++, aa += da) {
      Da4 = pow(growfac(aa)/growfac(0.9999),4.0);
      sig2 = 3.1*sqrt(Da4);//assuming k_sm = 0.1 h/Mpc, as in Jonathan's email
      Az = A0*Da0/(growfac(aa)/growfac(0.9999))*(1.0+3.1*Da0*Da0)/sig2; //Normalization for matching GI on large scales
      b1 = b_source(aa);
      for (j = 0; j < N_k; j++){
        kk = IA_JB_terms[j][0]*cosmology.coverH0;
        sk = exp(-pow(kk/ksm,2.0));
        table_P[i][j] = log(Az*Az*(Pdelta(kk,aa)*sk*sk+b1/pow(cosmology.coverH0,3.0)*(sk*(IA_JB_terms[j][4]*Da4+IA_JB_terms[j][5]*Da4+58./105.*sig2*p_lin(kk,aa))+b1*IA_JB_terms[j][10]*Da4)));
      }
    }
    update_cosmopara(&C);
//    printf("P_II_JB tabulated\n");
  }
  double lgk = log(k);
  if (lgk <= lgkmin || lgk >= lgkmax ||amax >= a || amin <= a){ return pow(cosmology.Omega_m*nuisance.c1rhocrit_ia/a,2.0)*Pdelta(k,a)*pow(f_red_LF(survey.m_lim,a),2.0);}
  return exp(interpol2d(table_P, Ntable.N_a, amin, amax, da, a, N_k, lgkmin, lgkmax, dlgk, lgk, 0.0, 0.0));
}


double P_dI_JB(double k, double a){
  static cosmopara C;
  
  static double **table_P =0;
  static double da =0., amin = 0., amax =1.;
  static double dlgk = 0., lgkmin = 0., lgkmax = 0.;
  static int N_k = 500;
  
  if (recompute_cosmo3D(C)){
    int i,j;
    double kk,aa,Da4,b1,Az,sig2,ksm,A0,a0,Da0;;
    ksm = 0.1/cosmology.coverH0;
//    printf("calculating P_II_JB\n");
    amin = 1./(0.99+redshift.shear_zdistrpar_zmax);
    amax= 1./(1.001+redshift.shear_zdistrpar_zmin);
    da = (amax - amin)/(Ntable.N_a-1.);
    
    lgkmin = log(IA_JB_terms[0][0]*cosmology.coverH0);
    lgkmax = log(IA_JB_terms[N_k-1][0]*cosmology.coverH0);
    dlgk = (lgkmax - lgkmin)/(N_k-1.);
    
    if (table_P==0){
      table_P = create_double_matrix(0, Ntable.N_a-1, 0, N_k-1);
    }
    a0 = 1./1.3;
    Da0 = growfac(a0)/growfac(0.9999);
    A0 = A_IA_Joachimi(a0)/(1.0+3.1*Da0*Da0)*cosmology.Omega_m*nuisance.c1rhocrit_ia;

    aa= amin;
    for (i=0; i<Ntable.N_a; i++, aa += da) {
      Da4 = pow(growfac(aa)/growfac(0.9999),4.0);
      sig2 = 3.1*sqrt(Da4);//assuming k_sm = 0.1 h/Mpc, as in Jonathan's email
      Az = A0*Da0/(growfac(aa)/growfac(0.9999))*(1.0+3.1*Da0*Da0)/sig2; //Normalization for matching GI on large scales
      b1 = b_source(aa);
      for (j = 0; j < N_k; j++){
        kk = IA_JB_terms[j][0]*cosmology.coverH0;
        table_P[i][j] = log(Az*(Pdelta(kk,aa)*exp(-pow(kk/ksm,2.0))+b1/pow(cosmology.coverH0,3.0)*(IA_JB_terms[j][4]*Da4+IA_JB_terms[j][5]*Da4+58./105.*sig2*p_lin(kk,aa))));
      }
    }
    update_cosmopara(&C);
  }
  double lgk = log(k);
  if (lgk <= lgkmin || lgk >= lgkmax ||amax >= a || amin <= a){ return cosmology.Omega_m*nuisance.c1rhocrit_ia/a*Pdelta(k,a)*f_red_LF(survey.m_lim,a);}

  return exp(interpol2d(table_P, Ntable.N_a, amin, amax, da, a, N_k, lgkmin, lgkmax, dlgk, lgk, 0.0, 0.0));
}
/*=========================================================*/
/*============= Intrinsic Alignment II-Term ==============*/
/*=========================================================*/


double int_for_C_II(double a, void *params)
{
  double res, ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
   ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  res= W_source(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_II(k,a);
  return res;
}
     
double C_II_nointerp(double s, int ni, int nj) 
{
  if(ni!=nj && redshift.shear_photoz==0) return 0.0;
  if (abs(ni-nj) >1 && redshift.shear_photoz==3) return 0.0;

  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_low_precision(int_for_C_II,(void*)array,fmax(amin_source(ni),amin_source(nj)),fmin(amax_source_IA(ni),amax_source_IA(nj)),NULL,1000); //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slow down for ni != nj with redshift.shear_photoz ==3
}


double int_for_C_II_lin(double a, void *params)
{
  double res, ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
   ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  res= W_source(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_II_lin(k,a);
  return res;
}
     
double C_II_lin_nointerp(double s, int ni, int nj) 
{
  if(ni!=nj && redshift.shear_photoz==0) return 0.0;
  if (abs(ni-nj) >1 && redshift.shear_photoz==3) return 0.0;

  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_low_precision(int_for_C_II_lin,(void*)array,fmax(amin_source(ni),amin_source(nj)),fmin(amax_source_IA(ni),amax_source_IA(nj)),NULL,1000); //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slow down for ni != nj with redshift.shear_photoz ==3
}


double int_for_C_II_JB(double a, void *params)
{
  double res, ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_source(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_II_JB(k,a);
  return res;
}

double C_II_JB_nointerp(double s, int ni, int nj)
{
  if(ni!=nj && redshift.shear_photoz==0) return 0.0;

  if (abs(ni-nj) >1 && redshift.shear_photoz==3) return 0.0;

  double array[3] = {(double) ni, (double) nj,s};
  return int_gsl_integrate_low_precision(int_for_C_II_JB,(void*)array,fmax(amin_source(ni),amin_source(nj)),fmin(amax_source_IA(ni),amax_source_IA(nj)),NULL,1000); //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slow down for ni != nj with redshift.shear_photoz ==3
}

/*=========================================================*/
/*============= gI routines ==============*/
/*=========================================================*/

double int_for_C_gI(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = W_gal(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI(k,a);
  return res;
}

double C_gI_nointerp(double s, int ni, int nj)
{
  if (amax_source_IA(nj)< amin_lens(nj)) return 0.; // source bin and lens bin don't overlap, no gI term
  double array[3] = {(double) ni, (double) nj,s};
  return -int_gsl_integrate_low_precision(int_for_C_gI,(void*)array,fmax(amin_source(nj),amin_lens(ni)),fmin(amax_source_IA(nj),amax_lens(ni)),NULL,1000);  //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slowdown for ni != nj with redshift.shear_photoz ==3
}

double int_for_C_gI_lin(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = W_gal(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI_lin(k,a);
  return res;
}

double C_gI_lin_nointerp(double s, int ni, int nj)
{
  if (amax_source_IA(nj)< amin_lens(nj)) return 0.; // source bin and lens bin don't overlap, no gI term
  double array[3] = {(double) ni, (double) nj,s};
  return -int_gsl_integrate_low_precision(int_for_C_gI_lin,(void*)array,fmax(amin_source(nj),amin_lens(ni)),fmin(amax_source_IA(nj),amax_lens(ni)),NULL,1000);  //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slowdown for ni != nj with redshift.shear_photoz ==3
}

double int_for_C_gI_JB(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = W_gal(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI_JB(k,a);
  return res;
}

double C_gI_JB_nointerp(double s, int ni, int nj)
{
  if (amax_source_IA(nj)< amin_lens(nj)) return 0.; // source bin and lens bin don't overlap, no gI term
  double array[3] = {(double) ni, (double) nj,s};
  return -int_gsl_integrate_low_precision(int_for_C_gI_JB,(void*)array,fmax(amin_source(nj),amin_lens(ni)),fmin(amax_source_IA(nj),amax_lens(ni)),NULL,1000);  //use low precision to avoid slowdown for ni != nj with redshift.shear_photoz ==3
}
/*=========================================================*/
/*============= GI routines ==============*/
/*=========================================================*/

double int_for_C_GI(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = (W_source(a,ar[0])*W_kappa(a,fK,ar[1])+W_source(a,ar[1])*W_kappa(a,fK,ar[0]))*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI(k,a);
  return res;
}

double C_GI_nointerp(double s, int ni, int nj) 
{
  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return -int_gsl_integrate_low_precision(int_for_C_GI,(void*)array,amin_source(j),amax_source_IA(k),NULL,1000);
}


double int_for_C_GI_lin(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = (W_source(a,ar[0])*W_kappa(a,fK,ar[1])+W_source(a,ar[1])*W_kappa(a,fK,ar[0]))*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI_lin(k,a);
  return res;
}

double C_GI_lin_nointerp(double s, int ni, int nj) 
{
  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return -int_gsl_integrate_low_precision(int_for_C_GI_lin,(void*)array,amin_source(j),amax_source_IA(k),NULL,1000);
}


double int_for_C_GI_JB(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res= (W_source(a,ar[0])*W_kappa(a,fK,ar[1])+W_source(a,ar[1])*W_kappa(a,fK,ar[0]))*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI_JB(k,a);
  return res;
}

double C_GI_JB_nointerp(double s, int ni, int nj) 
{
  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return -int_gsl_integrate_low_precision(int_for_C_GI_JB,(void*)array,amin_source(j),amax_source_IA(k),NULL,1000);
}

/******** faster shear/ggl + IA routines *********/

double int_for_C_ggl_IA(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= W_gal(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}
double int_for_C_ggl_IA_Az(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= W_gal(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double int_for_C_ggl_IA_mpp(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  res= W_HOD(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  return res*(gbias.b1_function(1./a-1.,(int)ar[0])*Pdelta(k,a))*dchi_da(a)/fK/fK;
  
//  return res*(gbias.b1_function(1./a-1.,(int)ar[0])*Pdelta(k,a)+P_gm_rm(k,a))*dchi_da(a)/fK/fK;
}

double int_for_C_ggl_IA_mpp_b2(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  res= W_HOD(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  return res*(gbias.b1_function(1./a-1.,(int)ar[0])*Pdelta(k,a)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)))*dchi_da(a)/fK/fK;
}

double int_for_C_ggl_IA_Az_b2(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= W_HOD(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  
  return res*(b1*Pdelta(k,a)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)))*dchi_da(a)/fK/fK;
}

double int_for_C_shear_shear_IA(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);
  norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm+wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double int_for_C_shear_shear_IA_Az(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0])*nuisance.A_z[(int)ar[0]];
  ws2 = W_source(a,ar[1])*nuisance.A_z[(int)ar[1]];
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm+wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double int_for_C_shear_shear_IA_mpp(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  res= ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm+wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double C_shear_shear_IA(double s, int ni, int nj)
{
 double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}

  switch(like.IA){
    case 1: return int_gsl_integrate_medium_precision(int_for_C_shear_shear_IA,(void*)array,amin_source(j),amax_source(k),NULL,1000);
    case 3: return int_gsl_integrate_medium_precision(int_for_C_shear_shear_IA_Az,(void*)array,amin_source(j),amax_source(k),NULL,1000);
    case 4: return int_gsl_integrate_medium_precision(int_for_C_shear_shear_IA_mpp,(void*)array,amin_source(j),amax_source(k),NULL,1000);
    default: { printf("IA.c: C_shear_shear_IA does not support like.IA = %d\nEXIT\n", like.IA);
               exit(1);
    }
  }
}

double C_ggl_IA(double s, int nl, int ns)
{
  double array[3] = {(double) nl, (double) ns,s};
  switch(like.IA){
    case 1: return int_gsl_integrate_medium_precision(int_for_C_ggl_IA,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
    case 3: if (gbias.b2[nl]) return int_gsl_integrate_low_precision(int_for_C_ggl_IA_Az_b2,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
            return int_gsl_integrate_low_precision(int_for_C_ggl_IA_Az,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
    case 4: if (gbias.b2[nl]) return int_gsl_integrate_low_precision(int_for_C_ggl_IA_mpp_b2,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
            return int_gsl_integrate_low_precision(int_for_C_ggl_IA_mpp,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
    default: printf("IA.c: C_ggl_IA does not support like.IA = %d\nEXIT\n", like.IA); exit(1);
  }
}
/********** tabulated versions of C_ggl_IA and C_shear_shear_IA******/
/************for calculation of angular corelation functions ********/
double C_ggl_IA_tab(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
  //  printf("calculating C_ggl_IA  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      sig = create_double_vector(0,tomo.ggl_Npowerspectra-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double res,llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;

      sig[k] = 1.;
      osc[k] = 0;
      res = C_ggl_IA(500.,ZL(k),ZS(k));
      if (res < 0){sig[k] = -1.;}
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i] = C_ggl_IA(exp(llog),ZL(k),ZS(k));
        if (res*sig[k] <0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0){
        for(i = 0; i < Ntable.N_ell; i++){
            res = table[k][i];
            table[k][i] = log(sig[k]*res);}
      }
      
    }

    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  int k = N_ggl(ni,nj);
  double f1 = 0.;
  if(test_zoverlap(ni,nj) && osc[k] ==0 ){f1 = sig[k]*exp(interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));}
  if(test_zoverlap(ni,nj) && osc[k] ==1 ){f1 = interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_shear_shear_IA_tab(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(C_shear_shear_IA(exp(llog),Z1(k),Z2(k)));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}