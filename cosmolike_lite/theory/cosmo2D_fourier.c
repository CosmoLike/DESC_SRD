#include "pt.c"

double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
double W_gal(double a, double nz); //complete weight for galaxy statistics
double W_HOD(double a, double nz); //galaxy weigth without bias factor (for projecting P_gg instead of P_nl)

double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj
double C_cl_lin_nointerp(double l, int ni, int nj);// galaxy clustering linear power spectrum of galaxies in bins ni
double C_cl_tomo_nointerp(double l, int ni, int nj);
double C_cl_HOD(double l, int ni);  //galaxy clustering power spectrum of galaxies in bin ni, using HOD model

double C_gl_tomo(double l, int ni, int nj);  //G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj);
double C_gl_HOD_tomo(double l, int ni, int nj);  //G-G lensing power spectrum from HOD model, lens bin ni, source bin nj

double C_shear_tomo(double l, int ni, int nj); //shear tomography power spectra
double C_shear_tomo_nointerp(double l, int ni, int nj);
/**********************************************************/

double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double dchi_da(double a){
  return 1./(a*a*hoverh0(a));
}
double W_kappa(double a, double fK, double nz){
  double wkappa = 1.5*cosmology.Omega_m*fK/a*g_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= (1.+MG_Sigma(a));
  }
  return wkappa;
}
double W_gal(double a, double nz){
  return gbias.b1_function(1./a-1.,(int)nz)*pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
}
double W_HOD(double a, double nz){
  return pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

/*********** Limber approximation integrands for angular power spectra *************/
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra

double int_for_C_cl_tomo_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  double s4 = PT_sigma4(k);
  res=W_HOD(a,ar[0])*W_HOD(a,ar[1])*dchi_da(a)/fK/fK;
  //if (b1*b1*Pdelta(k,a)+g4*(b1*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)) <0.)printf("%e %e %e  %e %e %e  %e\n",k/cosmology.coverH0,1./a-1.,b1,Pdelta(k,a),g4*PT_d2d2(k),g4*2.*s4, b1*b1*Pdelta(k,a)+g4*(b1*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)));
  res= res*(b1*b1*Pdelta(k,a)+g4*(b1*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)+b1*bs2*PT_d1s2(k)+0.5*b2*bs2*(PT_d2s2(k)-4./3.*s4)+.25*bs2*bs2*(PT_s2s2(k)-8./9.*s4)));
  return res;
}

double int_for_C_cl_tomo(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a);
  return res;
}
double int_for_C_cl_lin(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*p_lin(k,a);
  return res;
}
double int_for_C_cl_HOD(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_HOD(a,ar[0])*dchi_da(a)/fK/fK;
  if (res !=0){res= res*P_gg(k,a,(int)ar[0]);}
  return res;
}
double int_for_C_gl_tomo_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*(b1*Pdelta(k,a)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)));
  return res;
}
double int_for_C_gl_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_gl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_gal(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a);
  return res;
}

double int_for_C_gl_HOD_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  if (res !=0){res= res*P_gm(k,a, (int) ar[0]);}
  return res;
}

double int_for_C_shear_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_kappa(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a); 
  return res;
}

/*********** angular power spectra - without look-up tables ******************/

double C_cl_tomo_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ static int init =-1;
  double array[3] = {1.0*ni,1.0*nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
    if (ni != nj){
        if (init ==-1){
          printf("\nCalled C_cl(l,z1=%d,z2=%d) with non-linear bias parameters set.\n",ni,nj);
          printf("Cross-clustering beyond linear bias for cross-tomography bins not yet supported.\n");
          printf("Use linear bias only for z1!=z2 clustering\n\n");
          init = 1;}
          return int_gsl_integrate_high_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
        }
    return int_gsl_integrate_medium_precision(int_for_C_cl_tomo_b2,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
  }
  else if (ni == nj){  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);}
  return int_gsl_integrate_high_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
}

double C_cl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
  double array[3] = {1.0*ni,1.0*nj,l};
  return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
}

double C_gl_tomo_nointerp(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  double array[3] = {(double)ni,(double)nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
    return int_gsl_integrate_low_precision(int_for_C_gl_tomo_b2,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  }
  return int_gsl_integrate_low_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
}

double C_shear_tomo_nointerp(double l, int ni, int nj) //shear tomography power spectra of source galaxy bins ni, nj
{
  double array[3] = {(double) ni, (double) nj,l};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_medium_precision(int_for_C_shear_tomo,(void*)array,amin_source(j),amax_source(k),NULL,1000);
}

/*********** angular power spectra - with look-up tables ******************/

double C_cl_tomo(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.clustering_Nbin){
    printf("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range\nEXIT\n",ni,nj); exit(1);
  }

  if (recompute_clustering(C,G,N,ni,nj))
  {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
  
    double llog,res;
    int i,j,l1,l2,k;

    for (k=0; k<tomo.clustering_Nbin; k++) {
      for (j=k; j<tomo.clustering_Nbin; j++) {
        l1 = k*tomo.clustering_Nbin + j;
        l2 = j*tomo.clustering_Nbin + k;
        llog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          res = C_cl_tomo_nointerp(exp(llog),k,j);
        //  if (res < 0){printf("negative clustering power spectrum C_cl(%e,%d)- error\n", exp(llog),k); res =0.;}
          table[l1][i]= log(C_cl_tomo_nointerp(exp(llog),k,j));
          table[l2][i]=table[l1][i];
        }
      }
    }
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  double f1 = exp(interpol(table[ni*tomo.clustering_Nbin + nj], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}

double C_cl_HOD(double l, int ni)  //galaxy clustering power spectrum of galaxies in bin ni, using HOD model
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin){
    printf("Invalid tomography option in C_gl_HOD(l,%d)\nEXIT\n",ni);
    exit(EXIT_FAILURE);
  }
  if (recompute_clustering(C,G,N,ni,ni)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    double llog,array[3];
    int i,k;

    for (k=0; k<tomo.clustering_Nbin; k++){
      array[0]=(double) k;
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds){
        array[2] = exp(llog);
        table[k][i] =log(int_gsl_integrate_low_precision(int_for_C_cl_HOD,(void*)array,amin_lens(k),amax_lens(k),NULL,1000));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  
  double f1 = exp(interpol(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}
double C_gl_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k][i]= log(C_gl_tomo_nointerp(exp(llog),ZL(k),ZS(k)));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = exp(interpol(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.,1.));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_gl_HOD_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double llog,array[3];
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      array[0]=(double) ZL(k); array[1]=(double) ZS(k);
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(int_gsl_integrate_low_precision(int_for_C_gl_HOD_tomo,(void*)array,amin_source(ZS(k)),amax_lens(ZL(k)),NULL,1000));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = exp(interpol(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.,1.));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_shear_tomo(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
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
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog),Z1(k),Z2(k)));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}
