// lens efficiencies
double g_cmb (double a);//lens efficiency for CMB lensing
double g_tomo(double a, int zbin); // lens efficiency of galaxies in tomography bin zbin

//shear routines
double zdistr_photoz(double zz,int j); //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
//double zdistr_gaussian_ztrue (double zz, int j);

//clustering routines
double pf_photoz(double zz,int j); //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

//specific photo-z models

double sigma_zphot_shear (double z, int nz); //compute photo-z scatter (sigma(z)) from nuisance.sigma_zphot_shear parameters
double bias_zphot_shear (double z, int nz);  //compute photo-z bias (b_{zph}(z)) from nuisance.bias_zphot_shear parameters

double sigma_zphot_clustering (double z, int nz);
double bias_zphot_clustering (double z, int nz);


//lens + source galaxy counts
double int_nsource(double z, void* param);
double int_nlens(double z, void* param);
double nsource(int j); //returns n_gal for shear tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double nlens(int j); //returns n_gal for clustering tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double zmean(int j);//mean true redshift of (clustering/lens) galaxies within redshift bin j
double zmean_source(int j); //mean true redshift of source galaxies in tomography bin j


/***** redshift integration boundaries **********/
double amin_source(int i);
double amax_source(int i);
double amax_source_IA(int i);
double amin_lens(int i);
double amax_lens(int i);

//routines for association of a pair redshift bin numbers and power spectrum tomography bin
int test_kmax(double l, int zl); //valid l + lens/clustering redshift combination?
int test_zoverlap_cov(int zl, int zs); //restrict NG covaraince computation to source behind lens configurations
int test_zoverlap(int zl, int zs); //valid lens + source redshift bin combination?
int test_zoverlap_c (int zc, int zs); //valid cluster + source redshift bin combination?
int N_ggl(int zl, int zs);//(z_l,z_s) -> N_tomo_ggl
int ZL(int Nbin);//N_tomo_ggl -> z_l
int ZS(int Nbin);//N_tomo_ggl -> z_s
int N_cgl(int zl, int zs);//(z_c,z_s) -> N_tomo_cgl
int ZC(int Nbin);//N_tomo_cgl -> z_c
int ZSC(int Nbin);//N_tomo_cgl -> z_s
int N_shear (int z1, int z2);//(z_1,z_2)-> N_tomo_shear, assuming z1<=z2
int Z1(int Nbin);//N_tomo_shear -> z1, assuming z1<=z2
int Z2(int Nbin);//N_tomo_shear -> z2, assuming z1<=z2
void write_gglensing_zbins(char *surveyname);
double ggl_efficiency(int zl, int zs);
///
double g_bg (double a, int nzlens);//no longer supported - declaration only to prevent compile errors


/********** integration boundary routines *************/
double amin_source(int i){
  if (i == -1 || redshift.shear_photoz == 1 || redshift.shear_photoz == 2|| redshift.shear_photoz == 4){return 1./(redshift.shear_zdistrpar_zmax+1.);}
  if (redshift.shear_photoz == 0){ return 1./(1+tomo.shear_zmax[i]);}
  return 1./(1+fmin(tomo.shear_zmax[i] + 5.*nuisance.sigma_zphot_shear[i] + fabs(nuisance.bias_zphot_shear[i]),redshift.shear_zdistrpar_zmax));
}
double amax_source(int i){
return 1./(1.+fmax(redshift.shear_zdistrpar_zmin,0.001));
}

double amax_source_IA(int i){
  if (i == -1 || redshift.shear_photoz == 1 || redshift.shear_photoz == 2|| redshift.shear_photoz == 4){return 1./(1.+fmax(redshift.shear_zdistrpar_zmin,0.001));}
  if (redshift.shear_photoz == 0){ return 1./(1.+fmax(tomo.shear_zmin[i],0.001));}
  return 1./(1+fmax(tomo.shear_zmin[i] - 5.*nuisance.sigma_zphot_shear[i] - fabs(nuisance.bias_zphot_shear[i]),0.001));
}
double amin_lens(int i){
  if (i == -1 || redshift.clustering_photoz == 1 || redshift.clustering_photoz == 2){return 1./(redshift.clustering_zdistrpar_zmax+1.);}
  if (redshift.clustering_photoz == 0){ return 1./(1+tomo.clustering_zmax[i]);}
  if (redshift.clustering_photoz == 4){ return 1./(1+tomo.clustering_zmax[i]+2.*fabs(nuisance.bias_zphot_clustering[i]));}
  return 1./(1+fmin(tomo.clustering_zmax[i] + 5.*nuisance.sigma_zphot_clustering[i] + fabs(nuisance.bias_zphot_clustering[i]),redshift.clustering_zdistrpar_zmax));
}
double amax_lens(int i){
  if (i == -1 || redshift.clustering_photoz == 1 || redshift.clustering_photoz == 2){return 1./(1.+fmax(redshift.clustering_zdistrpar_zmin,0.001));}
  if (redshift.clustering_photoz == 0){ return 1./(1.+fmax(tomo.clustering_zmin[i],0.001));}
  if (redshift.clustering_photoz == 4){ return 1./(1+fmax(tomo.clustering_zmin[i]-2.*fabs(nuisance.bias_zphot_clustering[i]),0.001));}
  return 1./(1+fmax(tomo.clustering_zmin[i] - 5.*nuisance.sigma_zphot_clustering[i]-fabs(nuisance.bias_zphot_clustering[i]),0.001));
}

/************ redshift overlap tests, allowed tomography combinations **********/
int test_kmax(double l, int zl){ //test whether the (l,zl) bin is in the linear clustering regime - return 1 if true, 0 otherwise
  static double chiref[10] = {-1.};
  if (chiref[0] < 0){
    int i;
    for (i = 0; i < tomo.clustering_Nbin; i++){
      chiref[i] = chi(1./(1.+0.5*(tomo.clustering_zmin[i]+tomo.clustering_zmax[i])));
    }
  }
  double R_min = like.Rmin_bias; //set minimum scale to which we trust our bias model, in Mpc/h
  double kmax = constants.twopi/R_min*cosmology.coverH0;
  if((l+0.5)/chiref[zl] < kmax){ return 1;}
  return 0;
}





int test_zoverlap(int zl, int zs){ //test whether source bin zs is behind lens bin zl
  if (ggl_efficiency(zl,zs) > survey.ggl_overlap_cut) {return 1;}
  if (redshift.shear_photoz < 4 && tomo.clustering_zmax[zl] <= tomo.shear_zmin[zs]){return 1;}
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz != 4 && tomo.clustering_zmax[zl] < zmean_source(zs)){return 1;}
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz == 4 && zmean(zl)+0.1 < zmean_source(zs)){return 1;}
  return 0;
}
int test_zoverlap_cov(int zl, int zs){ //test whether source bin zs is behind lens bin zl
  if (redshift.shear_photoz < 4 && tomo.clustering_zmax[zl] <= tomo.shear_zmin[zs]){return 1;}
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz != 4 && tomo.clustering_zmax[zl] < zmean_source(zs)){return 1;}
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz == 4 && zmean(zl)+0.1 < zmean_source(zs)){return 1;}
  return 0;
}
int test_zoverlap_c(int zc, int zs){ //test whether source bin zs is behind lens bin zl
  if (redshift.shear_photoz < 4 && tomo.cluster_zmax[zc] <= tomo.shear_zmin[zs]){return 1;}
  if (redshift.shear_photoz == 4 && tomo.cluster_zmax[zc] < zmean_source(zs)){return 1;}
  return 0;
}

int N_ggl(int zl, int zs){
  static int N[10][10] = {-42};
  if (N[0][0] < 0){
    int i, j,n;
    n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[i][j] = n; n++;}
        else{N[i][j] = -1;}
      }
    }
  }
  return N[zl][zs];
}

void write_gglensing_zbins(char *surveyname){
  FILE *F1;
  char filename[400];
  sprintf(filename,"%s.gglensing_zbins",surveyname);
  printf("filename gglensing %s\n",filename);
  F1 = fopen(filename,"w");
  fprintf(F1,"#zl zs  accept?  overlap metric\n");
  fprintf(F1,"#pipeline supplied ggl_overlap_cut = %.3f\n",survey.ggl_overlap_cut);
  for (int i = 0; i < tomo.clustering_Nbin; i++){
    for(int j = 0; j<tomo.shear_Nbin;j++){
      fprintf(F1,"%d   %d   %d    %.3f\n",i,j,test_zoverlap(i,j),ggl_efficiency(i,j));
    }
  }
  fclose(F1);
}
int ZL(int Nbin){
  static int N[100] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[n] = i; n++;}
      }
    }
  }
  return N[Nbin];
}
int ZS(int Nbin){
  static int N[100] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[n] = j; n++;}
      }
    }
  }
  return N[Nbin];
}

int N_cgl(int zc, int zs){
  static int N[10][10] = {-42};
  if (N[0][0] < 0){
    int i, j,n;
    n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap_c(i,j)){N[i][j] = n; n++;}
        else{N[i][j] = -1;}
      }
    }
  }
  return N[zc][zs];
}
int ZC(int Nbin){
  static int N[100] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap_c(i,j)){N[n] = i; n++;}
      }
    }
  }
  return N[Nbin];
}
int ZSC(int Nbin){
  static int N[100] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap_c(i,j)){N[n] = j; n++;}
      }
    }
  }
  return N[Nbin];
}
int N_shear (int z1, int z2){ //find shear tomography bin number N_shear of tomography combination (z1,z2)
  static int N[10][10] = {-42};
  if (N[0][0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  return N[z1][z2];
}
int Z1(int Nbin){// find z1 of tomography combination (z1,z2) constituting shear tomography bin Nbin
  static int N[55] = {-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[n] = i;
        n++;
      }
    }
  }
  return N[Nbin];
}

int Z2(int Nbin){ // find z2 of tomography combination (z1,z2) constituting shear tomography bin Nbin
  static int N[55]={-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[n] = j;
        n++;
      }
    }
  }
  return N[Nbin];
}


/******************** routines for redshift distributions, including photo-zs (optional) ********************/
///////////////// Start shear routines /////////////////////
double sigma_zphot_shear (double z, int nz){
  return (1.+z)*nuisance.sigma_zphot_shear[nz];
}
double bias_zphot_shear (double z, int nz){
  return (1.+z)*nuisance.bias_zphot_shear[nz];
}

double zdistr_histo_n(double z,  void *params) // return nz(z,j) based on redshift file with structure z[i] nz[0][i] .. nz[tomo.shear_Nbin-1][i]
{
  double *array = (double*)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max,zhisto_min,dz;
  
  if (tab==0){
    double *z_v;
    int i,k,zbins;
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    tab=create_double_matrix(0,tomo.shear_Nbin-1,0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.shear_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){
        break;
      }
      for (k = 0; k < tomo.shear_Nbin; k++){
        fscanf(ein,"%le",&tab[k][i]);
      }
    }
    fclose(ein);
    dz = (z_v[i-1]-z_v[0])/(1.*i -1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    zbins = i;
    // now, set tomography bin boundaries
    for (k = 0; k < tomo.shear_Nbin; k++){
      double max = tab[k][0];
      for (i = 1; i < zbins; i++){
        if (tab[k][i]> max){max = tab[k][i];}
      }
      i = 0;
      while (tab[k][i] <1.e-8*max && i < zbins-2){i++;}
      tomo.shear_zmin[k] = z_v[i];
      i = zbins-1;
      while (tab[k][i] <1.e-8*max && i > 0){i--;}
      tomo.shear_zmax[k] = z_v[i];
      printf("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f\n",k,tomo.shear_zmin[k],k,tomo.shear_zmax[k]);
    }
    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.shear_zmax[tomo.shear_Nbin-1] || zhisto_min > tomo.shear_zmin[0]){
      printf("zhisto_min = %e,zhisto_max = %e\n",zhisto_min,zhisto_max);
      printf("tomo.shear_zmin[0] = %e, tomo.shear_zmax[N-1] = %e\n", tomo.shear_zmin[0],tomo.shear_zmax[tomo.shear_Nbin-1]);
      printf("Error in redshift.c:zdistr_histo_n: %s parameters incompatible with tomo.shear bin choice\nEXIT!\n",redshift.shear_REDSHIFT_FILE);
      exit(1);
    }
  }
  
  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)array[0]][(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}

double zdistr_histo_1(double z, void *params) //return nz(z) based on redshift file with one redshift distribution
{
  static double *tab =0;
  FILE *ein;
  
  static double zhisto_max,zhisto_min,dz;
  
  if (tab==0){
    double *z_v,space1,space2;
    int i,zbins;
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    tab=create_double_vector(0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.shear_REDSHIFT_FILE,"r");
    
    for (i=0;i<zbins;i++){
      fscanf(ein,"%le %le %le %le\n",&z_v[i],&space1,&space2,&tab[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
    }
    fclose(ein);
    dz = (z_v[i-1]-z_v[0])/(1.*i -1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    redshift.shear_zdistrpar_zmin = zhisto_min;
    redshift.shear_zdistrpar_zmax = zhisto_max;
    free_double_vector(z_v,0,zbins-1);
  }
  
  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}


double zdistr_photoz(double zz,int j) //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table = 0;
  static double da = 0.0;
  static double zhisto_max,zhisto_min;
  static nuisancepara N;
  static int zbins =2000;
  if (table ==0){ //force zmin, zmax, to match supplied file
    FILE *ein;
    double *z_v, space;
    int i,k,nz, N=3;
    if (redshift.shear_photoz ==4) {N = tomo.shear_Nbin;}
    nz = line_count(redshift.shear_REDSHIFT_FILE);
    z_v=create_double_vector(0, nz-1);
    ein=fopen(redshift.shear_REDSHIFT_FILE,"r");
    for (i=0;i<nz;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < N; k++){fscanf(ein,"%le",&space);}
    }
    fclose(ein);
    redshift.shear_zdistrpar_zmin = fmax(z_v[0],0.01);
    redshift.shear_zdistrpar_zmax = z_v[i-1] +(z_v[i-1]-z_v[0])/(1.*i-1.);
    //tomo.shear_zmin/max will be set when calling zdistr_histo_n
    //for (int i = 0; i < tomo.shear_Nbin; i++){
    //  tomo.shear_zmin[i] =redshift.shear_zdistrpar_zmin;
    //  tomo.shear_zmax[i] =redshift.shear_zdistrpar_zmax;
    //}
  }
  if (recompute_zphot_shear(N) || table==0){
    update_nuisance(&N);
    if (table ==0){
      table   = create_double_matrix(0, tomo.shear_Nbin, 0, zbins-1);
    }
    zhisto_max =redshift.shear_zdistrpar_zmax;
    zhisto_min = redshift.shear_zdistrpar_zmin;
    da = (zhisto_max-zhisto_min)/(1.0*zbins);
    double array[4], NORM[11],norm,zi,x1,x2,eta,outfrac;
    //the outlier fraction (outfrac) should be specified externally. This is a temporary hack.
    int i,k;
    for (i = 0; i< tomo.shear_Nbin; i++){
      array[0] =tomo.shear_zmin[i];
      array[1] =tomo.shear_zmax[i];
      norm = 0.0;
      outfrac=0.05;
      //eta = 0.68; //this is the approximate mixing parameter for the pseudo voigt for equal FWHM L and G distributions. Note than when using this formulation, the non-equivalent profile re-scalings must be accounted for below.
      eta = 0.736; //this is the approximate mixing parameter for the pseudo voigt for equally rescaled L and G distributions. Note that this leads to a slightly broader profile than when using equal FWHM L and G distributions (at a given gaussian sigma).
      // Currently, the input sigma parameter has not been rescaled. It still gives the width for the Gaussian component of the psuedo Voigt profile. In general, for a given FWHM, the corresponding sigma is smaller for the Voigt profile. For equal FWHM in the component profiles, eta = 0.68 and the effective FWHM is 1.64 that of the Guassian case. For equivalent rescalings, [which is the currently implemented case] eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.
      switch (redshift.shear_photoz){
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
          norm = int_gsl_integrate_medium_precision(zdistr_histo_1, (void*)array, array[0],array[1],NULL, 1024);
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            table[i+1][k] = 0;
            if (zi >= array[0] && zi <= array[1]){ table[i+1][k] = zdistr_histo_1(zi,(void*)array)/norm;}
          }
          break;
        case 1: //Pseudo Voigt (Lorentzian + Gaussian) photo-zs
          norm = 0.;
          if (sigma_zphot_shear(0.,i) == 0.){
            printf("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz =1, but nuisance.sigma_zphot_shear[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
            x2 = (array[1]-zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
      table[i+1][k] = zdistr_histo_1(zi,(void*)array)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))));
      // this creates a pseudo Voigt profile by adding a Gaussian and Lorentian (convolved with a tophat) with the correct weights. See, eg, the Wikipedia article on Voigt profiles for an explanation of where the numbers come from.
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k]/= norm;}
          break;
        case 2: //Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
          norm = 0.;

          if (sigma_zphot_shear(0.,i) == 0.){
            printf("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz =2, but nuisance.sigma_zphot_shear[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
            x2 = (array[1]-zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
      //table[i+1][k] = zdistr_histo_1(zi,(void*)array)*(1.0-outfrac)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1)))) + outfrac;
      //the outliers are currently assumed to have a flat distribution from zmin to zmax.
      // I think the above implementation is not correctly normalized for the outliers. Instead do it at the normalization step. This should probably be updated to allow for a leakage matrix.
            table[i+1][k] = zdistr_histo_1(zi,(void*)array)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))));
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k] = outfrac/da + (1.-outfrac)*table[i+1][k]/norm;}
          break;
        case 3: //Gaussian photo-zs
          norm = 0.;
          if (sigma_zphot_shear(0.,i) == 0.){
            printf("Source galaxy photo-z model underdetermined!\nredshift.shear_photoz =3, but nuisance.sigma_zphot_shear[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
            x2 = (array[1]-zi + bias_zphot_shear(zi,i))/(sqrt(2.)*sigma_zphot_shear(zi,i));
            table[i+1][k] = 0.5*zdistr_histo_1(zi,(void*)array)*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k]/= norm;}
          break;
        case 4: // histogram file contains n(z) estimates for each bin
          array[0] = 1.0*i;
          //call zdistr_histo_n once to update bin boundaries
          zdistr_histo_n(0.,(void*) array);
          norm = int_gsl_integrate_medium_precision(zdistr_histo_n, (void*)array, tomo.shear_zmin[i],tomo.shear_zmax[i],NULL, 1024)*(1.-nuisance.bias_zphot_shear[i]);
          if (norm == 0){
            printf("redshift.c:zdistr_photoz:norm(nz=%d)=0\nEXIT\n",i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            table[i+1][k] = zdistr_histo_n(zi -nuisance.bias_zphot_shear[i],(void*)array)/norm;
          }
          break;        
        default:
          printf("redshift.shear_photoz = %d not supported in this cosmolike version\n",redshift.shear_photoz);
          exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    norm = 0;
    for (i = 0; i <tomo.shear_Nbin; i++){norm += NORM[i];}
    
    for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
      table[0][k] = 0;
      for (i = 0; i <tomo.shear_Nbin; i++){table[0][k]+= table[i+1][k]*NORM[i]/norm;}
    }   
  }
  if (j >= tomo.shear_Nbin){
    printf("redshift.c: zdistr_photoz(z,%d) outside tomo.shear_Nbin range\n", j);
    exit(1);
  }
  if (zz >= zhisto_max || zz < zhisto_min) return 0.0;
  return table[j+1][(int)floor((zz-zhisto_min)/da)];
}



////////////// End shear routines //////////////////////
////////////// End shear routines //////////////////////

double sigma_zphot_clustering (double z, int nz){
  return (1.+z)*nuisance.sigma_zphot_clustering[nz];
}
double bias_zphot_clustering (double z, int nz){
  return (1.+z)*nuisance.bias_zphot_clustering[nz];
}
double pf_histo(double z, void *params) //return pf(z) based on redshift file with one redshift distribution
{
  static double *tab = 0;
  FILE *ein;

  static double zhisto_max,zhisto_min,dz;
  
  if (tab==0){
    double *z_v,space1,space2;
    int i,zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab=create_double_vector(0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein,"%le %le %le %le\n",&z_v[i],&space1,&space2,&tab[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
    }
    fclose(ein);
    dz = (z_v[i-1]-z_v[0])/(1.*i-1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin-1] || zhisto_min > tomo.clustering_zmin[0]){
      printf("Error in redshift.c:pf_histo.c: %s parameters incompatible with tomo.clustering bin choice\nEXIT!\n",redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }
  
  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}

 
double pf_histo_n(double z,  void *params) //return pf(z,j) based on redshift file with structure z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
{
  double *array = (double*)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max,zhisto_min,dz;
  
  if (tab==0){
    double *z_v;
    int i,k,zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab=create_double_matrix(0,tomo.clustering_Nbin-1,0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < tomo.clustering_Nbin; k++){
        fscanf(ein," %le",&tab[k][i]);
      }
    }
    fclose(ein);
    dz =(z_v[i-1]-z_v[0])/(1.*i-1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    // now, set tomography bin boundaries
    for (k = 0; k < tomo.clustering_Nbin; k++){
      double max = tab[k][0];
      for (i = 1; i < zbins; i++){
        if (tab[k][i]> max){max = tab[k][i];}
      }
      i = 0;
      while (tab[k][i] <1.e-8*max && i < zbins-2){i++;}
      tomo.clustering_zmin[k] = z_v[i];
      i = zbins-1;
      while (tab[k][i] <1.e-8*max && i > 0){i--;}
      tomo.clustering_zmax[k] = z_v[i];
      printf("tomo.clustering_zmin[%d] = %.3f,tomo.clustering_zmax[%d] = %.3f\n",k,tomo.clustering_zmin[k],k,tomo.clustering_zmax[k]);
    }

    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin-1] || zhisto_min > tomo.clustering_zmin[0]){
      printf("%e %e   %e %e\n",zhisto_min,tomo.clustering_zmin[0],zhisto_max,tomo.clustering_zmax[tomo.clustering_Nbin-1]);
      printf("Error in redshift.c:pf_histo_n.c: %s parameters incompatible with tomo.clustering bin choice\nEXIT!\n",redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }
  
  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)array[0]][(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}
double pf_photoz(double zz,int j) //returns n(ztrue, j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table = 0;
  static double da = 0.0;
  static double zhisto_max,zhisto_min;
  static int zbins = 2000;
  static nuisancepara N;
  if (table == 0){ //force zmin, zmax to match supplied file
    FILE *ein;
    double *z_v, space;
    int i,k,nz, N = 3;
    if (redshift.clustering_photoz ==4) {N = tomo.clustering_Nbin;}
    nz = line_count(redshift.clustering_REDSHIFT_FILE);
    z_v=create_double_vector(0, nz-1);
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<nz;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < N; k++){fscanf(ein,"%le",&space);}
    }
    fclose(ein);
    redshift.clustering_zdistrpar_zmin = fmax(z_v[0],0.01);
    redshift.clustering_zdistrpar_zmax = z_v[i-1]+(z_v[i-1]-z_v[0])/(1.*i-1.);   
   // printf("redshift.c: pf_photoz: multi-histo case, forcing redshift.clustering_zdistrpar_zmin = %.3f, redshift.clustering_zdistrpar_zmax = %.3f\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
   // for (int i = 0; i < tomo.clustering_Nbin; i++){
   //   tomo.clustering_zmin[i] =redshift.clustering_zdistrpar_zmin;
   //   tomo.clustering_zmax[i] =redshift.clustering_zdistrpar_zmax;
   // }
  }

  if (recompute_zphot_clustering(N) || table==0){
    update_nuisance(&N);
    if (table ==0){
      table   = create_double_matrix(0, tomo.clustering_Nbin, 0, zbins-1);
    }
    zhisto_max =redshift.clustering_zdistrpar_zmax;
    zhisto_min = redshift.clustering_zdistrpar_zmin;
    da = (zhisto_max-zhisto_min)/(1.0*zbins);
    double array[4], NORM[11],zi,norm,x1,x2,eta,outfrac;
    int i,k;
    for (i = 0; i< tomo.clustering_Nbin; i++){
      array[0] =tomo.clustering_zmin[i];
      array[1] =tomo.clustering_zmax[i];
      norm = 0.0;
      outfrac=0.05;
      //eta = 0.68; //this is the approximate mixing parameter for the pseudo voigt for equal FWHM L and G distributions. Note than when using this formulation, the non-equivalent profile re-scalings must be accounted for below.
      eta = 0.736; //this is the approximate mixing parameter for the pseudo voigt for equally rescaled L and G distributions. Note that this leads to a slightly broader profile than when using equal FWHM L and G distributions (at a given gaussian sigma).
      // Currently, the input sigma parameter has not been rescaled. It still gives the width for the Gaussian component of the psuedo Voigt profile. In general, for a given FWHM, the corresponding sigma is smaller for the Voigt profile. For equal FWHM in the component profiles, eta = 0.68 and the effective FWHM is 1.64 that of the Guassian case. For equivalent rescalings, [which is the currently implemented case] eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.
      switch(redshift.clustering_photoz){
        case 0: // no photo-zs, split 'true' n(z) histogram in bins
          norm = int_gsl_integrate_medium_precision(pf_histo, (void*)array, array[0],array[1],NULL, 1024);
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            table[i+1][k] = 0.;
            if (zi >= array[0] && zi <= array[1]){ table[i+1][k] = pf_histo(zi,(void*)array)/norm;}
          }
          break;
          case 1: //Pseudo Voigt (Lorentzian + Gaussian) photo-zs
          norm = 0.;
          if (sigma_zphot_clustering(0.,i) == 0.){
            printf("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz =1, but nuisance.sigma_zphot_clustering[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
            x2 = (array[1]-zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
      table[i+1][k] = pf_histo(zi,(void*)array)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))));
      // this creates a pseudo Voigt profile by adding a Gaussian and Lorentian (convolved with a tophat) with the correct weights. See, eg, the Wikipedia article on Voigt profiles for an explanation of where the numbers come from.
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k]/= norm;}
          break;
        case 2: //Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat) photo-zs
          norm = 0.;

          if (sigma_zphot_clustering(0.,i) == 0.){
            printf("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz =2, but nuisance.sigma_zphot_clustering[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
            x2 = (array[1]-zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
      //table[i+1][k] = pf_histo(zi,(void*)array)*(1.0-outfrac)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1)))) + outfrac;
      //the outliers are currently assumed to have a flat distribution from zmin to zmax.
      // I think the above implementation is not correctly normalized for the outliers. Instead do it at the normalization step. This should probably be updated to allow for a leakage matrix.
            table[i+1][k] = pf_histo(zi,(void*)array)*(eta*0.3183*(atan(x2)-atan(x1)) + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))));
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k] = outfrac/da + (1.-outfrac)*table[i+1][k]/norm;}
          break;
        case 3: //Gaussian photo-zs
          if (sigma_zphot_clustering(0.,i) == 0.){
            printf("Lens galaxy photo-z model underdetermined!\nredshift.clustering_photoz =3, but nuisance.sigma_zphot_clustering[%d] not set\nEXIT!\n", i);
            exit(1);
          }
          norm = 0.;
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
            x1 =(array[0] -zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
            x2 = (array[1]-zi + bias_zphot_clustering(zi,i))/(sqrt(2.)*sigma_zphot_clustering(zi,i));
            table[i+1][k] = 0.5*pf_histo(zi,(void*)array)*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
            norm += table[i+1][k]*da;
          }
          for (k = 0;k<zbins; k++){table[i+1][k]/= norm;}
          break;
        case 4: // histogram file contains n(z) estimates for each bin
          array[0] = 1.0*i;
          pf_histo_n(0.,(void*) array);
          norm = int_gsl_integrate_medium_precision(pf_histo_n, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024);
          if (norm == 0){
            printf("redshift.c:pf_photoz:norm(nz=%d)=0\nEXIT\n",i);
            exit(1);
          }
          for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){ table[i+1][k] = pf_histo_n(zi-nuisance.bias_zphot_clustering[i],(void*)array)/norm;}
          break;
        default:
          printf("redshift.clustering_photoz = %d not supported in this cosmolike version\n",redshift.clustering_photoz);
          exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    norm = 0;
    for (i = 0; i <tomo.clustering_Nbin; i++){norm += NORM[i];}
    for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
      table[0][k] = 0; 
      for (i = 0; i <tomo.clustering_Nbin; i++){table[0][k]+= table[i+1][k]*NORM[i]/norm;}
    } 
  }
/*  if (j >= tomo.clustering_Nbin){
    printf("redshift.c: pf_photoz(z,%d) outside tomo.clustering_Nbin range\n", j);
    exit(1);
  } */
  if (zz >= zhisto_max || zz < zhisto_min) {return 0.0;}
  return table[j+1][(int)floor((zz-zhisto_min)/da)];
}

/*********** routines calculating the number of source and lens galaxies per bin ****************/
double int_nsource(double z, void* param){
  return zdistr_photoz(z,-1);
}
double nsource(int j) //returns n_gal for shear tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table = 0;
  
  double array[3];
  int i;
  if (table ==0|| table[0][0] != survey.n_gal){
    array[0] = zdistr_photoz(0.,0);
    table   = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    if (redshift.shear_photoz == 4){
      printf("redshift.shear_photoz = 4, using tabulated tomo.n_source =");
      for (i = 0; i< tomo.shear_Nbin; i++){
        printf(" %e,", tomo.n_source[i]);
        table[i+1][0] = tomo.n_source[i];
        if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.){
          printf("\n\n!!!!!!!!!\n\nredshift.shear_photoz = 4, tomo.n_source[%d] = %e, EXIT\n\n!!!!!!!!!\n",i,tomo.n_source[i]);
          exit(0);
        }
      }
      printf("\n");
    } else{
      double res,norm;
      norm = int_gsl_integrate_medium_precision(int_nsource, NULL, redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax,NULL, 1024);
      for (i = 0; i< tomo.shear_Nbin; i++){
        res = int_gsl_integrate_medium_precision(int_nsource, (void*)array, tomo.shear_zmin[i],tomo.shear_zmax[i],NULL, 1024);
        table[i+1][0] = res/norm*survey.n_gal;
        if (table[i+1][0] < 0.01||table[i+1][0] > 100){printf("redshift.c:nsource: nsource(%d) = %e\n outside expected range.\nEXIT\n",i,table[i+1][0]);
          exit(1);
        }
      }
    }
    table[0][0] = survey.n_gal;
  }
  return table[j+1][0];
}

double int_nlens(double z, void* param){
  return pf_photoz(z,-1);
}
double nlens(int j) //returns n_gal for clustering tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table =0;
  double array[3];
  int i;
  if (table ==0 ||table[0][0] != survey.n_lens){
    array[0] = pf_photoz(0.,0);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    if (redshift.clustering_photoz == 4){
      printf("redshift.clustering_photoz = 4, using tabulated tomo.n_lens =");
      for (i = 0; i< tomo.clustering_Nbin; i++){
        printf(" %e,", tomo.n_lens[i]);
        table[i+1][0] = tomo.n_lens[i];
              }
      printf("\n");
    } 
    else{
      double res,norm;
      array[0] = tomo.clustering_zmin[0];
      array[1] = tomo.clustering_zmax[tomo.clustering_Nbin-1];
      norm = int_gsl_integrate_medium_precision(int_nlens, NULL, array[0],array[1],NULL, 1024);
      for (i = 0; i< tomo.clustering_Nbin; i++){
        array[0] = tomo.clustering_zmin[i];
        array[1] = tomo.clustering_zmax[i];
        res = int_gsl_integrate_medium_precision(int_nlens, (void*)array, array[0],array[1],NULL, 1024);
        table[i+1][0] = res/norm*survey.n_lens;
        if (table[i+1][0] < 0.001||table[i+1][0] > 100){
          printf("redshift.c:nlens: nlens(%d) = %e\n outside expected range.\nEXIT\n",i,table[i+1][0]);
          exit(1);
        }
      }
    }
    table[0][0] = survey.n_lens;
  }
  return table[j+1][0];
}


double int_for_zmean(double z, void *params){
  double *array = (double*)params;
  return z*pf_photoz(z,(int)array[0]);
}

double zmean(int j){ //mean true redshift of galaxies in tomography bin j
  static double **table = 0;
  if (table ==0){
    double array[1];
    array[0] = pf_photoz(0.,0);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    for (int i = 0; i< tomo.clustering_Nbin; i++){
     array[0]  = 1.0*i;
     table[i][0] = int_gsl_integrate_low_precision(int_for_zmean, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024); 
    }
  }
  return table[j][0];
}

double int_for_zmean_source(double z, void *params){
  
  double *array = (double*)params;
  return z*zdistr_photoz(z,(int)array[0]);
}
double zmean_source(int j){ //mean true redshift of source galaxies in tomography bin j
  static double **table = 0;
  if (table ==0){

    double array[1];
    array[0] = zdistr_photoz(0.,0);
    int i;
    table   = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    for (i = 0; i< tomo.shear_Nbin; i++){
     array[0]  = 1.0*i;
     table[i][0] = int_gsl_integrate_medium_precision(int_for_zmean_source, (void*)array, tomo.shear_zmin[i],tomo.shear_zmax[i],NULL, 1024);
    }
  }
  return table[j][0];
}

double int_for_ggl_efficiency(double z, void *params){
  
  double *array = (double*)params;
  int zs = (int)array[1];
  return pf_photoz(z,(int)array[0])*g_tomo(1./(1.+z),zs)*(1.+z)*f_K(chi(1./(1.+z)));
}
double max_g_tomo(int zs){
  double g, max = 0;
  for (double z = 0.; z < tomo.shear_zmax[zs]; z+=tomo.shear_zmax[zs]/50.){
    g = g_tomo(1./(1.+z),zs)*(1.+z)*f_K(chi(1./(1.+z)));
    if (g>max){max =g;}
  }
  return max;
}

double ggl_efficiency(int zl, int zs){
  static double **table = 0;
  if (table ==0){
    double array[2];
    double init = pf_photoz(0,0);
    init = zdistr_photoz(0,0);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, tomo.shear_Nbin);
    for (int i = 0; i< tomo.clustering_Nbin; i++){
     array[0]  = 1.0*i;
     for (int j = 0; j < tomo.shear_Nbin; j++){ 
        array[1] = 1.0*j;
       table[i][j] = int_gsl_integrate_medium_precision(int_for_ggl_efficiency, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024)/max_g_tomo(j);
      }
    }
  }
  return table[zl][zs];
}
/************ integrands for bin-averaged lens efficiencies **************/

double int_for_g_tomo(double aprime,void *params)
{
  double chi1, chi_prime,val;
  double *ar = (double *) params;
  int zbin= (int) ar[0];
  chi1 = chi(ar[1]);
  chi_prime = chi(aprime);

  val=zdistr_photoz(1./aprime-1.,zbin)*f_K(chi_prime-chi1)/f_K(chi_prime)/(aprime*aprime);
  return val;
}

double g_tomo(double a, int zbin) // for tomography bin zbin
{
  static nuisancepara N;
  static cosmopara C;
  
  static double **table = 0;
  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];
  if (table ==0 || recompute_zphot_shear(N) || recompute_expansion(C)){
    if (table==0) table   = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a-1);
    da = (0.999999-1./(redshift.shear_zdistrpar_zmax+1.))/(Ntable.N_a-1);
    for (j=-1;j<tomo.shear_Nbin;j++) {
      array[0]=(double) j; //if j=-1, no tomography is being done
      aa = 1./(redshift.shear_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a;i++,aa+=da) {
        array[1] = aa;
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),aa,NULL,4000); 
      }      
    }  
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a<=1./(redshift.shear_zdistrpar_zmax+1.) || a>1.0-da) return 0.0;
  return interpol(table[zbin+1], Ntable.N_a, 1./(redshift.shear_zdistrpar_zmax+1.), 0.999999, da, a, 1.0, 1.0); //zbin =-1 is non-tomography
}




double g_cmb (double a){
  static cosmopara C;
  
  static double chi_cmb = 0.;
  double array[1];
  if (recompute_expansion(C)){
    chi_cmb = int_gsl_integrate_medium_precision(int_for_chi,(void*)array, 1./(1091.), 1.,NULL,2000);
    update_cosmopara(&C);
  }
  return f_K(chi_cmb-chi(a))/f_K(chi_cmb);
}

double g_bg (double a, int nzlens){
  printf("g_bg option no longer supported\n'");
  return 0.;
}


/*
double p_zphot_shear(double zp, double zs){ //p(z_{phot} | z_{true}) for source sample
  static double init = -123.;
  static double **table_z = 0;
  static double zmin = 0.01, zmax = 1.99, dz = 0.02;
  double res;
  if (init != 1.0){
    if (!table_z) {table_z = create_double_matrix(0, Ntable.N_zp-1, 0, Ntable.N_zs-1);}
    FILE *ein;
    int i,j;
    double a2,a3,a4,z,z1;
    ein = fopen("../theory/distr_zphot_zspec_100","r");
    EXIT_MISSING_FILE(ein, "p_zphot_shear", "../theory/distr_zphot_zspec_100")
    for (i = 0; i< Ntable.N_zs; i++){
      for (j = 0; j< Ntable.N_zp; j++){
        fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
        table_z[j][i] = a3;
      }
    }
    fclose(ein);
    init = 1.0;
  }
  
  if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0.;}
  res = interpol2d(table_z,Ntable.N_zp, zmin, zmax,dz,zp,Ntable.N_zs, zmin, zmax,dz,zs,1.0,1.0)/dz;
  return res;
}

double p_zspec_shear(double zs, double zp){  //p(z_{true} | z_{phot}) for source sample
  static double init = -123.;
  static double **table_z = 0;
  static double zmin = 0.01, zmax = 1.99, dz = 0.02;
  double res;
  if (init != 1.0){
    if (!table_z) {table_z = create_double_matrix(0, Ntable.N_zs-1, 0, Ntable.N_zp-1);}
    FILE *ein;
    int i,j;
    double a2,a3,a4,z,z1;
    ein = fopen("../theory/distr_zphot_zspec_100","r");
      EXIT_MISSING_FILE(ein, "p_zspec_shear", "../theory/distr_zphot_zspec_100")

    for (i = 0; i< Ntable.N_zs; i++){
      for (j = 0; j< Ntable.N_zp; j++){
        fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
        table_z[i][j] = a4;
      }
      
    }
    fclose(ein);
    init = 1.0;
  }
  if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0.;}
  res = interpol2d(table_z,Ntable.N_zs, zmin, zmax,dz,zs,Ntable.N_zp, zmin, zmax,dz,zp,1.0,1.0)/dz;
  return res;
}

double p_zphot_clustering(double zp, double zs){ //p(z_{phot} | z_{true}) for lens sample
  static double init = -123.;
  static double **table_z = 0;
  static double zmin = 0.01, zmax = 1.99, dz = 0.02;
  double res;
  if (init != 1.0){
    if (!table_z) {table_z = create_double_matrix(0, Ntable.N_zp-1, 0, Ntable.N_zs-1);}
    FILE *ein;
    int i,j;
    double a2,a3,a4,z,z1;
    ein = fopen("../theory/distr_zphot_zspec_100","r");
    EXIT_MISSING_FILE(ein, "p_zphot_clustering", "../theory/distr_zphot_zspec_100")
    for (i = 0; i< Ntable.N_zs; i++){
      for (j = 0; j< Ntable.N_zp; j++){
        fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
        table_z[j][i] = a3;
      }
    }
    fclose(ein);
    init = 1.0;
  }
  
  if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0.;}
  res = interpol2d(table_z,Ntable.N_zp, zmin, zmax,dz,zp,Ntable.N_zs, zmin, zmax,dz,zs,1.0,1.0)/dz;
  return res;
}

double p_zspec_clustering(double zs, double zp){  //p(z_{true} | z_{phot}) for lens sample
  static double init = -123.;
  static double **table_z = 0;
  static double zmin = 0.01, zmax = 1.99, dz = 0.02;
  double res;
  if (init != 1.0){
    if (!table_z) {table_z = create_double_matrix(0, Ntable.N_zs-1, 0, Ntable.N_zp-1);}
    FILE *ein;
    int i,j;
    double a2,a3,a4,z,z1;
    ein = fopen("../theory/distr_zphot_zspec_100","r");
    EXIT_MISSING_FILE(ein, "p_zspec_clustering", "../theory/distr_zphot_zspec_100")
    
    for (i = 0; i< Ntable.N_zs; i++){
      for (j = 0; j< Ntable.N_zp; j++){
        fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
        table_z[i][j] = a4;
      }
      
    }
    fclose(ein);
    init = 1.0;
  }
  if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0.;}
  res = interpol2d(table_z,Ntable.N_zs, zmin, zmax,dz,zs,Ntable.N_zp, zmin, zmax,dz,zp,1.0,1.0)/dz;
  return res;
}
*/