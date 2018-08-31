void set_cosmological_parameters_to_Planck_WP();
void set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
void set_cosmological_parameters_to_Coyote();
void set_cosmological_parameters_to_WMAP_7years_BAO_SN();
void set_cosmological_parameters_to_SATO();
void set_cosmological_parameters_to_OWLS();
void set_cosmological_parameters_to_MICE();
void set_cosmological_parameters_to_MILLENIUM();
void set_cosmological_parameters_to_Tinker();
void set_cosmological_parameters_to_BCC();
void set_cosmological_parameters_to_chincilla();
void set_cosmological_parameters_to_();
void set_cov_parameters_to_();
void set_survey_parameters_to_();
void set_survey_parameters_to_DES_COV_Y1();
void set_survey_parameters_to_DES_Y1();
void set_survey_parameters_to_CFHTLS();
void set_survey_parameters_to_HSC();
void set_survey_parameters_to_DES();
void set_survey_parameters_to_DES_Tully_Fisher();
void set_survey_parameters_to_LSST_Tully_Fisher();
void set_survey_parameters_to_LSST();
void set_survey_parameters_to_WFIRST_Tully_Fisher();
void set_survey_parameters_to_WFIRST();
void set_survey_parameters_to_WFIRST_extended();
void set_survey_parameters_to_Euclid();
void set_survey_parameters_to_SATO();
void set_survey_parameters_to_DES_BCC_noise();
void set_survey_parameters_to_DES_BCC_nonoise();
void set_survey_parameters_to_DES_COV_BCC();

void set_HOD_redmagic();


// covariance parameters

void set_cov_parameters_to_(char *covparamfile, int output)
{

  char line[256];
  int iline=0;

  FILE* input = fopen(covparamfile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[128];

    iline++;
    if(line[0] == '#') continue;

    sscanf(line, "%128s : %128s", name, val);
    if(strcmp(name, "tmin")==0)
    {
      sscanf(val, "%lf", &covparams.tmin);
      covparams.tmin*=constants.arcmin;
      if(output==1)
      {
        printf("tmin %f \n",covparams.tmin);
      }
      continue;
    }
    else if(strcmp(name, "tmax")==0)
    {
      sscanf(val, "%lf", &covparams.tmax);
      covparams.tmax*=constants.arcmin;
      if(output==1)
      {
        printf("tmax %f \n",covparams.tmax);
      }
      continue;
    }
    else if(strcmp(name, "ntheta")==0)
    {
      sscanf(val, "%d", &covparams.ntheta);
      if(output==1)
      {
        printf("ntheta %d \n",covparams.ntheta);
      }
      continue;
    }
    else if(strcmp(name, "ng")==0)
    {
      sscanf(val, "%d", &covparams.ng);
      if(output==1)
      {
        printf("ng %d \n",covparams.ng);
      }
      continue;
    }
    else if(strcmp(name, "outdir")==0)
    {
      sprintf(covparams.outdir,"%s",val);
      if(output==1)
      {
        printf("outdir %s \n",covparams.outdir);
      }
      continue;
    }
    else if(strcmp(name, "filename")==0)
    {
      sprintf(covparams.filename,"%s",val);
      if(output==1)
      {
        printf("filename %s \n",covparams.filename);
      }
      continue;
    }
    else if(strcmp(name, "ss")==0)
    {
      sprintf(covparams.ss,"%s",val);
      if(output==1)
      {
        printf("ss %s \n",covparams.ss);
      }
      continue;
    }
    else if(strcmp(name, "ls")==0)
    {
      sprintf(covparams.ls,"%s",val);
      if(output==1)
      {
        printf("ls %s \n",covparams.ls);
      }
      continue;
    }
    else if(strcmp(name, "ll")==0)
    {
      sprintf(covparams.ll,"%s",val);
      if(output==1)
      {
        printf("ll %s \n",covparams.ll);
      }
      continue;
    }
  }
}

//various cosmological models

void set_cosmological_parameters_to_(char *cosmofile, int output)
{
  char line[256];
  int iline=0;

  FILE* input = fopen(cosmofile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[128];

    iline++;
    if(line[0] == '#') continue;

    sscanf(line, "%64s : %64s", name, val);
    if(strcmp(name, "Omega_m")==0)
    {
      sscanf(val, "%lf", &cosmology.Omega_m);
      if(output==1)
      {
        printf("Omega_m %f \n",cosmology.Omega_m);
      }
      continue;
    }
    else if(strcmp(name, "Omega_v")==0)
    {
      sscanf(val, "%lf", &cosmology.Omega_v);
      if(output==1)
      {
        printf("Omega_v %f \n",cosmology.Omega_v);
      }
      continue;
    }
    else if(strcmp(name, "sigma_8")==0)
    {
      sscanf(val, "%lf", &cosmology.sigma_8);
      if(output==1)
      {
        printf("sigma_8 %f \n",cosmology.sigma_8);
      }
      continue;
    }
    else if(strcmp(name, "n_spec")==0)
    {
      sscanf(val, "%lf", &cosmology.n_spec);
      if(output==1)
      {
        printf("n_spec %f \n",cosmology.n_spec);
      }
      continue;
    }
    else if(strcmp(name, "w0")==0)
    {
      sscanf(val, "%lf", &cosmology.w0);
      if(output==1)
      {
        printf("w0 %f \n",cosmology.w0);
      }
      continue;
    }
    else if(strcmp(name, "wa")==0)
    {
      sscanf(val, "%lf", &cosmology.wa);
      if(output==1)
      {
        printf("wa %f \n",cosmology.wa);
      }
      continue;
    }
    else if(strcmp(name, "omb")==0)
    {
      sscanf(val, "%lf", &cosmology.omb);
      if(output==1)
      {
        printf("omb %f \n",cosmology.omb);
      }
      continue;
    }
    else if(strcmp(name, "h0")==0)
    {
      sscanf(val, "%lf", &cosmology.h0);
      if(output==1)
      {
        printf("h0 %f \n",cosmology.h0);
      }
      continue;
    }
    else if(strcmp(name, "coverH0")==0)
    {
      sscanf(val, "%lf", &cosmology.coverH0);
      if(output==1)
      {
        printf("coverH0 %f \n",cosmology.coverH0);
      }
      continue;
    }
    else if(strcmp(name, "rho_crit")==0)
    {
      sscanf(val, "%lf", &cosmology.rho_crit);
      if(output==1)
      {
        printf("rho_crit %f \n",cosmology.rho_crit);
      }
      continue;
    }
    else if(strcmp(name, "f_NL")==0)
    {
      sscanf(val, "%lf", &cosmology.f_NL);
      if(output==1)
      {
        printf("f_NL %f \n",cosmology.f_NL);
      }
      continue;
    }
    else if(strcmp(name, "pdelta_runmode")==0)
    {
      sprintf(pdeltaparams.runmode,"%s",val);
      if(output==1)
      {
        printf("runmode %s \n",pdeltaparams.runmode);
      }
      continue;
    }
  }
}


/////////// Survey parameters //////////////

void set_survey_parameters_to_(char *surveyfile, int output)
{

  char line[256];
  int iline=0,i;

  FILE* input = fopen(surveyfile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[256];

    iline++;

    if(line[0] == '#') continue;

    sscanf(line, "%64s : %64s", name, val);
    if(strcmp(name, "area")==0)
    {
      sscanf(val, "%lf", &survey.area);
      if(output==1)
      {
        printf("area %f \n",survey.area);
      }
      continue;
    }
    else if(strcmp(name, "m_lim")==0)
    {
      sscanf(val, "%lf", &survey.m_lim);
      if(output==1)
      {
        printf("mlim %f \n",survey.m_lim);
      }
      continue;
    }
    else if(strcmp(name, "name")==0)
    {
      sprintf(survey.name,"%s",val);
      if(output==1)
      {
        printf("sheartomo %s \n",survey.name);
      }
      continue;
    }
    else if(strcmp(name, "Kcorrect_File")==0)
    {
      sprintf(survey.Kcorrect_File,"%s",val);
      if(output==1)
      {
        printf("Kcorrect_File %s \n",survey.Kcorrect_File);
      }
      continue;
    }
    else if(strcmp(name, "shear_REDSHIFT_FILE")==0)
    {
      sprintf(redshift.shear_REDSHIFT_FILE,"%s",val);
      if(output==1)
      {
        printf("shear_REDSHIFT_FILE %s \n",redshift.shear_REDSHIFT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "clustering_REDSHIFT_FILE")==0)
    {
      sprintf(redshift.clustering_REDSHIFT_FILE,"%s",val);
      if(output==1)
      {
        printf("clustering_REDSHIFT_FILE %s \n",redshift.clustering_REDSHIFT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "sourcephotoz")==0)
    {
      sprintf(survey.sourcephotoz,"%s",val);
      if(output==1)
      {
        printf("sourcephotoz %s \n",survey.sourcephotoz);
      }
      continue;
    }
    else if(strcmp(name, "lensphotoz")==0)
    {
      sprintf(survey.lensphotoz,"%s",val);
      if(output==1)
      {
        printf("lensphotoz %s \n",survey.lensphotoz);
      }
      continue;
    }
    else if(strcmp(name, "galsample")==0)
    {
      sprintf(survey.galsample,"%s",val);
      if(output==1)
      {
        printf("galsample %s \n",survey.galsample);
      }
      continue;
    }
    else if(strcmp(name, "source_tomobins")==0)
    {
      sscanf(val, "%d", &tomo.shear_Nbin);
      if(output==1)
      {
        printf("number of source tomo bins %d \n",tomo.shear_Nbin);
      }
      continue;
    }
    else if(strcmp(name, "lens_tomobins")==0)
    {
      sscanf(val, "%d", &tomo.clustering_Nbin);
      tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
      if(output==1)
      {
        printf("number of lens tomo bins %d \n",tomo.clustering_Nbin);
      }
      continue;
    }
    else if(strcmp(name, "sigma_e")==0)
    {
      sscanf(val, "%lf", &survey.sigma_e);
      if(output==1)
      {
        printf("sigmae %f \n",survey.sigma_e);
      }
      continue;
    }
    else if(strcmp(name, "ggl_overlap_cut")==0)
    {
      sscanf(val, "%lf", &survey.ggl_overlap_cut);
      if(output==1)
      {
        printf("ggl_overlap_cut %f \n",survey.ggl_overlap_cut);
      }
      continue;
    }
    else if(strcmp(name, "source_n_gal")==0)
    {
      if(strcmp(survey.sourcephotoz, "multihisto")==0)
      {
        i=0;
        for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
        {
          double var =0;
          if (i<tomo.shear_Nbin)
          {
            sscanf(p, "%lf", &var);
            tomo.n_source[i]=var;
            if(output==1)
            {
              printf("tomo.n_source[%d]=%f \n",i,tomo.n_source[i]);
            }
          }
/*          if (i>0)
          {
            sscanf(p, "%lf", &var);
            tomo.n_source[i-1]=var;
            printf("tomo.n_source[%d]=%f \n",i-1,tomo.n_source[i-1]);
          }*/
          survey.n_gal+=var;
          i++;
        }
      }
      else
      {
        sscanf(val, "%lf", &survey.n_gal);
      }
      if(output==1)
      {
        printf("ngal %f \n",survey.n_gal);
      }
      continue;
    }
    else if(strcmp(name, "lens_n_gal")==0)
    {
      if(strcmp(survey.lensphotoz, "multihisto")==0)
      {
        i=0;
        for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
        {
          double var =0.;
          if (i<tomo.clustering_Nbin)
          {
            sscanf(p, "%lf", &var);
            tomo.n_lens[i]=var;
            if(output==1)
            {
              printf("tomo.n_lens[%d]=%f \n",i,tomo.n_lens[i]);
            }
          }
/*          if (i>0)
          {
            sscanf(p, "%lf", &var);
            tomo.n_lens[i-1]=var;
            printf("tomo.n_lens[%d]=%f \n",i-1,tomo.n_lens[i-1]);
          }*/
          survey.n_lens+=var;
          i++;
        }
      }
      else
      {
        sscanf(val, "%lf", &survey.n_lens);
      }
      if(output==1)
      {
        printf("nlens %f \n",survey.n_lens);
      }
      continue;
    }
    else if(strcmp(name, "lens_tomogbias")==0)
    { 
      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        gbias.b[i]=var;
        if(output==1)
        {
          printf("gbias[%d]=%f \n",i,gbias.b[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "lens_zphot_sigma")==0)
    { 
      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        nuisance.sigma_zphot_clustering[i]=var;
        if(output==1)
        {
          printf("lens sigma[%d]=%f \n",i,nuisance.sigma_zphot_clustering[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "source_zphot_sigma")==0)
    { 

      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        nuisance.sigma_zphot_shear[i]=var;
        if(output==1)
        {
          printf("source sigma[%d]=%f \n",i,nuisance.sigma_zphot_shear[i]);
        }
        i++;
      }
    }
  }

  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;

}

void set_cosmological_parameters_to_Planck_WP()
{
  cosmology.Omega_m   = 0.315;
  cosmology.Omega_v   = 0.685;
  cosmology.sigma_8   = 0.829;
  cosmology.n_spec    = 0.9603;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.04868;
  cosmology.h0=0.673;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to Planck_WP\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP()
{
  cosmology.Omega_m   = 0.3156;
  cosmology.Omega_v   = 0.6844;
  cosmology.sigma_8   = 0.831;
  cosmology.n_spec    = 0.9645;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0491685;
  cosmology.h0=0.6727;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to Planck_15\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_Joe()
{
  cosmology.Omega_m   = 0.3;
  cosmology.Omega_v   = 0.7;
  cosmology.sigma_8   = 0.8;
  cosmology.n_spec    = 0.96;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.044;
  cosmology.h0=0.68;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to Planck_15\n");
  cosmology.f_NL = 0.0;
}


void set_cosmological_parameters_to_CFHTLens_cov()
{
  cosmology.Omega_m   = 0.279;
  cosmology.Omega_v   = 0.721;
  cosmology.sigma_8   = 0.817;
  cosmology.n_spec    = 0.96;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.046;
  cosmology.h0=0.701;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to CFHTLens covariance\n");
  cosmology.f_NL = 0.0;
}


void set_cosmological_parameters_to_Coyote()
{
  cosmology.Omega_m   =0.281;
  cosmology.Omega_v   = 0.719;
  cosmology.sigma_8   = 0.75;
  cosmology.n_spec    = 0.95;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0459;
  cosmology.h0=0.7;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to Coyote\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_BCC()
{
  cosmology.Omega_m   =0.23;
  cosmology.Omega_v   = 0.77;
  cosmology.sigma_8   = 0.83;
  cosmology.n_spec    = 0.96;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.042;
  cosmology.h0=0.72;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to BCC\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_chincilla()
{
  cosmology.Omega_m   =0.286;
  cosmology.Omega_v   = 0.714;
  cosmology.sigma_8   = 0.82;
  cosmology.n_spec    = 0.96;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.05;
  cosmology.h0=0.7;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("#Cosmology set to chincilla\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_MICE()
{
  cosmology.Omega_m   = 0.25;
  cosmology.Omega_v   = 0.75;
  cosmology.sigma_8   = 0.8;
  cosmology.n_spec    = 0.95;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.044;
  cosmology.h0=0.7;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to MICE\n");
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_WMAP_7years_BAO_SN()
{
  cosmology.Omega_m   = 0.272;
  cosmology.Omega_v   = 0.728;
  cosmology.sigma_8   = 0.807;
  cosmology.n_spec    = 0.961;
  
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0454463597;
  cosmology.h0=0.703;
  cosmology.coverH0= 2997.92458; 

  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set to WMAP_7years_BAO_SN\n");
  cosmology.f_NL = 0.0;
}


void set_cosmological_parameters_to_OWLS()
{
  cosmology.Omega_m   = 0.238;
  cosmology.Omega_v   = 0.762;
  cosmology.sigma_8   = 0.74;
  cosmology.n_spec    = 0.951;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0418;
  cosmology.h0=0.73;
  cosmology.coverH0= 2997.92458; 
  cosmology.rho_crit = 7.4775e+21;
  //printf("Cosmology set to OWLS\n");
}


void set_cosmological_parameters_to_SATO()
{
  cosmology.Omega_m   = 0.238;
  cosmology.Omega_v   = 0.762;
  cosmology.sigma_8   = 0.76;
  cosmology.n_spec    = 0.958;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.042;
  cosmology.h0=0.732;
  cosmology.coverH0= 2997.92458; 
  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set  to Sato\n");
}

void set_cosmological_parameters_to_MILLENIUM()
{
  cosmology.Omega_m   = 0.25;
  cosmology.Omega_v   = 0.75;
  cosmology.sigma_8   = 0.9;
  cosmology.n_spec    = 1.0;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.045;
  cosmology.h0=0.73;
  cosmology.coverH0= 2997.92458;
  cosmology.rho_crit = 7.4775e+21;
  printf("Cosmology set  to Sato\n");
}


void set_survey_parameters_to_DES_Y1()
{
  survey.area   = 1000.;
  survey.n_gal   = 6.08;
  survey.sigma_e   = 0.272259672041*sqrt(2.0);
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.0;
  sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
  sprintf(survey.name,"DES_Y1");
}


void set_survey_parameters_to_DES_SV()
{
  survey.area   = 154.;
  survey.n_gal   = 6.08;
  survey.sigma_e   = 0.272259672041*sqrt(2.0);
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.0;
  sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
  sprintf(survey.name,"DES_SV");
}

void set_survey_parameters_to_DES()
{
    survey.area   = 5000.0;
    survey.n_gal   = 10.0;
    survey.sigma_e   = 0.37;
    survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
    survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
    survey.m_lim=24.0;
    sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
    sprintf(survey.name,"DES");
}

//check Chang 2013 for details
void set_survey_parameters_to_LSST()
{
  survey.area   = 18000.0;
  survey.n_gal   = 26.0;
  survey.sigma_e   = 0.37;  
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.name,"LSST");
}


// from HSC white paper Oct 2012 (latest on dec1 2015)
void set_survey_parameters_to_HSC()
{
  survey.area = 1400.0;  // sq. deg.
  survey.n_gal = 20.0;   // gal per sq arcmin
  survey.sigma_e = 0.22*sqrt(2.0); // "0.22 is RMS shear per component"
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor = 1.0/constants.arcmin/constants.arcmin;
  survey.m_lim = 24.5;  // the photo-z errors depend on the choice for this number
  sprintf(survey.name,"HSC");
}

void set_survey_parameters_to_Euclid()
{
  survey.area   = 15000.0;
  survey.n_gal   = 26.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
  sprintf(survey.name,"Euclid");
}

void set_survey_parameters_to_smallULDB()
{
  survey.area   = 1000.0;
  survey.n_gal   = 20.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
  sprintf(survey.name,"Small ULDB");
}

void set_survey_parameters_to_mediumULDB()
{
  survey.area   = 3382.0;
  survey.n_gal   = 20.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
  sprintf(survey.name,"Medium ULDB");
}

void set_survey_parameters_to_largeULDB()
{
  survey.area   = 7993.0;
  survey.n_gal   = 20.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
  sprintf(survey.name,"Large ULDB");
}

void set_survey_parameters_to_CFHTLS()
{
	survey.area   = 154.0;
	survey.n_gal   = 13.3;
	survey.sigma_e   = 0.42;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_WFIRST() // according to WFIRST AFTA report April 2015
{
  survey.area   = 2200.0;
  survey.n_gal   = 45.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=28.0;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
  sprintf(survey.name,"WFIRST");
}

void set_survey_parameters_to_WFIRST_extended()
{
  survey.area   = 10000.0;
  survey.n_gal   = 45.0;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=28.0;
  sprintf(survey.Kcorrect_File,"../../zdistris/k+e.dat");
}

void set_survey_parameters_to_WFIRST_Tully_Fisher()
{
  survey.area   = 2500.0;
  survey.n_gal   = 8.0;
  survey.sigma_e   = 0.03;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_LSST_Tully_Fisher()
{
	survey.area   = 15000.0;
	survey.n_gal   = 1.1;
	survey.sigma_e   = 0.03;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}


void set_survey_parameters_to_DES_BCC_noise()
{
	survey.area   = 10000.0;
	survey.n_gal   = 10.0;
	survey.sigma_e   = 0.7;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_DES_BCC_nonoise()
{
	survey.area   = 10000.0;
	survey.n_gal   = 10.0;
	survey.sigma_e   = 0.0;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}


void set_HOD_redmagic()
{
  like.Rmin_bias = 0.01;// trust HOD down to ~ 10 kpc/size of galaxies
  redm.hod[0] = 12.1;
  redm.hod[1] = 0.40;
  redm.hod[2] = 13.7;
  redm.hod[3] = 11.9;
  redm.hod[4] = 1.0;
  redm.cg =1.0; //galaxy concentration
  redm.fc =0.19; //central occupation
  redm.b_a =0.; //assembly bias
  redm.parameterization =0; //Zehavi parameterization
  
}

