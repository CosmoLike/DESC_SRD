//these needs to be called during init at fiducial cosmology
void init_BAO_DES();
void init_BAO_LSST();
void init_BAO_WFIRST();
void init_BAO_BOSS();
double dist_BAO(double z);


typedef struct{
  int N;
  double z[10];
  double data[10];
  double sigma[10];
}baopara;
baopara BAO;


void init_BAO_DES(){
  int i;
  BAO.N = 2;
  BAO.z[0] = 0.8;
  BAO.z[1] = 1.0;
  for(i = 0; i< BAO.N; i++){
    BAO.data[i] = dist_BAO(BAO.z[i]);
    BAO.sigma[i] = 0.03;
  }
}

void init_BAO_WFIRST(){
  int i;
  BAO.N = 8;
  BAO.z[0] = 1.25;
  BAO.z[1] = 1.5;
  BAO.z[2] = 1.75;
  BAO.z[3] = 2.0;
  BAO.z[4] = 2.25;
  BAO.z[5] = 2.5;
  BAO.z[6] = 2.75;
  BAO.z[7] = 3.0;
  for(i = 0; i< BAO.N; i++){
    BAO.data[i] = dist_BAO(BAO.z[i]);
    BAO.sigma[i] = 0.02;
  }
  printf("LSST BAO initialized\n");
}

void init_BAO_LSST(){
  int i;
  BAO.N = 4;
  BAO.z[0] = 1.25;
  BAO.z[1] = 1.5;
  BAO.z[2] = 1.75;
  BAO.z[3] = 2.0;
  for(i = 0; i< BAO.N; i++){
    BAO.data[i] = dist_BAO(BAO.z[i]);
    BAO.sigma[i] = 0.03;
  }
  printf("LSST BAO initialized\n");
}

void init_BAO_BOSS(){
  int i;
  BAO.N = 1;
  BAO.z[0] = 0.32; 
  BAO.z[1] = 0.57;
  BAO.sigma[0] = 0.02;    
  BAO.sigma[1] = 0.01;
  for(i = 0; i< BAO.N; i++){
    BAO.data[i] = dist_BAO(BAO.z[i]);
  }
}


//added factor of h as BAO measures scale in Mpc, not Mpc/h
double dist_BAO(double z){
  return pow(1.+z,-1.)*f_K(chi(1./(1+z)))*cosmology.h0;
}

