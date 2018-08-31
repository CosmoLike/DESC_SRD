#ifndef FASTPTMODE
#include "pt.h"
double K_CH0(double k_mpch){
  return k_mpch*cosmology.coverH0;
}
double K_MPCH(double k_ch0){
  return k_ch0/cosmology.coverH0;
}
//interpolate pretabulated FPT output, ignoring some cosmology dependence
double PT_d1d2(double k_coverH0){ //interpolate PT_AB[0] - Pd1d2
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
    logkmin =log(1.e-4);
    logkmax = log(1.e+3);
    dlgk = (logkmax-logkmin)/(1.0*N);
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(PT_AB[0], N, logkmin, logkmax, dlgk,lgk, 0.,0.)*pow(cosmology.sigma_8/0.80,4.0);
}

double PT_d2d2(double k_coverH0){ //interpolate PT_AB[1]
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
    logkmin =log(1.e-4);
    logkmax = log(1.e+3);
    dlgk = (logkmax-logkmin)/(1.0*N);
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(PT_AB[1], N, logkmin, logkmax, dlgk,lgk, 0.,0.)*pow(cosmology.sigma_8/0.80,4.0);
}

double PT_d1s2(double k_coverH0){ //interpolate PT_AB[2]
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
    logkmin =log(1.e-4);
    logkmax = log(1.e+3);
    dlgk = (logkmax-logkmin)/(1.0*N);
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(PT_AB[2], N, logkmin, logkmax, dlgk,lgk, 0.,0.)*pow(cosmology.sigma_8/0.80,4.0);
}

double PT_d2s2(double k_coverH0){ //interpolate PT_AB[3]
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
    logkmin =log(1.e-4);
    logkmax = log(1.e+3);
    dlgk = (logkmax-logkmin)/(1.0*N);
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(PT_AB[3], N, logkmin, logkmax, dlgk,lgk, 0.,0.)*pow(cosmology.sigma_8/0.80,4.0);
}

double PT_s2s2(double k_coverH0){ //interpolate PT_AB[4]
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
    logkmin =log(1.e-4);
    logkmax = log(1.e+3);
    dlgk = (logkmax-logkmin)/(1.0*N);
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(PT_AB[4], N, logkmin, logkmax, dlgk,lgk, 0.,0.)*pow(cosmology.sigma_8/0.80,4.0);
}

double PT_sigma4(double k_coverH0){
  static int N = 0;
  // only call FASTPT if cosmology changed since last call
  if (N ==0){
    printf("FASTPT not linked, using pretabulated output at fixed cosmology.\n");
    printf("Compile cosmolike with -DFASTMODE for full calculation.\n\n");
    N = 70;
  }
  return PT_AB[5][0]*pow(cosmology.sigma_8/0.80,4.0);
}

#else

#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void write_plin(char *filebase);
void get_FPT_bias(void);
void get_FPT_bias_command_line(void);
void FPT_bias_interface(void);
double PT_d1d2(double k_coverH0);
double PT_d2d2(double k_coverH0);
double PT_d1s2(double k_coverH0);
double PT_d2s2(double k_coverH0);
double PT_s2s2(double k_coverH0);
double PT_sigma4(double k_coverH0);

void get_FPT_bias(void){
	//FPT_bias_interface();
	//command line version if compiling C and python together doesn't work
	//CREATES LOTS OF I/O - DO NOT USE ON A CLUSTER 
	get_FPT_bias_command_line();
}
double K_CH0(double k_mpch){
  return k_mpch*cosmology.coverH0;
}
double K_MPCH(double k_ch0){
  return k_ch0/cosmology.coverH0;
}
void FPT_input(double k[FPT.N], double P[FPT.N]){
  int i;
  double dlgk,lgk;
  if ((int)log10(FPT.k_max/FPT.k_min)*FPT.N_per_dec != FPT.N){
    printf("pt.c:FPT_input: inconsistent k-range and number of bins for FPT\n");
    printf("FPT.k_min=%e, FPT.k_max=%e, FPT.N_per_dec=%d; FPT.N=%d\nEXIT\n",FPT.k_min,FPT.k_max,FPT.N_per_dec,FPT.N);
    exit(1);
  }
  dlgk = log(10.)/(double) FPT.N_per_dec;
  for (lgk = log(FPT.k_min), i = 0; i< FPT.N; i++,lgk+=dlgk){
    k[i] = exp(lgk);
    P[i] = p_lin(K_CH0(exp(lgk)),1.0)*pow(cosmology.coverH0,3.0);
    //printf("%e %e\n",k[i],P[i]);
 } 
}
void write_plin(char *filebase){
  //create filename from filebas, processid, Omega_m and A_ia value to ensure unique file names
  double k[FPT.N], P[FPT.N];
  FPT_input(k, P);
  sprintf(FPT.Plin_FILE,"%s.%d_%e_%e",filebase,getpid(),cosmology.Omega_m,nuisance.A_ia);
  FILE *F;
  F = fopen(FPT.Plin_FILE,"w");
  fprintf(F,"# k[h/Mpc] P_lin(k,z=1.0) [Mpc^3/h^3]\n");
  for (int i = 0; i< FPT.N; i++){
    fprintf(F,"%e %e\n",k[i],P[i]);
  }
  fclose(F);
}
void get_FPT_bias_command_line(void){
  static cosmopara C; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    // update local cosmology
    update_cosmopara(&C); 
    char command[500];
    FILE *F;
    int i;

    if (FPT.tab_AB ==0){ //if table doesn't exit yet, create it
      FPT.tab_AB = create_double_matrix(0, FPT.N_AB-1, 0, FPT.N-1);
    }

    //write input power spectrum for current cosmology
    write_plin("P_lin");
    //now do FASTPT calculation
    sprintf(command,"python bias_ek.py %s -o",FPT.Plin_FILE);
   // printf("%s\n",command);
    int ret = system(command);
    if (ret){
      printf("like_fourier.c:get_FPT_bias_command_line: FAST-PT call failed\n%s\nEXIT\n",command);
      exit(1);
    }
    //read in FASTPT output
    F = fopen(FPT.Plin_FILE,"r");
    if (!F){
      printf("like_fourier.c:get_FPT_bias_command_line: can't read FPT.Plin_FIlE %s\nEXIT\n",FPT.Plin_FILE);
      exit(1);
    }
    // read in line by line, copy A, B values into FPT.tab_AB
    double Pk_unit_conversion = pow(cosmology.coverH0,3.0);
    for (i =0; i < FPT.N; i++){
      double k,P,A,B,C,D,E,G;
      ret = fscanf(F,"%le %le %le %le %le %le %le %le",&k,&P,&A,&B,&C,&D,&E,&G);
      if (ret != 8){printf("like_fourier.c:get_FPT_bias_AB: formatting error in column %d of %s\nEXIT\n", i,FPT.Plin_FILE);}
      FPT.tab_AB[0][i] = A/Pk_unit_conversion;
      FPT.tab_AB[1][i] = B/Pk_unit_conversion;
      FPT.tab_AB[2][i] = C/Pk_unit_conversion;
      FPT.tab_AB[3][i] = D/Pk_unit_conversion;
      FPT.tab_AB[4][i] = E/Pk_unit_conversion;
      FPT.tab_AB[5][i] = G/Pk_unit_conversion;
    }
    fclose(F);
    remove(FPT.Plin_FILE);
//    for (int i=0; i<FPT.N; i+= 10){
//            printf("%d %e %e\n",i,FPT.tab_AB[0][i],FPT.tab_AB[1][i]);
//    }
  }
}

// Convert a double pointer into a numpy array
PyObject * arrayFromDoublePointer(int n, double *x)
{
    int ndim = 1;
    npy_intp size = n;
    PyObject * array = PyArray_SimpleNewFromData(ndim, &size, NPY_DOUBLE, x);
    return array;
}

void setup_numpy(){
    import_array();
}

void FPT_bias_interface(void){
  static cosmopara C;
  static int init = 1; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    // update local cosmology
    update_cosmopara(&C); 
    if (FPT.tab_AB ==0){ //if table doesn't exit yet, create it
      FPT.tab_AB = create_double_matrix(0, FPT.N_AB-1, 0, FPT.N-1);
    }
    if (init){
      // Start python
      Py_Initialize();

      // Setup numpy - mysterious but needed!
      setup_numpy();
      // add FPT.path to PYTHONPATH - only relevant on zodiac where system python path doesn't work
     /* DIR* dir = opendir(FPT.path);
      if (dir)
      {
          // FPT.path exists. //
        closedir(dir);
        PyObject* sysPath = PySys_GetObject((char*)"path");
        PyObject* FPTPath = PyString_FromString(FPT.path);
        PyList_Append(sysPath, FPTPath);
        Py_DECREF(FPTPath);
      }*/
      init =0;
    }
    // Name of module and function to import and the number of args it has
    char * module_name_c = "interface_cosmolike_v2_1";
    char * function_name_c = "call_FPT_bias";
    int number_arguments = 2;
    // Convert C string to python string.
    PyObject * module_name_py = PyString_FromString(module_name_c);

    // Import the module
    PyObject * module = PyImport_Import(module_name_py);

    // Free module_name_py
    Py_DECREF(module_name_py);
    if (module == NULL) {
        fprintf(stderr,"Error: Module %s not found\n",module_name_c);
        fprintf(stderr,"add cosmolike/FASTPT_2_1 to PYTHONPATH\n");
        PyErr_Print();
        exit(1);
    }
    // Get the function out of the module
    PyObject * function = PyObject_GetAttrString(module, function_name_c);

    // Check that the function was found.
    if (function == NULL) {
        fprintf(stderr,"pt.c:FPT_bias_AB_interface: Error: Function %s not found in module %s\n",
            function_name_c,module_name_c);
        exit(1);
    }
    int N = FPT.N;
    double k[N], P[N];
    FPT_input(k, P);
    // Convert to numpy array.
    PyObject * arrayk = arrayFromDoublePointer(N, k);
    PyObject * arrayP = arrayFromDoublePointer(N, P);
    // Check we managed to create array
    if (arrayk == NULL || arrayP==NULL) {
        fprintf(stderr, "pt.c:FPT_bias_AB_interface: Error: Unable to create PyObject arrays\n");
        exit(1);
    }
    // Create the arguments to be passed into the function
    PyObject * arguments = PyTuple_New(number_arguments);    
    PyTuple_SetItem(arguments, 0, arrayk);
    PyTuple_SetItem(arguments, 1, arrayP);
    PyObject * result = PyObject_CallObject(function, arguments);
    
    // Free the arguments
    Py_DECREF(arguments);

    // Check that the function worked and print the python
    // traceback if it fails
    if (result == NULL) {
        printf("pt.c:FPT_bias_AB_interface: Error: call_FPT_bias failed\n");
        PyErr_Print();
        exit(1);
    }

    //copy A, B values into FPT.tab_AB
    double Pk_unit_conversion = pow(cosmology.coverH0,3.0);
    for (int i =0; i < FPT.N; i++){
      for (int j =0; j < FPT.N_AB; j++){
        double* A = PyArray_GETPTR2(result, j+2,i);
        FPT.tab_AB[j][i] = *A/Pk_unit_conversion;
      }
    }
    Py_DECREF(result);
/*    printf("//FAST_PT terms calculated at sigma_8 = %e\n",cosmology.sigma_8);
    printf("static double PT_AB[%d][%d]={",FPT.N_AB,FPT.N);
    for (int j =0; j < FPT.N_AB; j++){
      printf("{");
     for (int i =0; i < FPT.N-1; i++){
        printf("%.10e, ",FPT.tab_AB[j][i]);
      }
      printf("%.10e},\n",FPT.tab_AB[j][FPT.N-1]);
    }*/
  }
}

double PT_d1d2(double k_coverH0){ //interpolate FPT.tab_AB[0] - Pd1d2
  static cosmopara C;
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  // convert k from coverH0 (used in cosmolike) to h/Mpc units (used in FASTPT)
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[0], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d2d2(double k_coverH0){ //interpolate FPT.tab_AB[1]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[1], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d1s2(double k_coverH0){ //interpolate FPT.tab_AB[2]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    double k_unit_conversion = 1./cosmology.coverH0;
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[2], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d2s2(double k_coverH0){ //interpolate FPT.tab_AB[3]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[3], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_s2s2(double k_coverH0){ //interpolate FPT.tab_AB[4]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(K_MPCH(k_coverH0));
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[4], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_sigma4(double k_coverH0){
  static cosmopara C; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
  }
  double k = K_MPCH(k_coverH0);
  if (k < FPT.k_min || k >= FPT.k_max){return 0.;}
  return FPT.tab_AB[5][0];
}

#endif
