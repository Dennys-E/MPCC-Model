// This is the model for a grey atmosphere without consideration of wavelength
// dependent absorption and emission
// compile as
//gcc -c ascii.c
// gcc -o line_by_line.exe line_by_line.c cplkavg.c ascii.o

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "ascii.h"
#include "cplkavg.h"

//global values

#define _USE_MATH_DEFINES

#define nol 20 //number of layers
#define kappa 2./7.
#define p0 1000.0 //reference pressure at lower layer (ground interface)
#define SimTime (0.5*3.154e7) //Seconds in one year times factor
#define g 9.81 //m/s²
#define c_p 1004 //J/kg*K
#define sigma 5.6701e-8 // W/(m²K**4)
#define Gamma 9.8 //K/km
#define albedo 0.3
#define Solar_Constant 1361. //W/m²


//functions
double TtoTheta(double T, double p){
  return T*pow(p0/p, kappa);
}

double ThetaToT(double Theta, double p){
  return Theta*pow(p/p0, kappa);
}

int show(const void *v_0, int length){
  double *v = (double *) v_0;
  for(int i=0; i<length; i++){
    printf("i = %2d : %12.6f\n", i, v[i]);
  }
}

int compare_function(const void *a, const void *b) {
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x > *y) return -1;
  else if (*x < *y) return 1; return 0;
}

void write_file(char fname[], int len, double *p, double *var)
{
   FILE * fp;
   fp = fopen (fname,"w");
   for(int i = 0; i < len; i++){
       fprintf (fp, "%6.1f, %6.1f\n", p[i], var[i]);
   }
   fclose (fp);
}

void write_T_ground(double time, double T_ground){
  FILE * fp;
  fp = fopen("T_ground.txt", "w");
  fprintf(fp, "%6.1f, %12.6f\n", time, T_ground);
  fclose(fp);
}

double maximum(double *array, int len, double lower_bound){
  double temp = lower_bound;
  for (int i=0; i<len; i++){
    if(array[i]>temp){
      temp = array[i];
    }
  }
  return temp;
}

void transfer_radiation_test(double *L_up, double *L_down,
                        double *E_up, double *E_down, double *B,
                        double *T, double *wavelengths_files,
                        int nwvl_files, int nlyr_files,
                        double **tau_total, double dmu, int nmu){

  // Loop over wavelngths
  for(int iwvl=0; iwvl<nwvl_files; iwvl++){

    double dwvl_lower=0.;
    double dwvl_upper=0.;


    // Calculate lower and upper wavelength steps --> current wvl is centered
    if (iwvl==0){
      dwvl_upper = (wavelengths_files[1] - wavelengths_files[0])/2;
      dwvl_lower = dwvl_upper;
    }
    else if (iwvl==nwvl_files-1){
      dwvl_lower = (wavelengths_files[nwvl_files-1]-wavelengths_files[nwvl_files-2])/2;
      dwvl_upper = dwvl_lower;
    }
    else{
      dwvl_lower = (wavelengths_files[iwvl] - wavelengths_files[iwvl-1])/2;
      dwvl_upper = (wavelengths_files[iwvl+1] - wavelengths_files[iwvl])/2;
    }

    // loop over lyers in order to compute B(T[i])
    for(int i=0; i<nlyr_files; i++){
      // Integrate planck function for current wavelength step and all layers
      //printf("%d %f %f \n", iwvl, wavelengths_files[iwvl]-dwvl_lower, wavelengths_files[iwvl]+dwvl_upper);

      B[i] = cplkavg(wavelengths_files[iwvl]-dwvl_lower,
                     wavelengths_files[iwvl]+dwvl_upper, T[i]);
    }

    for (int mu_step=1; mu_step <= nmu; mu_step++){

      // Converting counting variable to actual mu (at mid point in interval)
      double mu = (mu_step - 0.5)*dmu;

      //Calculate the L_down transfer from toa to ground at each upper interface
        //(at a given angle)
      for(int i=0; i<=nol-1; i++){
        if (i==0){
          L_down[i] = 0;
        }
        else{
          L_down[i] = L_down[i-1]*exp(-tau_total[iwvl][i-1]/mu)
                      +  B[i]*(1. - exp(-tau_total[iwvl][i-1]/mu));
          E_down[i] = E_down[i] + 2*M_PI*L_down[i]*mu*dmu;
        }

      }

      //Calculate L_up transfer from ground to toa
      for(int i=nol-1; i>=0; i--){
        if (i==nol-1){
          L_up[i] = B[i];
          E_up[i] = E_up[i] + 2*M_PI*L_up[i]*mu*dmu;
        }
        else{
          L_up[i] = L_up[i+1]*exp(-tau_total[iwvl][i]/mu)
                      + B[i]*(1. - exp(-tau_total[iwvl][i]/mu));
          E_up[i] = E_up[i] + 2*M_PI*L_up[i]*mu*dmu;
        }
      }
    }
  }
}

// Takes in current information of temperature and params and propagates
// radiation through layers
void transfer_radiation(double *L_up, double *L_down,
                        double *E_up, double *E_down,
                        double *T,
                        double *tau, double dmu, int nmu){

  // Looping through mu
  for (int mu_step=1; mu_step <= nmu; mu_step++){

    // Converting counting variable to actual mu (at mid point in interval)
    double mu = (mu_step - 0.5)*dmu;

    //Calculate the L_down transfer from toa to ground at each upper interface
      //(at a given angle)
    for(int i=0; i<=nol-1; i++){
      if (i==0){
        L_down[i] = 0;
      }
      else{
        L_down[i] = L_down[i-1]*exp(-tau[i-1]/mu)
                      + (T[i-1]*T[i-1]*T[i-1]*T[i-1]*sigma/M_PI)
                      *(1. - exp(-tau[i-1]/mu));
        E_down[i] = E_down[i] + 2*M_PI*L_down[i]*mu*dmu;
      }

    }

    //Calculate L_up transfer from ground to toa
    for(int i=nol-1; i>=0; i--){
      if (i==nol-1){
        L_up[i] = (T[i]*T[i]*T[i]*T[i]*sigma/M_PI);
        E_up[i] = E_up[i] + 2*M_PI*L_up[i]*mu*dmu;
      }
      else{
        L_up[i] = L_up[i+1]*exp(-tau[i]/mu)
                    + (T[i]*T[i]*T[i]*T[i]*sigma/M_PI)*(1. - exp(-tau[i]/mu));
        E_up[i] = E_up[i] + 2*M_PI*L_up[i]*mu*dmu;
      }
    }
  }
}

// Simulation
//------------------------------------------------------------------------------
int main(){
  double T[nol], Theta[nol], p[nol],
         E_down[nol], E_up[nol], E_total[nol],
         L_down[nol], L_up[nol], tau[nol],
         TRC[nol], TRC_save[nol], B[nol];// Define layer structures
  double dp = p0/nol;
  //double dt = 60.;  //Define time step in s
  double dT_max = 5.; //K
  double tau_at = 1.;
  int nmu = 19;
  double dmu = 1./nmu;
  double E_sol = (1-albedo)*(Solar_Constant/4);
  double TRC_max;
  double dt;
  double dt_max = 60.*60.*10.; //maximum time step in s
  double dt_min = 10.;
  double t;
  double progress = 0.;


  //Fill in initial profiles
  //----------------------------------------------------------------------------
  for(int i=0; i<=nol-1; i++){
    T[i] = 288; //T profile
    tau[i] = tau_at/nol;
    p[i] = dp*(i+1./2.); //Moving to midlle of each layer (offset by dp/2)
    Theta[i] = TtoTheta(T[i], p[i]); //Converting to potential temperature
  }
  //----------------------------------------------------------------------------

  // Import gas-specific tau-files
  //----------------------------------------------------------------------------

  printf("%s...", "Importing data");
  int nwvl_files=0, nlyr_files=0;
  int status=0;

  sleep(1);

  printf("%s ", "CO2...");
  char tauCO2filename[FILENAME_MAX]="./lbl.arts/lbl.co2.asc";
  double *wvlCO2=NULL;
  double **tauCO2=NULL;
  ASCII_file2xy2D (tauCO2filename,
                            &nwvl_files, &nlyr_files, &wvlCO2, &tauCO2);

  printf("%s ", "CH4...");
  char tauCH4filename[FILENAME_MAX]="./lbl.arts/lbl.ch4.asc";
  double *wvlCH4=NULL;
  double **tauCH4=NULL;
  ASCII_file2xy2D (tauCH4filename,
                            &nwvl_files, &nlyr_files, &wvlCH4, &tauCH4);

  printf("%s ", "H2O...");
  char tauH2Ofilename[FILENAME_MAX]="./lbl.arts/lbl.h2o.asc";
  double *wvlH2O=NULL;
  double **tauH2O=NULL;
  ASCII_file2xy2D (tauH2Ofilename,
                            &nwvl_files, &nlyr_files, &wvlH2O, &tauH2O);

  printf("%s ", "N2O...");
  char tauN2Ofilename[FILENAME_MAX]="./lbl.arts/lbl.n2o.asc";
  double *wvlN2O=NULL;
  double **tauN2O=NULL;
  ASCII_file2xy2D (tauN2Ofilename,
                            &nwvl_files, &nlyr_files, &wvlN2O, &tauN2O);

  printf("%s ", "O3...");
  char tauO3filename[FILENAME_MAX]="./lbl.arts/lbl.o3.asc";
  double *wvlO3=NULL;
  double **tauO3=NULL;
  ASCII_file2xy2D (tauO3filename,
                            &nwvl_files, &nlyr_files, &wvlO3, &tauO3);


  printf("%s\n\n", "adding...");

  // Allocating memory for and defining wavelengths array in order to avoid
  // confusion. However, should be equal for all gases
  double *wavelengths_files;
  wavelengths_files = (double *) calloc(nwvl_files, sizeof(double));

  for(int iwvl=0; iwvl<nwvl_files; iwvl++){

    wavelengths_files[iwvl] = wvlCO2[iwvl]; //Arbitrary choice

    if(wvlCO2[iwvl]!=wavelengths_files[iwvl]){
      printf("%s\n", "Error: Imported files do not have equal wavl steps!");
      exit(0);
    }
    if(wvlCH4[iwvl]!=wavelengths_files[iwvl]){
      printf("%s\n", "Error: Imported files do not have equal wavl steps!");
      exit(0);
    }
    if(wvlH2O[iwvl]!=wavelengths_files[iwvl]){
      printf("%s\n", "Error: Imported files do not have equal wavl steps!");
      exit(0);
    }
    if(wvlN2O[iwvl]!=wavelengths_files[iwvl]){
      printf("%s\n", "Error: Imported files do not have equal wavl steps!");
      exit(0);
    }
    if(wvlO3[iwvl]!=wavelengths_files[iwvl]){
      printf("%s\n", "Error: Imported files do not have equal wavl steps!");
      exit(0);
    }
  }

  // Allocating memory for tau_total
  double **tau_total;
  tau_total = (double **) calloc(nwvl_files, sizeof(double*));

  for(int i = 0; i < nwvl_files; ++i) {
    tau_total[i] = (double *) calloc(nlyr_files, sizeof(double));
  }

  // Add up gas-specific taus to tau_total[iwvl][ilyr]
  for(int iwvl=0; iwvl<nwvl_files; iwvl++){

    for(int i=0; i<nlyr_files; i++){
      tau_total[iwvl][i] = tauCO2[iwvl][i]+tauCH4[iwvl][i];
                           +tauH2O[iwvl][i]+tauN2O[iwvl][i]+tauO3[iwvl][i];
    }
  }


  printf("%s\n\n", "...done");
  int answer;
  printf("%s\n", "Continue? yes=1/no=2");
  scanf("%d",&answer);

  if(answer==2){
    printf("%s\n", "Stopped");
    exit(0);
  }


  //----------------------------------------------------------------------------

  //Starting to time-step
  t = 0.;
  while(t<=SimTime){

    // Clearing E and L arrays
    for(int i=0; i<=nol-1; i++){
      E_down[i] = 0;
      E_up[i] = 0;
      E_total[i] = 0;
    }

    //Schwarzschild:
    //transfer_radiation(L_up, L_down, E_up, E_down, T, tau, dmu, nmu);

    transfer_radiation_test(L_up, L_down,
                            E_up, E_down, B,
                            T, wavelengths_files,
                            nwvl_files, nlyr_files,
                            tau_total, dmu, nmu);


    //i-loop to compute Temperature Rate of Change (TRC) in each layer
    for(int i=nol-1; i>=0; i--){
      if (i==nol-1){
        E_total[i] = E_down[i] + E_sol - E_up[i];
      }
      else{
        E_total[i] = E_down[i]-E_down[i+1] + E_up[i+1]-E_up[i];
      }
      //Temperature rate of change by absolute value, in order to find max later
      TRC_save[i] = fabs((g*E_total[i])/(100.*c_p*dp));
      //Actual TRC
      TRC[i] = ((g*E_total[i])/(100.*c_p*dp));
    } //i-loop ends here


    //Find maximum absolute value TRC in order to optimize dt
    TRC_max = maximum(TRC_save, nol, 0.);
    dt = dT_max/TRC_max;

    //Set upper and lower bounds
    if (dt >= dt_max){
      dt = dt_max;
    }
    if (dt <= dt_min){
      dt = dt_min;
    }
    if (t+dt>SimTime){
      dt = SimTime-t;
    }

    //Time step adding dT
    for(int i=nol-1; i>=0; i--){
      T[i] = T[i] + TRC[i]*dt;
      Theta[i] = TtoTheta(T[i], p[i]); //Convert to potential temperature
    }

    //Sort Theta in descending order
    qsort(Theta, sizeof(Theta)/sizeof(*Theta),
          sizeof(*Theta), compare_function);

    //Convert Theta back to T
    for(int i=0; i<=nol-1; i++){
      T[i] = ThetaToT(Theta[i], p[i]);
    }

    write_file("T.txt", nol, p, T);
    write_file("Theta.txt", nol, p, Theta);

    t = t+dt;
    //Print progress in current tau step
    progress = 100.*(t/SimTime);
    //printf("%f\r", progress);
    printf("Running... %3.2f  |  dt = %3.2f seconds  |  t = %3.2f days  |  T_ground = %3.2f K \r",
            progress, dt, t/(60*60*24), T[nol-1]);
    fflush(stdout);
  }

  printf("\n...done");

  return 0;
}
