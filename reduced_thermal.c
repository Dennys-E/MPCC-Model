// Radiation-convection model that integrates the (longwave )wavelengths by relevance, weighted by an appropriate factor.
// Version combining linebyline.c (AtW) and line_by_line.c (Dennys Erdtmann) + repwvl_thermal approach.
// Course: Advanced Atmospheric Physics by prof. Mayer

/* Comments on compilation
You may need to provide the directory containing ascii.h on the compiler command line:
gcc -c myprogram.c -I.
You also need to link with ascii.o. For that purpose, compile the ascii library, 
gcc -c ascii.c
and link your program with ascii.o as well as repwvl_thermal package and netcdf:
gcc -o myprogram myprogram.o ascii.o repwvl_thermal.o -l netcdf -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "ascii.h"
#include "repwvl_thermal.h"

#define KAPPA 2.0 / 7.0  // (cp-cv)/cp
#define P0 1000.0        // hPa or mbar
#define SIMTIME 3.1536e7 // simulation time in seconds (1 year = 3.1536e7)
#define G 9.8065           // m/s^2
#define CP 1004.0        // heat capacity J/(kg*K)
#define NLAY 20
#define SIGMA 5.6701e-8 // W/(m^2 K^4)
#define NMU 20
#define TCONSTANT 288.0     // K
#define EABS 238.0          // W/m^2
#define H_P 6.62607015e-34    // Planck constant (J/Hz)
#define C_LIGHT 299792458.0       // speed of light
#define K_B 1.3806485279e-23 // Boltzmann constant (J/K)

// prototypes
int integrate_radiation(double *T, double dp, double dt, double **tauCO2, double *wvl, double *weight, int nwvl);
int radiative_transfer(double *T, double *Eupthermal, double *Edownthermal, int nlev, int nwvl, double *wvl, double *weight, double **tauCO2);
int schwarzschild(double *T, double *B, double *tau_abs, double *Eup, double *Edown, int nlev);
int irradiancetoT(double *Edowntot, double *Euptot, double *T, double dp, double dt);
int convection(double *T, double *p);
double TtoTheta(double T, double p);
int descending_function(const void *a, const void *b);
double ThetatoT(double Theta, double p);
double maximum(double *array, int len, double lower_bound);
int show(const void *v);
void initialize_double_array(double *array, size_t length, double value);
double Planck(double T, double w);
void write_file(char fname[], int len, double *p, double *var);
void write_Tsurface(double time, double Tsurface);
int memalloc2D(double **matrix, int nrows, int ncols);

int main()
{
    // First step in determining runtime of the program
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Creating a file to output surface T into
    FILE *fptr;
    fptr = fopen("output_surfaceT_reduced_thermal.dat", "w");
    fprintf(fptr, "time \t \t \t T_surf \n");

    // Allocation and initialization
    /*
    double T[NLAY] = {166.0, 177.0, 180.0, 182.0, 188.0, 199.0, 208.0, 217.0, 225.0, 232.0, 238.0, 244.0, 250.0, 255.0, 260.0, 265.0, 270.0, 274.0, 278.0, 283.0}; */
    
    double T[NLAY];
    for(int i=0; i<NLAY; i++)
        T[i]=288.;
    
    double p[NLAY];
    double dp = P0 / NLAY; // hPa
    for (int i = 0; i < NLAY; i++)
    {
        p[i] = dp * (i + 0.5);
    }
    
    int nwvl = 0;
    
    // Initialize wavelength array
    double *wvl = NULL;
    
    // assign pointers to species profiles
    double *zdata = NULL;
    double *pdata = NULL;
    double *Tdata = NULL;
    double *rhodata = NULL;
    double *H2O_VMR = NULL;
    double *CO2_VMR = NULL;
    double *O3_VMR = NULL;
    double *N2O_VMR = NULL;
    double *CO_VMR = NULL;
    double *CH4_VMR = NULL;
    double *O2_VMR = NULL;
    double *HNO3_VMR = NULL;
    double *N2_VMR = NULL;
    int nlev = 0; 
    int status = 0;
    double *plev=NULL, *plevPa=NULL, *Tlev=NULL, *rholev=NULL;
    
    // Read profiles that are saved in the test.atm file.
    // All data arrays here are defined at levels.
    
    // test.atm contains data atmospheric state
    char testfile[FILENAME_MAX] = "./test.atm";
    status = read_9c_file(testfile, &zdata, &pdata, &Tdata, &rhodata, &H2O_VMR, &O3_VMR, &CO2_VMR, &CH4_VMR, &N2O_VMR, &nlev);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, "test.atm");
      return status;
    }

    // allocate heap memory for species data
    O2_VMR   = calloc (nlev, sizeof(double));
    N2_VMR   = calloc (nlev, sizeof(double));
    HNO3_VMR = calloc (nlev, sizeof(double));
    CO_VMR   = calloc (nlev, sizeof(double));

    plevPa   = calloc (nlev, sizeof(double));
  
    for (int ilev=0; ilev<nlev; ilev++) {
    //Convert ppm values to fractions
      H2O_VMR [ilev] *= 1E-6;
      O3_VMR  [ilev] *= 1E-6;
      CO2_VMR [ilev] *= 1E-6;
      CH4_VMR [ilev] *= 1E-6;
      N2O_VMR [ilev] *= 1E-6;
    
      //Define (constant i.e. well mixed) profiles for species that are not included in the .atm file
      O2_VMR  [ilev] = 0.2095;
      N2_VMR  [ilev] = 0.7808;
      HNO3_VMR[ilev] = 0.0;
      CO_VMR  [ilev] = 0.0;

      //Converting pressure to hPa
      plevPa[ilev] = pdata[ilev] * 100.0;
    } 
    
    // tau-->optical thickness profile to be calculated. Weights assign "importance" to respective wvl
    double **tau = NULL;
    double *weight = NULL;

    // determine optical thickness from quantities at levels
    // option: 50 or 100 significant wavelengths, both are saved in lookup files
    read_tau ("./ReducedLookupFile_thermal_50wvls_corrected.nc", nlev, plevPa, Tdata, H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR,CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR, &tau, &wvl, &weight, &nwvl, 1);
    
    double absdT[NLAY];
    double absdT_max = 0;
    double dT_limit = 10;   // K
    double dt = 1;      // initialization (12 hours)
    double dt_max = 60*60*6; // (6 hours) maximum time step in s
    double dt_min = 1;

    // Dynamic time loop
    double progress = 0.;
    double t = 0.0;
    fprintf(fptr, "%9.0f \t \t %7.3f \n", t, T[NLAY - 1]);
    double oldT[NLAY], dT[NLAY]; // working on this for actual difference in T!
    while (t <= SIMTIME)
    {
        // Save each initial temperature array for dT calculation later
        for (int i = 0; i < NLAY; i++)
        {
            oldT[i] = T[i];
        }

        // Radiation
        // Note that tau[][] ans nwvl are the quantities already passed and computed by the repwvl approach
        integrate_radiation(T, dp, dt, tau, wvl, weight, nwvl);

        // Convection
        convection(T, p);

        // Calculate the change in temperature per layer for the timestep
        for (int i = 0; i < NLAY; i++)
        {
            dT[i] = T[i] - oldT[i];
        }

        // Calculation of dt for the dynamic timestepping
        for (int i = 0; i < NLAY; i++)
        {
            absdT[i] = fabs(dT[i]); // Absolute value of dT in order to find max below
        }
        absdT_max = maximum(absdT, NLAY, 0.0);
        dt = dT_limit / absdT_max;
        dt = dt > dt_max ? dt_max : dt; // Set upper and lower bounds
        dt = dt < dt_min ? dt_min : dt;
        t += dt;

        // Print into data file
        fprintf(fptr, "%9.0f \t \t %7.3f \n", t / 3600, T[NLAY - 1]);

        // Condition to end run early (needs readjustment I guess)
        /*
        if (absdT_max < 0.01)
        {
            printf("Ran into break.\n");
            break;
        }
        */

        //Print progress in current timestep
        progress = 100. * (t / SIMTIME);
        printf("Running the time loop... %3.2f  |  dt = %3.2f seconds  |  t = %3.2f days  |  Tsurface = %3.2f K \r", progress, dt, t / 86400, T[NLAY - 1]);
        fflush(stdout);
    }
    fclose(fptr);

    // Write file with the final layer-temperature structure
    FILE *fptr2;
    fptr2 = fopen("output_Tstructure_reduced_thermal.dat", "w");
    fprintf(fptr2, "Final time = %5.0f hours, T (R) in different layers (L):\n", t / 3600);
    for (int i = 0; i < NLAY; i++)
    {
        fprintf(fptr2, "%d \t \t \t %4.7f\n", i, T[i]);
    }
    fclose(fptr2);

    // Second step in determining runtime of the program
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Runtime in sec %7.2f\n", cpu_time_used);

    // Direct terminal output
    show(T);

    return 0;
}

int integrate_radiation(double *T, double dp, double dt, double **tau, double *wvl, double *weight, int nwvl)
{
    // Function including radiative integration over wavelength for different molecules and different layers
    // Uses molecule file data for the tau array, as computed by repwvl
    // I: T array, dp
    // O: new T array

    int nlev = NLAY + 1;
    double Edowntot[nlev], Euptot[nlev];

    // Initialization to avoid memory/reference issues
    initialize_double_array(Edowntot, nlev, 0);
    initialize_double_array(Euptot, nlev, 0);

    // Produce irradiances with the input temperature profile
    radiative_transfer(T, Euptot, Edowntot, nlev, nwvl, wvl, weight, tau);

    // New temperature array from irradiances
    irradiancetoT(Edowntot, Euptot, T, dp, dt);
    return 0;
}

int radiative_transfer(double *T, double *Eupthermal, double *Edownthermal, int nlev, int nwvl, double *wvl, double *weight, double **tau)
{
    // Integration over wavelengths with Eup and Edown
    // I: T array, wavelength array
    // O: Eupthermal array, Edownthermal array

    double B[NLAY], tau_abs[NLAY];
    double Eup[nlev], Edown[nlev]; // W/m^2

    // Wavelength loop
    for (int w = 0; w < nwvl; w++)
    {
        initialize_double_array(Eup, nlev, 0);
        initialize_double_array(Edown, nlev, 0);

        for (int i = 0; i < NLAY; i++)
        {
            tau_abs[i] = tau[w][i];
            // Compute Planck function depending on the wavelength band (=^= sigma*T^4 and cplkavg function)
            B[i] = Planck(T[i], wvl[w]) * weight[w];
        }

        // Compute irradiance from the temperature profile
        schwarzschild(T, B, tau_abs, Eup, Edown, nlev);

        // Integrate over different wavelength contributions
        for (int i = 0; i < nlev; i++)
        {
            Eupthermal[i] += Eup[i];
            Edownthermal[i] += Edown[i];
        }
    }
    return 0;
}

int schwarzschild(double *T, double *B, double *tau_abs, double *Eup, double *Edown, int nlev)
{
    // Function computing irradiance from radiance, independent of wavelengths
    // I: T array, B array, tau_abs array
    // O: Eup array, Edown array

    double alpha[NLAY];
    double mu[NMU];
    double dmu = 1.0 / NMU;
    double Lup[nlev], Ldown[nlev];

    // Directional distribution loop
    for (int k = 0; k < NMU; k++)
    {
        mu[k] = dmu * (k + 0.5); // incoming angles

        // Compute emissivity/absorptivity for each layer and angle
        for (int i = 0; i < NLAY; i++)
        {
            alpha[i] = 1.0 - exp(-tau_abs[i] / mu[k]);
        }

        // Compute the up and down radiance for each layer and angle
        Lup[nlev - 1] = B[NLAY - 1]; // assume emissivity from surface = 1; and that it is isotropic
        for (int j = nlev - 2; j >= 0; j--)
        {
            Lup[j] = Lup[j + 1] * (1 - alpha[j]) + alpha[j] * B[j]; // upward radiance passing through level/surface j
        }
        Ldown[0] = 0;
        for (int j = 0; j < nlev - 1; j++)
        {
            Ldown[j + 1] = Ldown[j] * (1 - alpha[j]) + alpha[j] * B[j]; // downward radiance through level/surface j+1
        }

        // Compute the up and down irradiance by adding the different radiance contributions depending on angle
        for (int i = 0; i < nlev; i++)
        {
            Eup[i] += 2 * M_PI * Lup[i] * mu[k] * dmu;
        }
        for (int i = 0; i < nlev; i++)
        {
            Edown[i] += 2 * M_PI * Ldown[i] * mu[k] * dmu;
        }
    }
    return 0;
}

int irradiancetoT(double *Edowntot, double *Euptot, double *T, double dp, double dt)
{
    // Function combining up and down irradiances into a net irradiance per layer, and then calculating the subsequent change of T.
    // I: Edowntot array, Euptot array
    // O: T array, dT array (--> for dynamic timestepping)

    double Enet[NLAY], dTrad[NLAY];

    for (int i = 0; i < NLAY; i++)
    {
        // Net Earth surface layer irradiance downwards put into the first atmospheric layer
        if (i == NLAY - 1)
        {
            Enet[i] = Edowntot[i] + EABS - Euptot[i];
        }
        else
        {
            Enet[i] = Edowntot[i] - Edowntot[i + 1] + Euptot[i + 1] - Euptot[i]; // net incoming irradiance in the layer
        }

        // Irradiance to temperature change
        dTrad[i] = G * Enet[i] * dt / (CP * dp * 100);
        T[i] += dTrad[i];
    }
    return 0;
}

int convection(double *T, double *p)
{
    // Function describing convective transfer of air masses: order in regard to potential temperature.
    // I/O: T array
    double Theta[NLAY];

    for (int i = 0; i < NLAY; i++)
    {
        Theta[i] = TtoTheta(T[i], p[i]);
    }
    qsort(Theta, sizeof(Theta) / sizeof(*Theta), sizeof(*Theta), descending_function);
    for (int i = 0; i < NLAY; i++)
    {
        T[i] = ThetatoT(Theta[i], p[i]);
    }
    return 0;
}

double TtoTheta(double T, double p)
{
    // Function computing potential temperature from temperature
    return T * pow(P0 / p, KAPPA);
}

int descending_function(const void *a, const void *b)
{
    // Function used for the potential temperature sorting
    double *x = (double *)a;
    double *y = (double *)b;
    if (*x > *y)
        return -1;
    else if (*x < *y)
        return 1;
    return 0;
}

double ThetatoT(double Theta, double p)
{
    // Function computing temperature from  temperature
    return Theta * pow(p / P0, KAPPA);
}

double maximum(double *array, int len, double lower_bound)
{
    double temp = lower_bound;
    for (int i = 0; i < len; i++)
    {
        if (array[i] > temp)
        {
            temp = array[i];
        }
    }
    return temp;
}

int show(const void *v)
{
    // Personalized print function
    double *V = (double *)v;
    for (int i = 0; i < NLAY; i++)
    {
        printf("i = %2d : %10.7f\n", i, V[i]);
    }
    return 0;
}

void initialize_double_array(double *array, size_t length, double value)
{
    // Function to initialize an array to a specified constant
    for (int i = 0; i < length; i++)
    {
        array[i] = value;
    }
}

double Planck(double T, double w) // T in K, w in nm
{
    // Function that computes the Planck contribution for a wavelength and temperature input.
    // Replaces both SIGMA*pow(T, 4) and cplkavg function

    w *= 1e-9;
    return (2.0 * H_P * C_LIGHT * C_LIGHT) / (w * w * w * w * w) / (exp((H_P * C_LIGHT) / (w * K_B * T)) - 1.0)/1.0e9; // W/(m2 nm sterad)
}

void write_file(char fname[], int len, double *p, double *var)
{
    FILE *fp;
    fp = fopen(fname, "w");
    for (int i = 0; i < len; i++)
    {
        fprintf(fp, "%6.1f, %6.1f\n", p[i], var[i]);
    }
    fclose(fp);
}

void write_Tsurface(double time, double Tsurface)
{
    FILE *fp;
    fp = fopen("time_Tsurface.dat", "w");
    fprintf(fp, "%6.1f, %12.6f\n", time, Tsurface);
    fclose(fp);
}

int memalloc2D(double **matrix, int nrows, int ncols)
{
    // Function allocating memory for matrix and initialization to zero
    // (possible use: later for tau_total)
    // O: matrix[nrows][ncols] (should be)

    matrix = (double **)calloc(nrows, sizeof(double *));

    for (int i = 0; i < nrows; ++i)
    {
        matrix[i] = (double *)calloc(ncols, sizeof(double));
    }
    return 0;
}