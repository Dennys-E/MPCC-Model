// Radiation-convection model that integrates the wavelengths linebyline.
// Version combining linebyline.c (AtW) and line_by_line.c (Dennys Erdtmann).
// Output is consistent with linebyline.c output.
// Course: Advanced Atmospheric Physics by prof. Mayer
// Amber te Winkel, Jan 2022

/* Comments on compilation
You may need to provide the directory containing ascii.h on the compiler command line:
gcc -c myprogram.c -I.
You also need to link with ascii.o. For that purpose, compile the ascii library, 
gcc -c ascii.c
and link your program with ascii.o:
gcc -o myprogram myprogram.o ascii.o
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "ascii.h"

#define KAPPA 2.0 / 7.0  // (cp-cv)/cp
#define P0 1000.0        // hPa or mbar
#define SIMTIME 3.1536e7 // simulation time in seconds (1 year = 3.1536e7)
#define G 9.81           // m/s^2
#define CP 1004.0        // heat capacity J/(kg*K)
#define NLAY 20
#define SIGMA 5.6701e-8 // W/(m^2 K^4)
#define NMU 20
#define TCONSTANT 288.0     // K
#define EABS 238.0          // W/m^2
#define H 6.62607015e-34    // Planck constant (J/Hz)
#define C 299792458.0       // speed of light
#define KB 1.3806485279e-23 // Boltzmann constant (J/K)
#define NWVL 100001         // Number of wavelengths

// prototypes
int linebyline(double *T, double dp, double dt, double **tauCO2, double *wvl, double *dwvl, int nwvl);
int radiative_transfer(double *T, double *Eintup, double *Eintdown, int nlev, int nwvl, double *wvl, double *dwvl, double **tauCO2);
int schwarzschild(double *T, double *B, double *tau_abs, double *Eup, double *Edown, int nlev);
int irradiancetoT(double *Edowntot, double *Euptot, double *T, double dp, double dt);
int convection(double *T, double *p);
double TtoTheta(double T, double p);
int descending_function(const void *a, const void *b);
double ThetatoT(double Theta, double p);
double maximum(double *array, int len, double lower_bound);
int show(const void *v);
void initialize_double_array(double *array, size_t length, double value);
double Planckfunction(double T, double w);
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
    fptr = fopen("output_surfaceT_lblcombined.dat", "w");
    fprintf(fptr, "time \t \t \t T_surf \n");

    // Allocation and initialization
    double T[NLAY] = {166.0, 177.0, 180.0, 182.0, 188.0, 199.0, 208.0, 217.0, 225.0, 232.0, 238.0, 244.0, 250.0, 255.0, 260.0, 265.0, 270.0, 274.0, 278.0, 283.0};
    double p[NLAY];
    double dp = P0 / NLAY; // hPa
    for (int i = 0; i < NLAY; i++)
    {
        p[i] = dp * (i + 0.5);
    }

    // Importing molecular optical thickness data and putting it all in tauCO2 (2D array)
    int nwvl = 0;
    int nlay = 0;
    double *wvl = NULL;
    double **tauCO2 = NULL;
    double **tauH2O = NULL;
    double **tauO3 = NULL;
    double **tauN2O = NULL;
    double **tauCH4 = NULL;

    char tauCO2file[FILENAME_MAX] = "./lbl.co2.asc";
    char tauH2Ofile[FILENAME_MAX] = "./lbl.h2o.asc";
    char tauO3file[FILENAME_MAX] = "./lbl.o3.asc";
    char tauN2Ofile[FILENAME_MAX] = "./lbl.n2o.asc";
    char tauCH4file[FILENAME_MAX] = "./lbl.ch4.asc";

    ASCII_file2xy2D(tauCO2file, &nwvl, &nlay, &wvl, &tauCO2);
    ASCII_file2xy2D(tauH2Ofile, &nwvl, &nlay, &wvl, &tauH2O);
    ASCII_file2xy2D(tauO3file, &nwvl, &nlay, &wvl, &tauO3);
    ASCII_file2xy2D(tauN2Ofile, &nwvl, &nlay, &wvl, &tauN2O);
    ASCII_file2xy2D(tauCH4file, &nwvl, &nlay, &wvl, &tauCH4);

    //Values are correctly taken from the data, check with:
    /*
    printf("wvl[0]: %.5f, tauCO2[0][0]: %.9f\n", wvl[0], tauCO2[0][0]);
    printf("wvl[0]: %.5f, tauH2O[0][0]: %.9f\n", wvl[0], tauH2O[0][0]);
    printf("wvl[0]: %.5f, tauO3[0][0]: %.9f\n", wvl[0], tauO3[0][0]);
    printf("wvl[0]: %.5f, tauN2O[0][0]: %.9f\n", wvl[0], tauN2O[0][0]);
    printf("wvl[0]: %.5f, tauCH4[0][0]: %.9f\n", wvl[0], tauCH4[0][0]);
    printf("nlay: %d, nwvl: %d\n", nlay, nwvl);
    */

    // Checks:
    if (nwvl != NWVL)
    {
        printf("Number of wavelengths in data does not correspond to the expected value.");
        exit(1);
    }
    if (nlay != NLAY)
    {
        printf("Number of layers in data does not correspond to the expected value.");
        exit(1);
    }

    double dwvl[nwvl];
    for (int w = 0; w < nwvl; w++)
    {
        // Add tau values together for the different molecules (store in tauCO2 for data storage efficiency)
        for (int i = 0; i < NLAY; i++)
        {
            tauCO2[w][i] = tauCO2[w][i] + tauH2O[w][i] + tauO3[w][i] + tauN2O[w][i] + tauCH4[w][i];
        }
        // Calculate the delta in wavelengths (dwvl) for later use
        if (w == 0)
        {
            dwvl[w] = (wvl[1] - wvl[0]) * 1e-9;
        }
        else if (w == nwvl - 1)
        {
            dwvl[w] = (wvl[nwvl - 1] - wvl[nwvl - 2]) * 1e-9;
        }
        else
        {
            dwvl[w] = (wvl[w + 1] - wvl[w - 1]) * 1e-9 / 2.0;
        }
    }

    double absdT[NLAY];
    double absdT_max = 0;
    double dT_limit = 10;   // K
    double dt = 43200;      // initialization (12 hours)
    double dt_max = 172800; // (48 hours) maximum time step in s
    double dt_min = 3600;

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
        linebyline(T, dp, dt, tauCO2, wvl, dwvl, nwvl);

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

        // Condition to end run early
        if (absdT_max < 0.01)
        {
            printf("Ran into break.\n");
            break;
        }

        //Print progress in current timestep
        progress = 100. * (t / SIMTIME);
        printf("Running the time loop... %3.2f  |  dt = %3.2f seconds  |  t = %3.2f days  |  Tsurface = %3.2f K \r", progress, dt, t / 86400, T[NLAY - 1]);
        fflush(stdout);
    }
    fclose(fptr);

    // Write file with the final layer-temperature structure
    FILE *fptr2;
    fptr2 = fopen("output_Tstructure_lblcombined.dat", "w");
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

int linebyline(double *T, double dp, double dt, double **tauCO2, double *wvl, double *dwvl, int nwvl)
{
    // Function including radiative line by line integration over wavelength for different molecules and different layers
    // Called linebyline due to similar previous versions for the window and grey atmosphere.
    // Uses molecule file data for the tau array
    // I: T array, dp
    // O: new T array

    int nlev = NLAY + 1;
    double Edowntot[nlev], Euptot[nlev];

    // Initialization to avoid memory/reference issues
    initialize_double_array(Edowntot, nlev, 0);
    initialize_double_array(Euptot, nlev, 0);

    // Produce irradiances with the input temperature profile
    radiative_transfer(T, Euptot, Edowntot, nlev, nwvl, wvl, dwvl, tauCO2);

    // New temperature array from irradiances
    irradiancetoT(Edowntot, Euptot, T, dp, dt);
    return 0;
}

int radiative_transfer(double *T, double *Eintup, double *Eintdown, int nlev, int nwvl, double *wvl, double *dwvl, double **tauCO2)
{
    // Integration over wavelengths with Eup and Edown
    // I: T array, wavelength array
    // O: Eintup array, Eintdown array

    double B[NLAY], tau_abs[NLAY];
    double Eup[nlev], Edown[nlev]; // W/m^2

    // Wavelength loop
    for (int w = 0; w < nwvl; w++)
    {
        initialize_double_array(Eup, nlev, 0);
        initialize_double_array(Edown, nlev, 0);

        for (int i = 0; i < NLAY; i++)
        {
            tau_abs[i] = tauCO2[w][i];
            // Compute Planck function depending on the wavelength band (=^= sigma*T^4 and cplkavg function)
            B[i] = Planckfunction(T[i], wvl[w]) * dwvl[w];
        }

        // Compute irradiance from the temperature profile
        schwarzschild(T, B, tau_abs, Eup, Edown, nlev);

        // Integrate over different wavelength contributions
        for (int i = 0; i < nlev; i++)
        {
            Eintup[i] += Eup[i];
            Eintdown[i] += Edown[i];
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

double Planckfunction(double T, double w) // T in K, w in nm
{
    // Function that computes the Planck contribution for a wavelength and temperature input.
    // Replaces both SIGMA*pow(T, 4) and cplkavg function

    double wvlmeter = w * 1e-9;
    double res = ((2.0 * H * C * C) / (wvlmeter * wvlmeter * wvlmeter * wvlmeter * wvlmeter)) * 1 / (exp((H * C) / (wvlmeter * KB * T)) - 1.0);
    return res;
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