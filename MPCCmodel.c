// Radiation-convection model that integrates the (shortwave and longwave) wavelengths by relevance, weighted by an appropriate factor.
// Version combining linebyline.c (AtW) and line_by_line.c (Dennys Erdtmann) + repwvl_thermal approach.
// Course: Advanced Atmospheric Physics by prof. Mayer

/* Compile with
gcc -Wall -o MPCCmodel.exe MPCCmodel.c ascii.o repwvl_thermal.o repwvl_solar.c -l netcdf -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "ascii.h"
#include "repwvl_thermal.h"
#include "repwvl_solar.h"

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
#define G_RAYLEIGH 0
#define M_AIR 28.97
#define U 1.6605e-27
#define CLOUD_TUNING_PARAM 1.135 // 1.192 Gets multiplied with cloud optical thickness tau_ext_cloud

// prototypes
int integrate_radiation(int nlev, double *T, double dp, double dt, double **tau_thermal, double *wvl_thermal, double *weight_thermal, int nwvl_thermal, double **tau_solar, double *wvl_solar, double *weight_solar, double *E0_solar, int nwvl_solar, double *tau_abs_thermal_cloud, int *cloud_presence, double *tau_ext_solar_cloud, double *omega0_solar_cloud, double *g_solar_cloud);
int thermal_radiative_transfer(double *T, double *Edownthermal, double *Eupthermal, int nlev, int nwvl_thermal, double *wvl_thermal, double *weight_thermal, double **tau_thermal, double *tau_abs_thermal_cloud, int *cloud_presence);
int schwarzschild(double *B, double *tau_abs, double *Eup, double *Edown, int nlev);
int solar_radiative_transfer(double *T, double dp, double *Edownsolar, double *Eupsolar, int nlev, int nwvl_solar, double *wvl_solar, double *weight_solar, double *E0_solar, double **tau_solar, int *cloud_presence, double *tau_ext_solar_cloud, double *omega0_solar_cloud, double *g_solar_cloud);
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
double **memalloc2D(double **matrix, int nrows, int ncols);
void eddington_v2 (double dtau, double g, double omega0, double mu0, double *t, double *r, double *rdir, double *sdir, double *tdir);

int main()
{    
    // First step in determining runtime of the program
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Creating a file to output surface T into
    /*FILE *fptr;
    fptr = fopen("output_surfaceT_MPCCmodel.dat", "w");
    fprintf(fptr, "time (hours) \t \t \t T_surf \n");
    */

    
    // Allocation and initialization
    double p[NLAY];
    double dp = P0 / NLAY; // hPa
    for (int i = 0; i < NLAY; i++)
    {
        p[i] = dp * (i + 0.5);
    }
    
    int status = 0, nlev = 0, nwvl_thermal = 0, nwvl_solar = 0;  
    double **tau_thermal = NULL, **tau_solar;
    double *wvl_thermal = NULL, *weight_thermal = NULL, *wvl_solar=NULL, *weight_solar = NULL,     *E0_solar = NULL;
    double *zlev = NULL, *plev = NULL, *Tlev = NULL, *rhodata = NULL, *plevPa = NULL, *Tlay =     NULL, *rholev = NULL;
    double *H2O_VMR = NULL,  *H2O_VMR_lay=NULL;
    double *CO2_VMR=NULL,  *CO2_VMR_lay=NULL;
    double *O3_VMR=NULL,   *O3_VMR_lay=NULL;
    double *N2O_VMR=NULL,  *N2O_VMR_lay=NULL;
    double *CO_VMR=NULL,   *CO_VMR_lay=NULL;
    double *CH4_VMR=NULL,  *CH4_VMR_lay=NULL;
    double *O2_VMR=NULL,   *O2_VMR_lay=NULL;
    double *HNO3_VMR=NULL, *HNO3_VMR_lay=NULL;
    double *N2_VMR=NULL,   *N2_VMR_lay=NULL;
    
    
    // Read atmospheric profile from test.atm (defined at levels)
    char testfile[FILENAME_MAX] = "./test.atm";
    status = read_9c_file(testfile, &zlev, &plev, &Tlev, &rhodata, &H2O_VMR, &O3_VMR, &CO2_VMR, &CH4_VMR, &N2O_VMR, &nlev);
    if (status != 0)
    {
        fprintf (stderr, "Error %d reading %s\n", status, "test.atm");
        return status;
    }
    if (nlev != NLAY + 1)
    {
        printf("Error: number of layers (nlev) from data does not correspond to set NLAY in model");
    }
    
    
    // Allocate heap memory
    plevPa   = calloc (nlev, sizeof(double));
    O2_VMR   = calloc (nlev, sizeof(double));
    N2_VMR   = calloc (nlev, sizeof(double));
    HNO3_VMR = calloc (nlev, sizeof(double));
    CO_VMR   = calloc (nlev, sizeof(double));

    for (int ilev = 0; ilev < nlev; ilev++)
    {
        //CO2_VMR [ilev] = 280;
        //CH4_VMR[ilev] = 0.6;
        // Convert ppm values to fractions
        H2O_VMR [ilev] *= 1E-6;
        O3_VMR  [ilev] *= 1E-6;
        CO2_VMR [ilev] *= 1E-6; // double CO2
        CH4_VMR [ilev] *= 1E-6;
        N2O_VMR [ilev] *= 1E-6;
    
        //Define (constant i.e. well mixed) profiles for species that are not included in the .atm file
        O2_VMR  [ilev] = 0.2095;
        N2_VMR  [ilev] = 0.7808;
        HNO3_VMR[ilev] = 0.0;
        CO_VMR  [ilev] = 0.0;
        
        // Converting pressure to hPa
        plevPa[ilev] = plev[ilev] * 100.0;
        
        printf("Level %d: Tlev %6.2f CO2_VMR %f CH4_VMR %f \n", ilev, Tlev[ilev], CO2_VMR[ilev], CH4_VMR[ilev]);
    } 
    
    // Set up for solar radiative effects (quantities defined at layers)
    Tlay  = calloc (NLAY, sizeof(double));
    H2O_VMR_lay  = calloc (NLAY, sizeof(double));
    CO2_VMR_lay  = calloc (NLAY, sizeof(double));
    O3_VMR_lay   = calloc (NLAY, sizeof(double));
    N2O_VMR_lay  = calloc (NLAY, sizeof(double));
    CO_VMR_lay   = calloc (NLAY, sizeof(double));
    CH4_VMR_lay  = calloc (NLAY, sizeof(double));
    O2_VMR_lay   = calloc (NLAY, sizeof(double));
    HNO3_VMR_lay = calloc (NLAY, sizeof(double));
    N2_VMR_lay   = calloc (NLAY, sizeof(double));

    for (int ilay=0; ilay<NLAY; ilay++)  {
        Tlay[ilay] = (Tlev[ilay] + Tlev[ilay+1])/2.0;
        H2O_VMR_lay [ilay] = (H2O_VMR[ilay]  + H2O_VMR[ilay+1])/2.0;
        CO2_VMR_lay [ilay] = (CO2_VMR[ilay]  + CO2_VMR[ilay+1])/2.0;
        O3_VMR_lay  [ilay] = (O3_VMR[ilay]   + O3_VMR[ilay+1])/2.0;
        N2O_VMR_lay [ilay] = (N2O_VMR[ilay]  + N2O_VMR[ilay+1])/2.0;
        CO_VMR_lay  [ilay] = (CO_VMR[ilay]   + CO_VMR[ilay+1])/2.0;
        CH4_VMR_lay [ilay] = (CH4_VMR[ilay]  + CH4_VMR[ilay+1])/2.0;
        O2_VMR_lay  [ilay] = (O2_VMR[ilay]   + O2_VMR[ilay+1])/2.0;
        HNO3_VMR_lay[ilay] = (HNO3_VMR[ilay] + HNO3_VMR[ilay+1])/2.0;
        N2_VMR_lay  [ilay] = (N2_VMR[ilay]   + N2_VMR[ilay+1])/2.0;
    }
    
    
    // Read in cloud properties interpolated to 50 repwvl data (thermal as well as solar)
    double *tau_ext_thermal_cloud = NULL,  *tau_ext_solar_cloud = NULL;
    double *omega0_thermal_cloud = NULL,  *omega0_solar_cloud = NULL;
    double *g_thermal_cloud = NULL,  *g_solar_cloud = NULL;
    
    char cloud_thermal_lookup[FILENAME_MAX] = "./ReducedLookupFile_thermal_50wvls_corrected.mie.delta.cldprp";
    status = read_4c_file(cloud_thermal_lookup, &wvl_thermal, &tau_ext_thermal_cloud, &omega0_thermal_cloud, &g_thermal_cloud, &nwvl_thermal);
    if (status != 0)
    {
        fprintf (stderr, "Error %d reading %s\n", status, "cloud thermal lookup file");
        return status;
    }
    if (nlev != NLAY + 1)
    {
        printf("Error: number of layers (nlev) from data does not correspond to set NLAY in model");
    }
    
    char cloud_solar_lookup[FILENAME_MAX] = "./ReducedLookupFile_solar_50wvls.mie.delta.cldprp";
    status = read_4c_file(cloud_solar_lookup, &wvl_solar, &tau_ext_solar_cloud, &omega0_solar_cloud, &g_solar_cloud, &nwvl_solar);
    if (status != 0)
    {
        fprintf (stderr, "Error %d reading %s\n", status, "cloud solar lookup file");
        return status;
    }
    if (nlev != NLAY + 1)
    {
        printf("Error: number of layers (nlev) from data does not correspond to set NLAY in model");
    }
    
    
    //printf("Cloud tuning parameter: %f \n", CLOUD_TUNING_PARAM);
    //printf("Cloud presence: ");

    // Defining layers at which clouds are present and tuning the cloud optical thicknesses
    int cloud_presence[NLAY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0};
    int n_cloud_lyr = 0;
    for (int i = 0; i < NLAY; i++)
    {
        n_cloud_lyr += cloud_presence[i];
        //printf("%d ", cloud_presence[i]);
    }
    //cloud_presence[17] *= 1.05;// lowest clouds: at about 1 km
    //printf("Cloud presence:");
    //show(cloud_presence);
    //exit(0);
    //printf("\n");

    double tau_abs_thermal_cloud[nwvl_thermal];
    for (int w = 0; w < nwvl_thermal; w++) // nwvl_thermal and nwvl_solar are both 50
    {
        tau_abs_thermal_cloud[w] = 0;
        if (n_cloud_lyr > 0){
            tau_ext_thermal_cloud[w] *= CLOUD_TUNING_PARAM/n_cloud_lyr;
            tau_ext_solar_cloud[w] *= CLOUD_TUNING_PARAM/n_cloud_lyr;
            tau_abs_thermal_cloud[w] = tau_ext_thermal_cloud[w] * (1-omega0_thermal_cloud[w]);
        }
    }

    
    // Dynamic time loop
    double progress = 0.;
    double t = 0.0; // seconds
    double oldTlay[NLAY], dTlay[NLAY], absdT[NLAY];
    double dT_limit = 10;   // K
    double dt = 1;      // initialization
    double dt_max = 60*60*6; // maximum time step (6 hours)
    double dt_min = 60*60;
    double t_temp = 0;
    double e_w[NLAY], e_w_old[NLAY], e_w_factor[NLAY];
    
    //fprintf(fptr, "%9.0f \t \t %7.3f \n", t, Tlay[NLAY - 1]); // print initial values
    
    // Initially retrieving the solar tau
    read_tau_solar("./ReducedLookupFile_solar_50wvls.nc", nlev, plevPa, Tlay,
          H2O_VMR_lay, CO2_VMR_lay, O3_VMR_lay, N2O_VMR_lay, CO_VMR_lay, CH4_VMR_lay, O2_VMR_lay, HNO3_VMR_lay, N2_VMR_lay,
          &tau_solar, &wvl_solar, &weight_solar, &E0_solar, &nwvl_solar, 0); // defined at layers

    while (t <= SIMTIME)
    {
        // Retrieving the tau's for thermal and solar:
        read_tau("./ReducedLookupFile_thermal_50wvls_corrected.nc", nlev, plevPa, Tlay, H2O_VMR_lay, CO2_VMR_lay, O3_VMR_lay, N2O_VMR_lay, CO_VMR_lay, CH4_VMR_lay, O2_VMR_lay, HNO3_VMR_lay, N2_VMR_lay, &tau_thermal, &wvl_thermal, &weight_thermal, &nwvl_thermal, 0); // defined at layers  
        if (t_temp>=60*60*5*24)
        {
            read_tau_solar("./ReducedLookupFile_solar_50wvls.nc", nlev, plevPa, Tlay, H2O_VMR_lay, CO2_VMR_lay, O3_VMR_lay, N2O_VMR_lay, CO_VMR_lay, CH4_VMR_lay, O2_VMR_lay, HNO3_VMR_lay, N2_VMR_lay, &tau_solar, &wvl_solar, &weight_solar, &E0_solar, &nwvl_solar, 0); // defined at layers
            t_temp = 0;
        }
        
        // Calculate vapor pressure/partial pressure of H2O
        for (int i = 0; i < NLAY; i++)
            if (t == 0.0) // initial timestep
                e_w_old[i] = 6.1094*exp(17.625*(Tlay[i]-273.15)/(Tlay[i]-273.15+243)); // hPa
            else
                e_w_old[i] = e_w[i];
        
        // Save each initial temperature array for (real) dT calculation later
        for (int i = 0; i < NLAY; i++)
        {
            oldTlay[i] = Tlay[i];
        } 

        // Radiation
        integrate_radiation(nlev, Tlay, dp, dt, tau_thermal, wvl_thermal, weight_thermal, nwvl_thermal, tau_solar, wvl_solar, weight_solar, E0_solar, nwvl_solar, tau_abs_thermal_cloud, cloud_presence, tau_ext_solar_cloud, omega0_solar_cloud, g_solar_cloud);
        
        // Convection
        convection(Tlay, p);

        // Calculate the change in temperature per layer for the timestep and include the water vapor feedback
        for (int i = 0; i < NLAY; i++)
        {
            dTlay[i] = Tlay[i] - oldTlay[i];
            e_w[i] =  6.1094*exp(17.625*(Tlay[i]-273.15)/(Tlay[i]-273.15+243)); // hPa
            e_w_factor[i] = e_w[i]/e_w_old[i];
            H2O_VMR_lay[i] *= e_w_factor[i];         
        }

        // Calculation of dt for the dynamic timestepping
        for (int i = 0; i < NLAY; i++)
        {
            absdT[i] = fabs(dTlay[i]); // Absolute value of dT in order to find max below
        }
        double absdT_max = maximum(absdT, NLAY, 0.0);
        dt = dT_limit / absdT_max;
        dt = dt > dt_max ? dt_max : dt; // Set upper and lower bounds
        dt = dt < dt_min ? dt_min : dt;
        t += dt;
        t_temp += dt;

        // Print into data file
        //fprintf(fptr, "%9.0f \t \t %7.3f \n", t / 3600, Tlay[NLAY - 1]);

        // Condition to end run early (needs readjustment I guess)
        if (absdT_max < 0.0001)
        {
            printf("Ran into break.\n");
            break;
        }

        // Print progress in current timestep
        progress = 100. * (t / SIMTIME);
        printf("Running the time loop... %3.2f  |  dt = %3.2f seconds  |  t = %3.2f days  |  Tsurface = %3.2f K \r", progress, dt, t / 86400, Tlay[NLAY - 1]);
        fflush(stdout);
    }
    //fclose(fptr);

    // Write file with the final layer-temperature structure
    FILE *fptr2 = fopen("output_Tstructure_MPCCmodel_clouds1.192_2CO2_vapor_alf.dat", "w");
    fprintf(fptr2, "Final time = %5.0f hours, T (R) in different layers (L):\n", t / 3600);
    for (int i = 0; i < NLAY; i++)
    {
        fprintf(fptr2, "%d \t \t \t %4.7f\n", i, Tlay[i]);
    }
    fclose(fptr2);

    // Second step in determining runtime of the program
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Runtime in sec %7.2f\n", cpu_time_used);

    // Direct terminal output
    show(Tlay);

    return 0;
}

int integrate_radiation(int nlev, double *T, double dp, double dt, double **tau_thermal, double *wvl_thermal, double *weight_thermal, int nwvl_thermal, double **tau_solar, double *wvl_solar, double *weight_solar, double *E0_solar, int nwvl_solar, double *tau_abs_thermal_cloud, int *cloud_presence, double *tau_ext_solar_cloud, double *omega0_solar_cloud, double *g_solar_cloud)
{
    // Function including radiative integration over wavelength for different molecules and different layers
    // Uses molecule file data for the tau_thermal array, as computed by repwvl
    // I: T array, dp
    // O: new T array

    double Edowntot[nlev], Euptot[nlev], Edownsolar[nlev], Eupsolar[nlev], Edownthermal[nlev], Eupthermal[nlev];

    // Initialization to avoid memory/reference issues
    initialize_double_array(Edowntot, nlev, 0);
    initialize_double_array(Euptot, nlev, 0);
    initialize_double_array(Edownsolar, nlev, 0);
    initialize_double_array(Eupsolar, nlev, 0);
    initialize_double_array(Edownthermal, nlev, 0);
    initialize_double_array(Eupthermal, nlev, 0);

    // Produce irradiances with the input temperature profile
    thermal_radiative_transfer(T, Edownthermal, Eupthermal, nlev, nwvl_thermal, wvl_thermal, weight_thermal, tau_thermal, tau_abs_thermal_cloud, cloud_presence);
    solar_radiative_transfer(T, dp, Edownsolar, Eupsolar, nlev, nwvl_solar, wvl_solar, weight_solar, E0_solar, tau_solar, cloud_presence, tau_ext_solar_cloud, omega0_solar_cloud, g_solar_cloud);
    
    for (int i = 0; i < nlev; i++)
    {
        Edowntot[i] = Edownthermal[i] + Edownsolar[i];
        Euptot[i] = Eupthermal[i] + Eupsolar[i];
    }
    //double E_TOA_net_down = Edowntot[0] - Euptot[0];
    //printf("E_TOA_net_down: %6.2f\n", E_TOA_net_down);
    
    //exit(0);
    
    // New temperature array from irradiances
    irradiancetoT(Edowntot, Euptot, T, dp, dt);
    
    return 0;
}

int thermal_radiative_transfer(double *T, double *Edownthermal, double *Eupthermal, int nlev, int nwvl_thermal, double *wvl_thermal, double *weight_thermal, double **tau_thermal, double *tau_abs_thermal_cloud, int *cloud_presence)
{
    // Integration over wavelengths with Eup and Edown
    // I: T array, wavelength array
    // O: Eupthermal array, Edownthermal array

    double B[NLAY], tau_abs[NLAY];
    double Eup[nlev], Edown[nlev]; // W/m^2

    // Wavelength loop
    for (int w = 0; w < nwvl_thermal; w++)
    {
        initialize_double_array(Eup, nlev, 0);
        initialize_double_array(Edown, nlev, 0);

        for (int i = 0; i < NLAY; i++)
        {
            tau_abs[i] = tau_thermal[w][i] + tau_abs_thermal_cloud[w]*cloud_presence[i];
            // Compute Planck function depending on the wavelength band (=^= sigma*T^4 and cplkavg function)
            B[i] = Planck(T[i], wvl_thermal[w]) * weight_thermal[w];
        }

        // Compute irradiance from the temperature profile
        schwarzschild(B, tau_abs, Eup, Edown, nlev);

        // Integrate over different wavelength contributions
        for (int i = 0; i < nlev; i++)
        {
            Eupthermal[i] += Eup[i];
            Edownthermal[i] += Edown[i];
        }
    }
    /*printf("\n");
    for (int i = 0; i < nlev; i++)
        if (i < 1)
            printf("Level %d: Edownthermal %6.2f \t Eupthermal %6.2f \n", i, Edownthermal[i], Eupthermal[i]);
        else if (i > 19){
            printf("Level %d: Edownthermal %6.2f \t Eupthermal %6.2f \n", i, Edownthermal[i], Eupthermal[i]);
        }
    printf("\n");*/
    

    
    return 0;
}

int schwarzschild(double *B, double *tau_abs, double *Eup, double *Edown, int nlev)
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
            Enet[i] = Edowntot[i] - Euptot[i];
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

int solar_radiative_transfer(double *T, double dp, double *Edownsolar, double *Eupsolar, int nlev, int nwvl_solar, double *wvl_solar, double *weight_solar, double *E0_solar, double **tau_solar, int *cloud_presence, double *tau_ext_solar_cloud, double *omega0_solar_cloud, double *g_solar_cloud)
{
    // Clear the relevant arrays before each use of the function
    double Edir[nlev], Edir_wvl[nlev];
    for(int i = 0; i < nlev; i++){
        Eupsolar[i] = 0.;
        Edownsolar[i] = 0.;
        Edir[i] = 0.;
        Edir_wvl[i] =0.;
     }
    
    // Compute tau_Rayleigh[w]
    double tau_Rayleigh[nwvl_solar];
    double g[nwvl_solar][NLAY];
    double sigma_to_tau = dp*100/(M_AIR*U*G);
    
    for(int w = 0; w < nwvl_solar; w++){
        tau_Rayleigh[w] = sigma_to_tau*(4.840E-31 * (550.0/wvl_solar[w])*(550.0/wvl_solar[w])*(550.0/wvl_solar[w])*(550.0/wvl_solar[w]));
    }

    // Add up tau contributions to tau_tot. This is the point where the imported cloud data is used.
    double tau_tot[nwvl_solar][NLAY]; // --> tau_ext
    double omega0[nwvl_solar][NLAY];  
    
    double tau_abs_solar_cloud[nwvl_solar], tau_sca_solar_cloud[nwvl_solar];
    for (int w = 0; w < nwvl_solar; w++)
    {
        tau_abs_solar_cloud[w] = tau_ext_solar_cloud[w]*(1-omega0_solar_cloud[w]);
        tau_sca_solar_cloud[w] = tau_ext_solar_cloud[w]*omega0_solar_cloud[w];
    }
       
    for(int i = 0; i < NLAY; i++){
        for(int w = 0; w < nwvl_solar; w++)
        {           
            tau_tot[w][i] = tau_Rayleigh[w]+(tau_abs_solar_cloud[w]+tau_sca_solar_cloud[w])*cloud_presence[i]+tau_solar[w][i];
            g[w][i] = (tau_Rayleigh[w]*G_RAYLEIGH + tau_sca_solar_cloud[w]*cloud_presence[i]*g_solar_cloud[w])/(tau_Rayleigh[w] + tau_sca_solar_cloud[w]*cloud_presence[i]);
            omega0[w][i] = (tau_Rayleigh[w] + tau_sca_solar_cloud[w]*cloud_presence[i])/tau_tot[w][i];
            
            if (tau_tot[w][i] > 100)
                tau_tot[w][i] = 100.;
        }
    }
    
    // In order to compute the eddington function for the adding routine, all necessary layer quantites have to be defined
    double t[NLAY], r[NLAY], rdir[NLAY], sdir[NLAY], tdir[NLAY];
    double R[NLAY], Tr[NLAY], Sdir[NLAY], Tdir[NLAY];

    double mu0 = 0.5; // 60 degree solar zenith angle

    for(int w = 0; w < nwvl_solar; w++){
        for(int i = 0; i < NLAY; i++){
            // This loop over wavelengths will be relevant for a while. E.g. all quantites with the suffix _wvl are relevant for one wavelength and will be added in the end of the loop.
            eddington_v2 (tau_tot[w][i], g[w][i], omega0[w][i], mu0, &t[i],
                          &r[i], &rdir[i], &sdir[i], &tdir[i]);
        }
        
        // Define boundary conditions and "initial values"
        Edir_wvl[0] = E0_solar[w]*mu0/2.; 
        R[0] = r[0]; 
        Tr[0] = t[0];  
        Tdir[0] = tdir[0]; 
        Sdir[0] = sdir[0]; 
        
        // Compute R, Tr (Tr instead of T in order to avoid confusion with temperature T) and Tdir
        for (int i = 0; i < NLAY-1; i++)
        {
            R[i+1] = r[i+1] + R[i]*t[i+1]*t[i+1]/(1-R[i]*r[i+1]); 
            Tr[i+1] = Tr[i]*t[i+1]/(1-R[i]*r[i+1]); 
            Tdir[i+1] = Tdir[i] * tdir[i+1]; 
        }
        
        // Using T, compute direct downward irradiance 
        for (int i = 1; i < nlev; i++){
            Edir_wvl[i] = Tdir[i-1] * Edir_wvl[0];
        }
        
        // Now the remaining components of Sdir can be computed
        for (int i = 0; i < NLAY-1; i++){
            Sdir[i+1] = (t[i+1]*Sdir[i]+Tdir[i]*rdir[i+1]*R[i]*t[i+1])/(1-R[i]*r[i+1])+Tdir[i]*sdir[i+1]; 
        }
        
        
        // Albedo calc
        double T_ice = 258;
        double T_no_ice = 293;
        double albedo_ice = 0.7;
        double albedo_no_ice = 0.1;
    
        double Ag = 0.12;
        
        if (T[NLAY-1] > T_no_ice)
            Ag = albedo_no_ice;
        else if (T[NLAY-1] < T_ice)
            Ag = albedo_ice;
        else {
            double x = 1 + (T[NLAY-1] - T_no_ice)/(T_no_ice - T_ice);
            Ag = albedo_ice + (albedo_no_ice - albedo_ice)*((1+0.477/2*x - 0.477/2*pow(x,3)));
        }
        
        double Edownsolar_wvl[nlev], Eupsolar_wvl[nlev];
        
        Edownsolar_wvl[0] = 0.;
        Edownsolar_wvl[NLAY] = Edir_wvl[0]*(Sdir[NLAY-1]+Tdir[NLAY-1]*R[NLAY-1]*Ag)/(1-R[NLAY-1]*Ag); 
        Eupsolar_wvl[NLAY] = Ag*(Edownsolar_wvl[NLAY]+Tdir[NLAY-1]*Edir_wvl[0]); 

        // Calculate the remaining components of Eupsolar and Edownsolar in current wavelength, going from bottom to TOA
        for (int i = NLAY; i>1; i--){
            Edownsolar_wvl[i-1]=(R[i-2]*t[i-1]*Eupsolar_wvl[i]+Edir_wvl[0]*Sdir[i-2]+Edir_wvl[i-1]*rdir[i-1]*R[i-2])/(1-R[i-2]*r[i-1]);
            Eupsolar_wvl[i-1]=(t[i-1]*Eupsolar_wvl[i]+Edir_wvl[0]*Sdir[i-2]*r[i-1]+Edir_wvl[i-1]*rdir[i-1])/(1-R[i-2]*r[i-1]);
        }
        Eupsolar_wvl[0] = t[0]*Eupsolar_wvl[1]; //technically + r[0]*Edownsolar_wvl[0] but second term is set to zero
        
        // Adding up _wvl contributions to the returned irradiances
        for (int i = 0; i < nlev; i++){
            Edir[i] += Edir_wvl[i]*weight_solar[w];
            Eupsolar[i] += Eupsolar_wvl[i]*weight_solar[w];
            Edownsolar[i] += (Edownsolar_wvl[i])*weight_solar[w];
        }
    } // end of wavelength loop

    for (int i=0; i<nlev; i++){
        Edownsolar[i] = Edownsolar[i] + Edir[i];
    }
    
    /*for (int i = 0; i < nlev; i++)
        if (i < 1)
            printf("Level %d: Edownsolar %6.2f \t Eupsolar %6.2f \n", i, Edownsolar[i], Eupsolar[i]);
        else if (i > 19){
            printf("Level %d: Edownsolar %6.2f \t Eupsolar %6.2f \n", i, Edownsolar[i], Eupsolar[i]);
        }
    //printf("\n");
        
    //for (int i = 0; i<nlev; i++)
      //  printf("Level %d: Edownsolar %6.2f \t Eupsolar %6.2f \n", i, Edownsolar[i], Eupsolar[i]);*/
      
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
    for (size_t i = 0; i < length; i++)
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

double **memalloc2D(double **matrix, int nrows, int ncols)
{
    // Function allocating memory for matrix and initialization to zero
    // (possible use: later for tau_total)
    // O: matrix[nrows][ncols] (should be)
    
    matrix = (double **)calloc(nrows, sizeof(double *));
    
    if (!matrix)
        return NULL;

    for (int i = 0; i < nrows; ++i)
    {
        matrix[i] = (double *)calloc(ncols, sizeof(double));
        if (!matrix[i])
        {
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}


/* calculate layer properties t, r, rdir, sdir, and tdir from            */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_v2 (double dtau, double g, double omega0, double mu0,
		   double *t, double *r, double *rdir, double *sdir, double *tdir)
{
  double alpha1=0, alpha2=0, alpha3=0, alpha4=0, alpha5=0, alpha6=0;
  double a11=0, a12=0, a13=0, a23=0, a33=0;
  double lambda=0, b=0, A=0;
  double denom=0;

  /* first, avoid omega0=1 because of instability */
  if (omega0 > 0.999999)
    omega0=0.999999;

  alpha1= (1.0-omega0)+0.75*(1.0-omega0*g);
  alpha2=-(1.0-omega0)+0.75*(1.0-omega0*g);
  
  lambda=sqrt(alpha1*alpha1-alpha2*alpha2);
  
  A=1.0/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));
  
  a11=A*2.0*lambda/alpha2;
  a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));
  
  b=0.5-0.75*g*mu0;
  alpha3=-omega0*b; 
  alpha4=omega0*(1-b);
  
  denom = (1.0/mu0/mu0-lambda*lambda);
  alpha5=((alpha1-1.0/mu0)*alpha3-alpha2*alpha4)/denom;
  alpha6=(alpha2*alpha3-(alpha1+1.0/mu0)*alpha4)/denom;
  
  a33=exp(-dtau/mu0);
  
  a13=alpha5*(1.0-(a11)*(a33))-alpha6*(a12);
  a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11));

  *t    = a11;
  *r    = a12;
  *tdir = a33;
  *rdir = a13 / mu0;
  *sdir = a23 / mu0;
}
