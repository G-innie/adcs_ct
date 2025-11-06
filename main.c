#include <stdlib.h>		// including several standard C header files. The #include directive tells the preprocessor to insert the contents of another file into the source code at the point where the #include directive is found. Include directives are typically used to include the C header files for C functions that are held outsite of the current source file.
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <adcs/adcs.h>		// including ISIS ADCS Next Gen Library header files
#include <adcs/linear_algebra.h>
#include <adcs/adcs_utils.h>
#include <adcs/adcs_determination.h>

#include <simfunc.h> //includes the functions called in the main, note: the definitions are in /headers/simfunc.c
#include <PDtoSSS.h>
#include <DisturbanceTorques.h> //includes the functions needed for the Disturbance Torques computation.
#include "selector.h"
#include "ftcmd.h"

#define ADCSLOOP_iMTQATTEMPTS 2		//In the C Programming Language, the #define directive allows the definition of macros within your source code. These macro definitions allow constant values to be declared for use throughout your code.
#define ADCSLOOP_MAGCORRPROCEDURE 0
#define Imax (750E-6 * 10E9)
#define SSS_BUFFER_SIZE (64)

ftcmd_t* imtq;		//ftcmd_t stantds for functional testing commands (https://kth.diva-portal.org/smash/get/diva2:1703720/FULLTEXT01.pdf)
ftcmd_t* sss;		// A pointer variable, with the name ftcmd_t. A pointer is a variable that stores the memory address of another variable as its value. A pointer variable points to a data type (like int) of the same type, and is created with the * operator.

char sss_cmd[SSS_BUFFER_SIZE];						// char stores single characters, int stores integers and float stores floating point numbers.
char sss_reply[SSS_BUFFER_SIZE];
unsigned int sss_buffer_size = SSS_BUFFER_SIZE;		// unsigned int represents an integer value without a sign
int SSS[8];

char iMTQ_Com[64];
char imtq_buffer[256];
unsigned int imtq_buffer_size=256;

int ch;

int running = 1;

static datetime_t now;								// A static variable preserves its previous value in its previous scope and is not initialized again in the new scope. 
static volatile unsigned int loopCounter;			// The volatile keyword in C is used to indicate to the compiler that a variable’s value may change unexpectedly. This is often the case when a variable is being accessed by multiple threads or when it represents hardware that is external to the computer.

// Storage for ADCS related parameters
static unsigned char Nmax = 13;						// Unsigned char is a character datatype where the variable consumes all the 8 bits of the memory and there is no sign bit (which is there in signed char). So it means that the range of unsigned char data type ranges from 0 to 255.
static double magnetometer_bias[3] = {				// double is a data type that is used to store high-precision floating-point data or numbers (up to 15 to 17 digits)
	0, 0, 0
}; // MISTmagm matlab model var_sim.m line 80
static bool biasest_enabled = true;
static double magnetometer_skew[3][3] = { { 1, 0, 0 },
					  { 0, 1, 0 },
					  { 0, 0, 1 } };
static double axismap_magneto[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
static double axismap_torquer[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
static double inertiatensor[3][3] = {
	{ 0.0370, 0, 0 },
	{ 0, 0.0510, 0 },
	{ 0, 0, 0.0205 }
}; //MIST inertial tensor from matlab
static double InverseInertia[3][3];

/*
static double inertiatensor_real[3][3] = {
	{ 0.0510, 0, 0 },
	{ 0, 0.0370, 0 },
	{ 0, 0, 0.0205 }
}; //MIST inertial tensor from matlab*/
static double sat_dip = 0.2;
static double cut_off = 0.01; // was 0
static double time_factor = 1;
static double threshold_angularrate = 0.1; // was 0.01
static double lambdabias = 1E-4; // was 0.1, 0.002 in demo
static double P0ValUARKF = 0.01; // was 0.1
static double QUARKF = 0.01;
static unsigned char KFSelection = 1;
static double UAKFProcessAdaptation = 0; //suggested value set to 0

// Noise Variables
// General Considerations:
	// 1) The code is very monolitic, consider transfering this piece to a function for better readability & memory considerations
	// 2) The variables for mean, variance and samplefreq should be global
	// 3) NOISE parameters could be given as text file input
int noiseCounter = 0;

//static double mean_mag;
//static double mean_sol;
static double variance_mag = 1E-9;
static double variance_sol = 1E-2;
static int NOISE_BOOL = 0; // 1 = noise activated, 0 = noise disactivated
//double whiteNoiseGauss[N];
//double sampledNoise[N];
double noise = 0.0; // default value, gets updated if NOISE_BOOL == 1
double temp;
double S,Z,U1,U2,u,v,fac;
int phase;
int seed;
int ii,jj;
double mean = 0.0; // DON'T CHANGE
double variance = 1.0; // DON'T CHANGE
int dummy;

//Input Variables
double mag_dipole_scaled[3]; //use to be long
double photodiodes[6];
static double mag_measurements[3];
double mag_in[3];
double sol_in[3];

#define Size 90000

//Variables for Albedo model
double R_earth = 6380;
double lat_span[8] = {-90,-73,-62,-28,21,65,79,90};
double albedo_coeff_span[8] = {0.7,0.65,0.47,0.22,0.22,0.42,0.57,0.61};
/*also refer to pd_maxalbedothreshold and pd_usealbedo*/
int albedoringsnumb[6] = {4,9,14,19,24,29}; //Vector with number of segments per ring. (6 rings, cf. "arings")
int unsigned segm = 0;
int unsigned arings = 6; //Number of rings
int sFlux = 1367; // Solar flux (W/m^2)
double astrounit = 1.495978707E8; //Mean sun distance in km
double sun_dist = 1.522180392575172E8;//distance from earth to the sun [Km] at a certain point in the orbit !!! further implementation: make this JD dependent.
double total_photodiodes[6];
double Albedo_PD_I[6]={0,0,0,0,0,0}; //sum of all albedo photodiodes currents.
double largest_element[Size]; // Max element of the albedo_threshold 

//Variables for Disturbance Torques computation (Explanations for these variables can be found in DisturbanceTorques.c)
double dnface[3][11] ={
		{1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1},
		{0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, -1, 1, -1, 1, -1, 0}
};
double darea[11]={0.02, 0.03, 0.01, 0.03, 0.03, 0.01, 0.03, 0.03, 0.03, 0.03, 0.01};
double drface[3][11]={
		{0.05, 0, 0, -0.05, 0, 0, 0.2, 0.2, -0.2, -0.2, 0.05},
		{0, 0.05, 0, 0, -0.05, 0, 0, 0, 0, 0, 0},
		{0.05, 0, 0.15, 0, 0, -0.15, 0.15, 0.15, 0.15, 0.15, -0.1}
};
double dskewomearth[3][3]={
		{0, -7.291E-5, 0},
		{7.291E-5,0,0},
		{0,0,0}
};
double muRc[3]= {0,0.05,0};
double cM[3] = {0,0,0};
double cD=2.1;
double sigS[3][11] = {
		{0.7500,0.7500,0.1500,0.7500,0.7500,0.1500,0.7500,0.1500,0.7500,0.1500,0.1500},
		{0.1700,0.1700,0.6900,0.1700,0.1700,0.6900,0.1700,0.6900,0.1700,0.6900,0.6900},
		{0.0800,0.0800,0.1600,0.0800,0.0800,0.1600,0.0800,0.1600,0.0800,0.1600,0.1600}
};
double sigIR[3][11] = {
		{1,1,1,1,1,1,1,1,1,1,1},
		{0,0,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0,0,0}
};
bool activation_flag = true; // True to include disturbance torques, this might require a reduction in integration time step.
bool rnd_dipole = false;
bool falsifier = false;
// Density variable for the exponential model
double Rho; 
// disturbance torques for the four runge/kutta timesteps
double tDT1[3]; // Total Disturbance torque
double TD[3];	// Aerodynamic Disturbance Torque
double TGG[3];	// Gravity Gradient Disturbance Torque
double TRD[3];  // Residual Dipole Disturbance Torque
double TS[3];   // Radiation Pressure Disturbance Torque    

//Matrices for text file writing
double Eclipse_flag_output[Size];
double quat_est_mat[Size][4];
double qReal_Body_mat[Size][4];
double quat_real_mat[Size][4];
double Dip_mat[Size][3];
double Tcon_mat[Size][3];
double B_mat[Size][3];
double Omega_mat[Size][3];
double Omega_sym_mat[Size][3];
double uSun_mat[Size][3];
double Csol_ECI[Size][3];
double CB_ECI[Size][3];
double dcmE2O[3][3]; // from ECI to OCF
double dcmO2E[3][3]; // from ECI to OCF  shouldn't this be from OCF to ECI
double uSun_ISIS[Size][3];

//Housekeping
int house_time[Size][6];
double house_inmag[Size][3];
double house_insol[Size][3];
int house_inmode[Size];
double house_indt[Size];
double house_calcpos[Size][3];
double house_calcvel[Size][3];
double house_calcsol[Size][3];
double house_calcmag[Size][3];
int house_calceclipse[Size];
double house_outq[Size][4];
double house_outvel[Size][3];
double house_outmagbias[Size][3];
double house_errorcovsum[Size];
int date[Size][6];


static bool pd_usealbedo = false; //true when albedo effect is included allows for albedo correction with UTILS library
//(Do not remove, used by the ISIS library as well)
static unsigned int pd_bias = 0;
static unsigned int pd_maxalbedothreshold = 1516169; // 1376763; 
//static unsigned int pd_maxalbedothreshold = 505;
static unsigned char eclipseflag;

double dt = 1; //time between meassuremnts
double actuation_duration = 1; //0.1 on isis demo change to matlab step

//MY ADDITION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// MODE SELECTION/SWITCHING PARAMETERS: 
// Modes to switch: (numbering as adcs_types.h file)
adcs_mode_tag_t MODE_VEC[] = {ADCS_MODE_OFF, ADCS_MODE_ARAM, ADCS_MODE_DETERMINATION, ADCS_MODE_NADIRPOINTING, ADCS_MODE_ZENITPOINTING, ADCS_MODE_BDOT, ADCS_MODE_BDOT_WITH_RATE_ESTIMATE, ADCS_MODE_ENGINEERING};

unsigned int flag = 0;						// integer to select current adcs mode


// ADCS_param TEXT FILE RELATED PARAMETERS. The structrure of this txt follow must be the following:
// Gain related parameters
#define adcs_param_count 13				// number of scalar adcs parameters to be specified.
double **adcs_param = NULL;	// initialization of matrix to store the different adcs related scalar parameters. A detailed list can be found on the doocument called ICD ADCS Ng library on Mist Cloud/Satellite/Attitude/Isis correspondence	
// adcs_param[i][0] = proportional gain (i stands for the different updated versions of this parameters along the simulation)
// adcs_param[i][1] = differential gain
// adcs_param[i][2] = bdot gain
// adcs_param[i][3] = initial error covariance value for the Kalman filter
// adcs_param[i][4] = noise figure for the normalized magnetic measurement
// adcs_param[i][5] = noise figure for the time differential of the normalized magnetic measurements
// adcs_param[i][6] = noise figure for the normalized solar measurements
// adcs_param[i][7] = noise figure for the time differential of the normalized solar measurements
// adcs_param[i][8] = quaternion process noise
// adcs_param[i][9] = angular rate process noise
static double gain_angularmomentum = 0.1;	// angular momentum gain
double *pointingvector = NULL;				// initialization of axis in the body frame to point at the inertial reference
double *switch_param = NULL;				// switching times of adcs parameters throughout the simulation 
unsigned int param_sets = 1;				// initialization of number of adcs parameter sets that will be used
unsigned int flag_adcs_upd = 0; 			// to select set of adcs parameters to update 

// Switch mode related parameters
int *adcs_modes = NULL;				// adcs modedouble *switch_modes = NULL;		// switching times of adcs modes throughout the simulation s used throughout the simulation
double *switch_modes = NULL;		// switching times of adcs modes throughout the simulation 
unsigned int modes_number = 1;		// initialization of number of adcs modes used during simulation

// TLE reset related parameters
static TLE_t **tlesets; 			// declare a pointer to a pointer to TLE_t
TLE_t tle;							// declare a static variable named tle of type TLE_t. The type TLE_t is defined in adcs_ephemeris.h
double *switch_TLE; 				// switching times of the TLE's used throughout the simulation 
unsigned int TLE_number = 1;		// initialization of number of adcs modes used during simulation
unsigned int flag_TLE_upd = 0; 		// to select set of TLEs to update (0-->initial set, 1-->second set, 2-->third set)
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

ADCS_HK_t adcs_hk;

//Output Variables
static double com_dip_out[3] = { 0 };
static double com_dip_imtq[3];
static double q_out[4] = { 0 };
static double ohm_out[3] = { 0 };
static double bias_est[3] = { 0 };

//Extra variables
char delim[] = ",";
unsigned int L = 1;
//Partial actuation auxiliary variables
unsigned int torque_actuation_counter;
unsigned char torque_actuation_flag;
unsigned int torque_actuation_desired_number;

void (*loopStepCallback)(adcs_task_tag_t) = NULL;

// Read ADCS_param text file (MY ADDITION) *************************************************
static int Retrieve_ADCS_Param(const char *filename, TLE_t ***tlesets, unsigned int *num_sets, double **switch_TLE, unsigned int *modes_number, int **adcs_modes, double **switch_modes, unsigned int *param_sets, double ***adcs_param, double **switch_param) {	// filename should be "..\\data\\ADCS_param.txt" 
	printf("ADCS parameters Loading started\n");
	printf("----------------------------------\n");
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Failed to open file ADCS_param.txt!\n");
        return 1;
    }

    // Read the number of TLE sets from the first line
	if (fscanf(file, "%d", num_sets) != 1) {
        printf("Failed to read the number of TLE sets.\n");
        fclose(file);
        return 1;
    }

    // Allocate memory for the array of TLE sets
    *tlesets = malloc(*num_sets * sizeof(TLE_t *));
	*switch_TLE = malloc(*num_sets * sizeof(double));
    if (*tlesets == NULL || *switch_TLE == NULL) {
        printf("Memory allocation of TLE sets or TLE switching times array failed\n");
        fclose(file);
        return 1;
    }

    // Read each line of TLE data and parse it into the corresponding TLE structure
    for (unsigned int i = 0; i < *num_sets; i++) {
        TLE_t *tleset = malloc(sizeof(TLE_t));
        if (tleset == NULL) {
            printf("Memory allocation of TLE set %d failed\n", i + 1);
            // Free memory allocated so far
            for (unsigned int j = 0; j < i; j++) {
                free((*tlesets)[j]);
            }
            free(*tlesets);
            fclose(file);
            return 1;
        }
        if (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
               &((*switch_TLE)[i]),&tleset->inclination, &tleset->RAAN, &tleset->e,
               &tleset->ohm, &tleset->m, &tleset->n, &tleset->BSTAR,
               &tleset->EPOCH_year, &tleset->EPOCH_day,
               &tleset->EPOCH_dayfraction) != 11) {
            printf("Failed to read TLE set %d.\n", i + 1);
            // Free memory allocated so far
            for (unsigned int j = 0; j <= i; j++) {
                free((*tlesets)[j]);
            }
            free(*tlesets);
            fclose(file);
            return 1;
        }
        (*tlesets)[i] = tleset;
    }

	// Read number of ADCS modes that will be used
    if (fscanf(file, "%d", modes_number) != 1) {
        printf("Failed to read the number of modes.\n");
        free(*tlesets);
		fclose(file);
        return 1;
    }

	// Allocate memory for the mode switching vector
	*adcs_modes = malloc(*modes_number * sizeof(int));
	*switch_modes = malloc(*modes_number * sizeof(double));
    if (*adcs_modes == NULL || *switch_modes == NULL) {
        printf("Memory allocation of modes or switching times failed.\n");
		free(*tlesets);
		free(*adcs_modes);
		free(*switch_modes);
        fclose(file);
        return 1;
    }

	// Read the modes and the switching modes times
    for (unsigned int i = 0; i < *modes_number; ++i) {
        if (fscanf(file, "%lf,%d", &((*switch_modes)[i]), &((*adcs_modes)[i])) != 2) {
            printf("Failed to read modes and switching times %d.\n", i + 1);
            free(*tlesets);
			free(*adcs_modes);
			free(*switch_modes);
        	fclose(file);
            return 1;
        }
    }

	// Read number of sets of ADCS parameters that will be used
    if (fscanf(file, "%d", param_sets) != 1) {
        printf("Failed to read the number ADCS parameters sets.\n");
        free(*adcs_modes);
        free(*switch_modes);
        free(*tlesets);
        fclose(file);
        return 1;
    }

	// Allocate memory for the adcs param matrix and the adcs param switching times
	*switch_param = malloc(*param_sets * sizeof(double));
    *adcs_param = malloc(*param_sets * sizeof(double *));
	if (*adcs_param == NULL || *switch_param == NULL) {
		printf("Memory allocation of adcs_param rows failed.\n");
		// Free previously allocated memory if necessary
		free(*adcs_modes);
		free(*switch_modes);
		free(*tlesets);
		free(*switch_param);
		free(*adcs_param);
		fclose(file);
		return 1;
	}

	// Allocate memory for each row of Adcs Parameters
	for (unsigned int i = 0; i < *param_sets; i++) {
		(*adcs_param)[i] = malloc(adcs_param_count * sizeof(double));
		if ((*adcs_param)[i] == NULL) {
			printf("Memory allocation of adcs_param row %d failed.\n", i);
			// Free previously allocated memory if necessary
			for (unsigned int j = 0; j < i; j++) {
				free((*adcs_param)[j]);
			}
			free(*adcs_modes);
			free(*switch_modes);
			free(*tlesets);
			free(*switch_param);
			free(*adcs_param);
			fclose(file);
			return 1;
		}
	}

	// Read the adcs parameters and the switching adcs parameters times
	for (unsigned int i = 0; i < *param_sets; ++i) {
		if (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &((*switch_param)[i]),
                &((*adcs_param)[i][0]), &((*adcs_param)[i][1]), &((*adcs_param)[i][2]), &((*adcs_param)[i][3]),
				 &((*adcs_param)[i][4]), &((*adcs_param)[i][5]), &((*adcs_param)[i][6]), &((*adcs_param)[i][7]),
				  &((*adcs_param)[i][8]), &((*adcs_param)[i][9]), &((*adcs_param)[i][10]), &((*adcs_param)[i][11]),
				  &((*adcs_param)[i][12])) != adcs_param_count + 1) {
			printf("Failed to read ADCS parameters and switching time for set %d.\n", i + 1);
			// Free previously allocated memory
            for (unsigned int j = 0; j < i; j++) {
                free((*adcs_param)[j]);
            }
            free(*adcs_param);
            free(*switch_param);
            free(*adcs_modes);
            free(*switch_modes);
            free(*tlesets);
            fclose(file);
            return 1;
		}
	}

    // Close the file
    fclose(file);
	printf("----------------------------------\n");
	printf("ADCS parameters Loaded Successfully \n");
	return 0; // Return success code
}
//*************************************************

static datetime_t MOCK_CurrentTime = {
	.year = 2017,
	.month = 6,
	.day = 21,
	.hours = 00,
	.minutes = 00,
	.seconds = 00.0
}; //MIST MATLAB date = launch date

// static datetime_t MOCK_CurrentTime2 = {
// 	.year = 2017,
// 	.month = 6,
// 	.day = 21,
// 	.hours = 2,
// 	.minutes = 13,
// 	.seconds = 20.0
// }; //MIST MATLAB date = launch date


static void MOCK_GetCurrentTime(datetime_t *now)
{
	memcpy(now, &MOCK_CurrentTime, sizeof(datetime_t));

	MOCK_CurrentTime.seconds =
		MOCK_CurrentTime.seconds + 1; //MATLAB runs with 0.5 second step

	if (MOCK_CurrentTime.seconds > 59) {
		MOCK_CurrentTime.minutes++;
		MOCK_CurrentTime.seconds = 0;
	}
	if (MOCK_CurrentTime.minutes > 59) {
		MOCK_CurrentTime.hours++;
		MOCK_CurrentTime.minutes = 0;
	}
	if (MOCK_CurrentTime.hours > 23) {
		MOCK_CurrentTime.day++;
		MOCK_CurrentTime.hours = 0;
	}
	if (MOCK_CurrentTime.day > 30) {
		MOCK_CurrentTime.month++;
		MOCK_CurrentTime.day = 1;
	}
	if (MOCK_CurrentTime.month > 12) {
		MOCK_CurrentTime.year++;
		MOCK_CurrentTime.month = 1;
	}
}

static void ADCSLoop_ProcessPDInput(double *pd_input, bool use_albedo,
					double *solvec_out)
{
	unsigned int adcs_pds[6] = { 0 };

	// X photodiodes
	adcs_pds[0] = pd_input[0];
	adcs_pds[1] = pd_input[1];

	// Y photodiodes
	adcs_pds[2] = pd_input[2];
	adcs_pds[3] = pd_input[3];

	// Z photodiodes
	adcs_pds[4] = pd_input[4];
	adcs_pds[5] = pd_input[5];

	// Convert photodiode inputs to a vector, optionally taking albedo into account
	if (use_albedo) {
		UTILS_sixPD2vectorWithAlbedo(adcs_pds, solvec_out);
	} else {
		UTILS_sixPD2vector(adcs_pds, solvec_out);
	}

	//**
	// * Process photodiode inputs and turns them into a solar vector
	// *
	// * @param[in] pd_input Raw photodiode inputs in BRF in the following order: X+, X-, Y+, Y-, Z+, Z-
	// * @param[in] pdselection Specification of which photodiode inputs to use and which one to set to zero
	// * @param[in] use_albedo If TRUE: take albedo into account when constructing solar vector
	// * @param[out] solvec_out Solar vector
	// */
}
/*
void transform_ECI_to_OCF(double rOrb[3], double vOrb[3], double uECI[3], double uOCF[3]) {
    // Compute the unit vectors of the OCF
    double x[3], y[3], z[3];
	LINALG_vector_normalise(3, rOrb, z);  // Unit(z) = rOrb / ||rOrb||
    LINALG_vector_cross_product(3, rOrb, 3, vOrb, y);  // y = rOrb x vOrb
    LINALG_vector_normalise(3, y, y);  // Unit(y) = y / ||y||
    LINALG_vector_cross_product(3, y, 3, z, x);  // x = y x z

    // Transformation matrix
    double m[3][3] = {{x[0], x[1], x[2]},
                      {y[0], y[1], y[2]},
                      {z[0], z[1], z[2]}};

    // Apply transformation
    LINALG_matrix_matXvec(3, 3, m, 3, uECI, uOCF);
}
*/
void ADCSLoop_Initialize(unsigned int flag_TLE_upd_local, unsigned int flag_adcs_upd_local)  // MY ADDITION: FUNCTION ARGUMENTS 
{
	// Check TLEs
	// MY ADDITION: select set of TLEs and ADCS parameters ********************
	// printf("\n\nflag_TLE_upd value: %d\n\n", flag_TLE_upd_local);
	
	double pointingvector[3] = {
    adcs_param[flag_adcs_upd_local][10], 
    adcs_param[flag_adcs_upd_local][11], 
    adcs_param[flag_adcs_upd_local][12]
	}; //the axis in the body frame to point at the inertial frame

	tle = *tlesets[flag_TLE_upd_local];

	// printf("\n*****************************************");
	// printf("\nSet of TLEs used: %d", flag_TLE_upd_local);
	// printf("\nSet of ADCS parameters used: %d", flag_adcs_upd_local);
	printf("\n***************************************** \n");
	printf("ADCSLoop_Initialize: Applying TLE reset for set %d\n", flag_TLE_upd_local);
    printf("New TLE Data: inclination = %f, RAAN = %f, EPOCH_year = %f, EPOCH_day = %f\n, EPOCH_dayfraction = %f\n" ,
           tle.inclination, tle.RAAN, tle.EPOCH_year, tle.EPOCH_day, tle.EPOCH_dayfraction);
	// printf("\n***************************************** \n");
	// //*******************************************************

	// Get current time
	//MOCK_GetCurrentTime(&now);
	// Initialize ADCS library
	ADCS_MAIN_initialise(&now, inertiatensor, &tle, Nmax, adcs_param[flag_adcs_upd_local][0],
				 adcs_param[flag_adcs_upd_local][1], adcs_param[flag_adcs_upd_local][2], gain_angularmomentum,
				 pointingvector, sat_dip, cut_off, time_factor,
				 adcs_param[flag_adcs_upd_local][3], adcs_param[flag_adcs_upd_local][4], adcs_param[flag_adcs_upd_local][5], adcs_param[flag_adcs_upd_local][6],
				 adcs_param[flag_adcs_upd_local][7], adcs_param[flag_adcs_upd_local][8], adcs_param[flag_adcs_upd_local][9],
				 threshold_angularrate, lambdabias, P0ValUARKF,
				 QUARKF, loopStepCallback, KFSelection,
				 UAKFProcessAdaptation);

	// Initialize PD values
	UTILS_set_photodiode_bias(pd_bias);
	UTILS_set_max_albedo_threshold(pd_maxalbedothreshold);

	// Enable / disable bias estimator
	// No need to check for false because call to ADCS_MAIN_initialise disables the bias estimator by default
	if (biasest_enabled == true) {
		DETERMINATION_bias_set_correction(1);
	}


}

BOOL WINAPI CtrlHandler(DWORD fdwCtrlType)
{
	switch (fdwCtrlType)
	{
		//Handle Ctrl-C Signal
		case CTRL_C_EVENT:
			printf("Ctrl-C event\n\n");
			running = 0;
			return TRUE;

		default:
			printf("Control something");
			return FALSE;
	}
}

//Main loop Start --------------------------------------------------------------------------------------------------------
int main(int argn, const char** argv)
{
	printf("ADCSNG test program\n\n");

if (SetConsoleCtrlHandler(CtrlHandler, TRUE))
	{
		printf("\n The control handler is installed\n");
		printf("\n -- You can end the simulation by pressing Ctrl-C");
		printf("\n	and this will save the data...\n");
		printf("\n(...waiting in a loop for events...)\n\n");

	}
	else
	{
		printf("\nERROR: Could not set control handler");
		return 1;
	}	
	
	time_t strtt = time(NULL); /* Start timer */
	
	Retrieve_ADCS_Param("..\\data\\ADCS_param.txt", &tlesets, &TLE_number, &switch_TLE , &modes_number, &adcs_modes, &switch_modes, &param_sets, &adcs_param, &switch_param);
	ADCSLoop_Initialize(flag_TLE_upd, flag_adcs_upd);  // FIRST INITIALIZATION // MY ADDITION: function argument input
	
	// The first argument corresponds to the port of the iMTQ simulator
	imtq = selector_init_port(0, argn, argv);
	// and the second to the port of the SSS
	sss = selector_init_port(1, argn, argv);

	if (selector_check(imtq) || selector_check(sss)) {
			// If any of the ports is open, we wait for the arduino to reset
			Sleep(2000);
	}

	HOUSEKEEP_display_frame();
	HOUSEKEEP_display_error_summary();

	//	HOUSEKEEP_initialise(&now);

	printf("Ends Initialization \n"
		   "---------------------\n");

	//Get Data from text file Param
	printf("Loading Text Files... \n\n");

	// unsigned int col = 3;
	// unsigned int colpS = 3; //if we use photodiodes this is set to 6

	unsigned int Parlength = 11;
	double Param[Parlength];
	// load Parameters file---------------------------------------------------------------------
	FILE *fp0 = fopen("..\\data\\Param.txt", "r");
	if (fp0 == NULL) {
		printf("Failed to open file Param.txt!\n");
		return 1;
	}
	for (unsigned int i = 0; i < Parlength; i++) {
		fscanf(fp0, "%lf", &Param[i]);
	}
	fclose(fp0);
	//--------------------------------------------------------------
	printf("Files Loaded Successfully \n\n");
	printf("ADCS_MAIN_adcs_service()\n");

	//Initial Sate
	double X[7] = { Param[2], Param[3], Param[4], Param[5],
			Param[6], Param[7], Param[8] };
	//double X[7]= {0,0,0,1,0,0,0}; /* Activate this line to manually input initial state vector */

	double stepsize = Param[0];
	unsigned int mod = Param[1];
	double endTime = Param[9];
	//double TotalData = (Time[row - 1] / stepsize) + 1;
	
	double BcurrentOCF[3], uSunCurrentOCF[3], rOrbCurrent[3],
		vOrbCurrent[3], B_field_eci[3], eci_solar_direction[3];
	
	double Eclipse_flag[1];
	double jD;
	double runningTime[] = { 0 };
	//double angleVec[3];
	double quat[4];
	double C_O2B[3][3];
	double qO2B[3];
	double uSunBody[3];
	double T_real[3];
	double wOrb[3];
	// double DxDt[7];

	unsigned int k = 0;
	double q[4];
	double qReal_ECI[4];
	double C_E2B[3][3];
	double dcmquat_Real[3][3];

	/*--------Albedo variables---------*/
	/*Variables needed for the calculations (dummies)*/
	double dum1[3],dum2[3],dum3[3],C_dum[3][3],C_dum1[3][3];
	double pr1 = 0;
	double pr2 = 0;
	double pr3 = 0;
	double dum4[3],PhDUM[6];
	/*------------------------------------*/
	/*-Find albedo points parameters, it's the view cone area's fragmentation.-*/
	for (unsigned int i=0;i<6;i++)
	{
		segm = segm + albedoringsnumb[i];
	}
	/*-For each ring we need the included area factor to calculate where the ring's middle line is w.r.t the rOrb and
	 *  the polar angle to identify it on that ring-*/
	double summa = 0;
	double A_ins_prov[arings];
	for (unsigned int i=0;i<arings;i++)
	{
		A_ins_prov[i] = ((double)albedoringsnumb[i]/2+summa)/segm; //area factor
		summa = summa + albedoringsnumb[i];
	}
	double A_inside[segm];
	int ct;
	double alpha[segm];
	//alpha and gamma are used to define the albedo vectors later
	for(unsigned int j=0;j<arings;j++)
	{
		bool pass = true;
		for(int i=0;i<albedoringsnumb[j];i++)
		{
			if (j==0)
			{
				ct = 0;
				pass = false;
			}
			else if(pass)
			{
				ct = ct + albedoringsnumb[j-1];
				pass = false;
			}
			A_inside[i+ct] = A_ins_prov[j];
			alpha[i+ct] = 2*M_PI*((double)i/albedoringsnumb[j]); //radians
		}
	}
	double cosgamma[segm],singamma[segm];
	double Albedo_PD_I_mat[segm][6]; //albedo photodiodes currents for all albedo points (direction and intensity)

	/*-----------------------*/

	double E[segm]; // Intensity vector (contains intensities of all 99 points)
	double b_r[3]; //albedo point direction vectors matrix
	/*(Non indexed values in the cycle are adjourned at every cycle step, noneed to allocate memory to a matrix of these)*/
	double b_r_unit[3];
	double a_r_unit[3]; //albedo points normal vector in relative frame
	double a_N_unit[3]; //albedo points normal vector in ECI frame (used to calculate latitude)
	double norm_i_b_r;
	double a_o_unit[3]; //albedo points normal vector in OCF frame
	double a_b_unit[3]; //albedo points normal vector in body frame
	double b_o_unit[3]; //albedo point direction vector in OCF frame
	double b_b_unit[3][segm]; //albedo points directions vectors matrix in body frame(99 vectors contained in columns,expressed in relative frame)
	double b_b[3][segm]; //albedo point direction vector matrix in body frame
	double b_o[3]; //albedo point direction vector matrix in OCF frame
	double C_OR[3][3] = {
			{0,1,0},
			{0,0,1},
			{1,0,0}

	}; //DCM R to OCF
	double latitude;
	double albedo_coeff;
	double norm_b_b;
	double Rel_sol_en = sFlux*(astrounit*astrounit)/(sun_dist*sun_dist);
	/*---------------------------------*/
	/*----------RK4 variables----------*/

	double T_real_int[3];
	/*--------------------------------*/
	// Calculate desired actuations number in one ADCS cycle
	torque_actuation_desired_number = mod * time_factor;

	/*
	// DON'T USE
	double t_temp1 = 0.0;
	//, t_temp2 = 0.0;  //Dummy variables to calculate step size	
	struct timeval st, et;			//Used for time measurments
	mingw_gettimeofday(&st, NULL);*/
	time_t strt_while = time(NULL);
	
	struct timeval st_dyn, et_dyn;
	int dyn_time = 0;

	//Simulation MAIN LOOP START  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double ohm[3]; // first three componentes of state vector 
	double quaternion4[4]; // last four components of state vector 
	
	// FILE *position_file = fopen("..\\data\\rOrbInloop.txt", "w");
	// FILE *magfile = fopen("..\\data\\magECI_igrf_ISIS.txt", "w");
	// FILE *sunfile = fopen("..\\data\\sunECI_ISIS.txt", "w");
	FILE *measurementM = fopen("..\\data\\mag_measurement_body.txt", "w");
	FILE *measurementS = fopen("..\\data\\uSun_measurement_body.txt", "w");
	while (running) {
		//for (unsigned int i = 0; i < 3; i++) 
			//angleVec[i] = X[i];
		//
		
		// attaining SGP4 parameters for such as position, velocity, magnetic field etc at every whole timestep
		if(loopCounter % mod == 0) // The modulus operator in C (%) calculates the remainder of the division of the left operand by the right operand
		{
		
			//double total_seconds_passed2 = now.seconds + now.minutes*60 + now.hours*3600; 
			//printf("Current value of k: %d \n", k);

			double eci_position_km[3], eci_velocity_kms[3];
			for (unsigned int i = 0; i < 6; i++) {
				date[k][0] = now.day;
				date[k][1] = now.month;
				date[k][2] = now.year;
				date[k][3] = now.hours;
				date[k][4] = now.minutes;
				date[k][5] = now.seconds;
			}

			// TLE SET SELECTION
				if (runningTime[0]>switch_TLE[flag_TLE_upd+1] && flag_TLE_upd<(TLE_number-1)){
					printf("\nUPDATING TLEs SET...<<<<<<<<<<<<<\n"); 
					flag_TLE_upd++;
					
					tle = *tlesets[flag_TLE_upd];
					ADCSLoop_Initialize(flag_TLE_upd, flag_adcs_upd);
					EPHEMERIS_sgp4_reset (&tle);
					printf("Current time: %04d-%02d-%02d %02d:%02d:%02d", now.year, now.month, now.day, now.hours, now.minutes, now.seconds);
					printf("\n*****************************************\n\n");
				}
			
			// static datetime_t now_copy;
			// memcpy(&now_copy, &now, sizeof(datetime_t));
			// if ((now.hours == 2) && (now.minutes < 23)) {
			// printf("Time: %04d-%02d-%02d %02d:%02d:%02d\n", 
       		// now.year, now.month, now.day, now.hours, now.minutes, now.seconds);
			// printf("Position ECI (km): [%f, %f, %f]\n", eci_position_km[0], eci_position_km[1], eci_position_km[2]);
			// }
			EPHEMERIS_sgp4(&now, eci_position_km, eci_velocity_kms); // call for SGP4 orbit propagation model
			// if ((now.hours == 2) && (now.minutes < 23)) {
			// printf("Time: %04d-%02d-%02d %02d:%02d:%02d\n", 
       		// now.year, now.month, now.day, now.hours, now.minutes, now.seconds);
			// printf("Position ECI (km): [%f, %f, %f]\n", eci_position_km[0], eci_position_km[1], eci_position_km[2]);
			// }
			//now = now_copy;
			// fprintf(position_file, "%f %f %f\n", eci_position_km[0], eci_position_km[1], eci_position_km[2]);
			// Get the solar direction in ECI frame
			// for (unsigned int i = 0; i < 3; i++) {
			// 	house_rECI2[k][i] = eci_position_km[i];
			// }
			//Calculate the Julian Date
			double julian_date = TIME_get_julian_date(&now);
			jD = julian_date;
			//double eci_solar_direction[3];
			
			//for (unsigned int j = 0; j < col; j++) {
				//eci_solar_direction[j] = uSun[k][j];
			//}
			
			EPHEMERIS_get_solar_direction(julian_date, J2000, now.year, eci_solar_direction);
			//fprintf(sunfile, "%f %f %f\n", eci_solar_direction[0], eci_solar_direction[1], eci_solar_direction[2]);
			
			// Check if the satellite is in eclipse
			// create a copy of position vector since check_eclipse alters the position vector
			double eci_position_km_copy[3];
			memcpy(eci_position_km_copy, eci_position_km, sizeof(eci_position_km));
			unsigned char is_eclipse = EPHEMERIS_check_eclipse(eci_position_km_copy, eci_solar_direction);
			is_eclipse ^= 1; // switch the value to be consistent with the disturbance.c torque 
		
			// Get the magnetic field in ECI frame in T	
			EPHEMERIS_igrf(&now, eci_position_km, B_field_eci);
			//fprintf(magfile, "%f %f %f\n", B_field_eci[0], B_field_eci[1], B_field_eci[2]);
			
			// going back to variables used in the rest of the code, may refine later and use less variables
			for (int i = 0; i < 3; i++) {
				BcurrentOCF[i] = B_field_eci[i];
				uSunCurrentOCF[i] = eci_solar_direction[i];
				rOrbCurrent[i] = eci_position_km[i];
				vOrbCurrent[i] = eci_velocity_kms[i];
			}
			Eclipse_flag[0] = (double)is_eclipse;
			//printf("eclipseFlag should be zero as with matlab : %f\n", Eclipse_flag[0]);

		}
		// setting initial conditions 
		if (*runningTime == 0 ) {
			q[3] = X[3]; 
			q[0] = X[4]; 
			q[1] = X[5]; 
			q[2] = X[6];
			//printf("initial position: [%f, %f, %f]\n", rOrbCurrent[0], rOrbCurrent[1], rOrbCurrent[2]);
			wOrbcalc(rOrbCurrent, vOrbCurrent, wOrb);  
			UTILS_q2dcm(q, C_O2B); // because initial condition is in OCF to body frame so you get dcm O2B
			ECI2OCF(rOrbCurrent, vOrbCurrent, dcmE2O);

			LINALG_matrix_matXvec(3, 3, C_O2B, 3, wOrb, ohm); // oscar bylund transformation if we still wanna use wOrb somehow
			for (int i = 0; i < 3; i++) {
				X[i] = X[i] + ohm[i]; // initial state angular velocities in bodyToInertial frame 
			}
			printf("initial angular velocities: [%f, %f, %f]\n", X[0], X[1], X[2]);
			// UTILS_dcm2q(dcmE2O_transpose, quaternion4); // quaternions representing the rotation bw OCF to ECI 
			LINALG_matrix_matXmat(3, 3, C_O2B, 3, 3, dcmE2O,
							C_E2B); // DCM of the rotation from ECI to body  = C_E2B

			UTILS_dcm2q(C_E2B, quaternion4); // quaternion representation of the rotation from ECI to body frame
			// also done in CASPER thesis `
			X[3] = quaternion4[3];
			X[4] = quaternion4[0];
			X[5] = quaternion4[1];
			X[6] = quaternion4[2]; 
		
			double normQuat = sqrt(pow(X[3], 2) + pow(X[4], 2) +
					   pow(X[5], 2) + pow(X[6], 2));
					   
			// unnecessary but good practice to norm quaternions 
			for (unsigned int i = 3; i < 7; i++) {
				X[i] = X[i] / normQuat;
			}
			printf("Initial orientation in ECI to body: [%f, %f, %f, %f]\n", quaternion4[3], quaternion4[0], quaternion4[1], quaternion4[2]);
		}
		
		// printf("print index k and time: %d, %f\n", k, runningTime[0]);
		Eclipse_flag_output[k] = *Eclipse_flag; // Output
		
		for (unsigned int i = 3; i < 7; i++) {
			quat[i - 3] = X[i]; // here X[i]
		}
		double quatturn[4] = { quat[1], quat[2], quat[3], quat[0] };
		//printf("quat in ECI to body: [%f, %f, %f, %f]\n", quatturn[3], quatturn[0], quatturn[1], quatturn[2]);
		UTILS_q2dcm(quatturn, C_E2B); // DCM of the rotation from ECI to body  = C_E2B
		ECI2OCF(rOrbCurrent, vOrbCurrent, dcmE2O);
		LINALG_matrix_transpose(3, 3, dcmE2O, dcmO2E);
		LINALG_matrix_matXmat(3, 3, C_E2B, 3, 3, dcmO2E, C_O2B);

		//		 printf("velocity orbit %f \n",vOrbCurrent[1]);
		LINALG_matrix_matXvec(3, 3, C_E2B, 3, B_field_eci,
					  mag_measurements);
		LINALG_matrix_matXvec(3, 3, C_E2B, 3, eci_solar_direction, uSunBody);
		
		if(loopCounter % mod == 0)
		{
			fprintf(measurementM, "%f %f %f\n", mag_measurements[0], mag_measurements[1], mag_measurements[2]);
			fprintf(measurementS, "%f %f %f\n", uSunBody[0], uSunBody[1], uSunBody[2]);
			
				}
		// Correct for magnetometer bias and skew and then transform to BRF
		UTILS_correct_mag_measurement(
			(const double *)magnetometer_bias,
			(const double(*)[3])magnetometer_skew, mag_measurements,
			ADCSLOOP_MAGCORRPROCEDURE, mag_measurements);
		LINALG_matrix_matXvec(3, 3, axismap_magneto, 3,
					  mag_measurements, mag_in);

		//Send "mag_in to iMTQ"

		//printf("\nmag_in \n%.9f %.9f %.9f\n", mag_in[0], mag_in[1], mag_in[2]);

		//-------------------------------------------- NOISE_BEGIN----------------

		if(NOISE_BOOL == 1){
			do{
					U1 = (double)rand()/RAND_MAX;
					U2 = (double)rand()/RAND_MAX;
					u = 2.0*U1-1.0;
					v = 2.0*U2-1.0;
					S = u*u+v*v;
				} while (S>=1);
			fac = sqrt(-2.0*log(S)/S);
			phase = rand()%2;
			if(phase)
				Z = v*fac;
			else
			{
				Z = u*fac;
			}
			phase = 1-phase;
			noise = Z*variance+mean;
			//printf("\n noise: %6.15f \n", noise);
		}

		//---------------------------------------------- NOISE_END---------------------

		//printf("\n noise: %6.15f \n", sampledNoise[noiseCounter]);
		//printf("\nmag: %d\n", noiseCounter);
		for(ii = 0; ii<3; ii++){
			//printf("\n original %d : %6.15f", ii, mag_in[ii]);
			mag_in[ii] = mag_in[ii] + noise*NOISE_BOOL*variance_mag;
			//printf("\n updated %d : %6.15f", ii, mag_in[ii]);
		}

		int mtm_x = mag_in[0] * 1E9;
		//printf("\nmtm_x %i\n", mtm_x);
		int mtm_y = mag_in[1] * 1E9;
		int mtm_z = mag_in[2] * 1E9;
  		//:mtm_c;x;y;z\n
  		snprintf(iMTQ_Com, sizeof(iMTQ_Com), "mtm_c;%i;%i;%i", mtm_x, mtm_y, mtm_z);
  		//printf("%s\n", iMTQ_Com);
  		WriteToUart(imtq, iMTQ_Com);
  		//imtq_buffer_size=256;
  		//WriteReceiveUart(imtq, "act_d", imtq_buffer, &imtq_buffer_size);
  		//printf("%s\n", imtq_buffer);
		//---------------------

		//Calculate Albedo contribution to photodiodes reading
		if (pd_usealbedo)
		{ /*-albedo points parameters used by the model are calculated once above.-*/

		/*--Position depenent parameters*/
		double r_dist = sqrt(pow(rOrbCurrent[0],2)+pow(rOrbCurrent[1],2)+pow(rOrbCurrent[2],2));
		double costeta = R_earth/r_dist;
		for(unsigned int i=0;i<segm;i++)
		{
			cosgamma[i] = 1-A_inside[i]*(1-costeta);
			singamma[i] = sqrt(1-pow(cosgamma[i],2));
		}
		/*-------------------------------------*/

		/*-Find direction vector and intenity from each albedo point-*/
		//Calcultions are made in a relative frame according to literature and the rotated to body

		for(unsigned int i=0;i<segm;i++)
		{
			a_r_unit[0] = cosgamma[i];
			a_r_unit[1] = singamma[i]*cos(alpha[i]);
			a_r_unit[2] = singamma[i]*sin(alpha[i]);

			b_r[0] = a_r_unit[0]*R_earth-r_dist;
			b_r[1] = a_r_unit[1]*R_earth;
			b_r[2] = a_r_unit[2]*R_earth;
			norm_i_b_r = sqrt(pow(b_r[0],2)+pow(b_r[1],2)+pow(b_r[2],2));

			b_r_unit[0] = b_r[0]/norm_i_b_r;
			b_r_unit[1] = b_r[1]/norm_i_b_r;
			b_r_unit[2] = b_r[2]/norm_i_b_r;

			LINALG_matrix_matXvec(3, 3, C_OR, 3, a_r_unit, a_o_unit);
			LINALG_matrix_matXvec(3, 3, C_O2B, 3, a_o_unit, a_b_unit);
			LINALG_matrix_matXvec(3, 3, C_OR, 3, b_r_unit, b_o_unit);
			LINALG_matrix_matXvec(3, 3, C_O2B, 3, b_o_unit, dum1);
			b_b_unit[0][i] = dum1[0];b_b_unit[1][i] = dum1[1];b_b_unit[2][i] = dum1[2];
			LINALG_matrix_matXvec(3, 3, C_OR, 3, b_r, b_o);
			LINALG_matrix_matXvec(3, 3, C_O2B, 3, b_o, dum2);
			b_b[0][i] = dum2[0];b_b_unit[1][i] = dum2[1];b_b_unit[2][i] = dum2[2];

			ECI2OCF(rOrbCurrent,vOrbCurrent, C_dum);
			LINALG_matrix_transpose(3, 3, C_dum,C_dum1);
			LINALG_matrix_matXvec(3, 3, C_dum1, 3, a_o_unit, a_N_unit);
			latitude = (180/M_PI)* asin(a_N_unit[2]); //updated at every cycle iteration, there's no latitude vector for every point.

			interp1(lat_span,8,albedo_coeff_span,&latitude,1,&albedo_coeff); //check interp vectors and function calling, i.e. pointers
			dum3[0]=-b_b_unit[0][i];dum3[1]=-b_b_unit[1][i];dum3[2]=-b_b_unit[2][i];
			pr1 =dot_product1(a_b_unit,dum3,3);
			pr2 =dot_product1(a_b_unit,uSunBody,3);
			if(pr2<0)
			{
				pr2 = 0;
			}
			pr3 = pow(R_earth,2);
			norm_b_b = sqrt(pow(b_b[0][i],2)+pow(b_b[1][i],2)+pow(b_b[2][i],2));
			E[i] = ((2*Rel_sol_en*pr3*albedo_coeff)/(segm*norm_b_b*norm_b_b))*(1-costeta)*pr1*pr2;
			/*This cycle outputs an intensity vector E[99] and two matrices containing the direction albedo vectors on the columns*/
		}
		/*---------------------------------------------*/

		/*-Find the PD currents due to every albedo vector -*/
		LINALG_matrix_zero(segm,6,Albedo_PD_I_mat);
		for (unsigned int i = 0; i < 6; i++)
		{
			Albedo_PD_I[i] = 0; // Resets for every point of the orbit
		}
		for (unsigned int i=0;i<segm;i++)
		{
			dum4[0]=b_b_unit[0][i];dum4[1]=b_b_unit[1][i];dum4[2]=b_b_unit[2][i];
			Albedo2PhD(dum4,E[i],Rel_sol_en,PhDUM,Imax);
			Albedo_PD_I_mat[i][0] = PhDUM[0];
			Albedo_PD_I_mat[i][1] = PhDUM[1];
			Albedo_PD_I_mat[i][2] = PhDUM[2];
			Albedo_PD_I_mat[i][3] = PhDUM[3];
			Albedo_PD_I_mat[i][4] = PhDUM[4];
			Albedo_PD_I_mat[i][5] = PhDUM[5];

			Albedo_PD_I[0] = Albedo_PD_I[0] + Albedo_PD_I_mat[i][0];
			Albedo_PD_I[1] = Albedo_PD_I[1] + Albedo_PD_I_mat[i][1];
			Albedo_PD_I[2] = Albedo_PD_I[2] + Albedo_PD_I_mat[i][2];
			Albedo_PD_I[3] = Albedo_PD_I[3] + Albedo_PD_I_mat[i][3];
			Albedo_PD_I[4] = Albedo_PD_I[4] + Albedo_PD_I_mat[i][4];
			Albedo_PD_I[5] = Albedo_PD_I[5] + Albedo_PD_I_mat[i][5];

		}
		////////// 2 Options: 
		// 1) Set albedo threshold equal to max photodiode reading largest_element[k];
		// 2) Set it to an average value calculated from Matlab
		largest_element[k] = Albedo_PD_I[0];
		for (unsigned int i = 1; i<6; i++)
		{
			if (largest_element[k] < Albedo_PD_I[i] )
			{
				largest_element[k] = Albedo_PD_I[i];
		 	}
		}
		// pd_maxalbedothreshold = 1376763.50886009; //largest_element[k];
		// UTILS_set_max_albedo_threshold(pd_maxalbedothreshold);

		}
		/*---------------------------------------------*/

		//Convert Sun vector to Photodiode Meassurement
		S2PD(uSunBody, photodiodes, Imax);
		//printf("\nsun 1: \n%5.6f %5.6f %5.6f",uSunBody[0],uSunBody[1],uSunBody[2]);

		//Find total photodiode reading
		for(unsigned int i=0;i<6;i++)
		{
			total_photodiodes[i] = (photodiodes[i] + Albedo_PD_I[i])* *Eclipse_flag;
			total_photodiodes[i] = total_photodiodes[i]*(1+noise*NOISE_BOOL*variance_sol);
		}

		//Calculate and send voltages expected on OBC ADC input to SSS
		PDCurrToADC(total_photodiodes, Imax, SSS);
		snprintf(sss_cmd, sizeof(sss_cmd), "setchans;%i;%i;%i;%i;%i;%i;%i;%i", SSS[0], SSS[1], SSS[2], SSS[3], SSS[4], SSS[5], SSS[6], SSS[7]);
		//printf("%s\n", sss_cmd);
		WriteReceiveUart(sss, sss_cmd, sss_reply, &sss_buffer_size);
		if (sss_reply[9] != '0' && selector_check(imtq)) {
			printf("setchans: error %c\n", sss_reply[9]);
		}
		
		//Photodiode conversion, input selection, and solar vector construction
		ADCSLoop_ProcessPDInput(total_photodiodes, pd_usealbedo, sol_in);
		for (unsigned int i = 0; i < 3; i++) {
					uSun_ISIS[k][i] = sol_in[i];
		}
		//printf("\nsun2: \n%5.6f %5.6f %5.6f", sol_in[0],sol_in[1],sol_in[2]);

		if (loopCounter % mod == 0) //Avoid 2 time call due to stepsize in solver
		{
			//Get Current time
		
			//MOCK_GetCurrentTime(&now);
			// now.secon
			//double total_seconds_passed2 = now.seconds + now.minutes*60 + now.hours*3600; 
			//printf("Current time: %f \n", total_seconds_passed2);

			if (selector_check(imtq)) {
				//Read commanded dipole from iMTQ
				int act_status = -1;
				
				// Timer for variable sleep time
				mingw_gettimeofday(&et_dyn, NULL); //End timer
				if (k != 0) {
						dyn_time = (int)(1000*((et_dyn.tv_sec - st_dyn.tv_sec) + (et_dyn.tv_usec - st_dyn.tv_usec)/1000000.0));
				}
				do {
					Sleep((int)((1000-dyn_time)/5));
					//Sleep(200);
					imtq_buffer_size=256;
					WriteReceiveUart(imtq, "act_d", imtq_buffer, &imtq_buffer_size);
					printf("RX: %s\n", imtq_buffer);
					sscanf(imtq_buffer, "act_d;%d;%le;%le;%le;%lf", &act_status,  &com_dip_imtq[0], &com_dip_imtq[1], &com_dip_imtq[2], &actuation_duration);
				} while (act_status != 0);
				mingw_gettimeofday(&st_dyn, NULL);
				printf("Received act_d: %d, %d, %d, %d\n",  (int)com_dip_imtq[0], (int)com_dip_imtq[1], (int)com_dip_imtq[2], (int)actuation_duration);
				WriteToUart(imtq, "clean_d");

				static double axismap_imtq[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

				// Scale received commanded dipoles from iMTQ Sim.
				for (unsigned int i = 0; i < 3; i++) {
					com_dip_imtq[i] = com_dip_imtq[i]*0.0001;
				}

				LINALG_matrix_matXvec(3, 3, axismap_imtq, 3, com_dip_imtq,
							  com_dip_out);				

				// Set actuation counter to zero since new torque has been calculated, i.e. new ADCS cycle started
				torque_actuation_counter = 0;

				// Storing variables in arrays then copied in .txt files, saved here to pick the value at integer time step
				for (unsigned int i = 0; i < 3; i++) {
					Dip_mat[k][i] = adcs_hk.output_comdip[i];
					Tcon_mat[k][i] = adcs_hk.output_torque[i];
					Csol_ECI[k][i] = adcs_hk.calc_soleci[i];
					CB_ECI[k][i] = adcs_hk.calc_mageci[i];
					Omega_mat[k][i] = adcs_hk.output_ohmest[i];
					B_mat[k][i] = adcs_hk.input_magbrf[i];
					uSun_mat[k][i] = adcs_hk.input_solbrf[i];
					Omega_sym_mat[k][i] = X[i]; // what state is that i.e surely from body to ECI angular velocities
				}
				
				//Convert qreal to ECI and Save estimated quaternion to compare with matlab (!! Scalar part is the last quat.)
				// since elsewise it is taken care of in the initiocal condition i.e. for intial time 
				if (*runningTime != 0) {
				q[0] = X[4];
				q[1] = X[5];
				q[2] = X[6];
				q[3] = X[3];
				}

				ECI2OCF(rOrbCurrent, vOrbCurrent, dcmE2O);

				UTILS_q2dcm(q, dcmquat_Real);
				LINALG_matrix_matXmat(3, 3, dcmquat_Real, 3, 3, dcmE2O,
							  C_E2B); // DCM of the rotation from ECI to body  
				UTILS_dcm2q(C_E2B, qReal_ECI); // quaternion representation of the rotation from ECI to body frame

				for (unsigned int i = 0; i < 4; i++) {
					qReal_Body_mat[k][i] = q[i];
					quat_real_mat[k][i] = qReal_ECI[i];
					quat_est_mat[k][i] = adcs_hk.output_qest[i];
				}
				k++;

			} else {

				// Call to ADCS library
				// printf("Eclipse without switching is sent: %d\n", eclipseflag);
				//printf("\n CURRENT TIME: %f \n", runningTime[0]); <<<<<<


				// MODE SELECTION (MY ADDITION) <<<< <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				if (runningTime[0]>switch_modes[flag+1] && flag<(modes_number-1)){
					printf("\nSWITCHING MODE...<<<<<<<<<<<<<<<<\n\n");  
					flag++;			
				}
		
			    int current_mode = adcs_modes[flag];  
	            adcs_mode_tag_t MODE = MODE_VEC[current_mode];  

				// // TLE SET SELECTION
				// if (runningTime[0]>switch_TLE[flag_TLE_upd+1] && flag_TLE_upd<(TLE_number-1)){
				// 	printf("\nUPDATING TLEs SET...<<<<<<<<<<<<<\n"); 
				// 	flag_TLE_upd++;

				// 	tle = *tlesets[flag_TLE_upd];
				// 	//EPHEMERIS_sgp4_reset(&tle);
				// 	ADCSLoop_Initialize(flag_TLE_upd, flag_adcs_upd);
				// 	printf("Current time: %04d-%02d-%02d %02d:%02d:%02d", now.year, now.month, now.day, now.hours, now.minutes, now.seconds);
				// 	printf("\n*****************************************\n\n");
				// }

				// ADCS PARAMETERS SET SELECTION 
				if (runningTime[0]>switch_param[flag_adcs_upd+1] && flag_adcs_upd<(param_sets-1)){
					printf("\nUPDATING ADCS PARAMETERS SET...<<<<<<<<<<<<<\n"); 
					flag_adcs_upd++;
					
					ADCSLoop_Initialize(flag_TLE_upd, flag_adcs_upd);
					printf("Current time: %04d-%02d-%02d %02d:%02d:%02d", now.year, now.month, now.day, now.hours, now.minutes, now.seconds);
					printf("\n*****************************************\n\n");
				}
				// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

				ADCS_MAIN_adcs_service(&now, mag_in, sol_in, MODE, dt,
							   com_dip_out, q_out, ohm_out,
							   bias_est, &eclipseflag,
							   actuation_duration);

				// Set actuation counter to zero since new torque has been calculated, i.e. new ADCS cycle started
				torque_actuation_counter = 0;

				// Storing variables in arrays then copied in .txt files, saved here to pick the value at integer time step
				for (unsigned int i = 0; i < 3; i++) {
					Dip_mat[k][i] = adcs_hk.output_comdip[i];
					Tcon_mat[k][i] = adcs_hk.output_torque[i];
					Csol_ECI[k][i] = adcs_hk.calc_soleci[i];
					CB_ECI[k][i] = adcs_hk.calc_mageci[i];
					Omega_mat[k][i] = adcs_hk.output_ohmest[i];
					B_mat[k][i] = adcs_hk.input_magbrf[i];
					uSun_mat[k][i] = adcs_hk.input_solbrf[i];
					Omega_sym_mat[k][i] = X[i];
				}
				
				//Convert qreal to ECI and Save estimated quaternion to compare with matlab (!! Scalar part is the last quat.)
				// just making sure real part is at the right location, therefore using q as bucket 
				q[0] = X[4];
				q[1] = X[5];
				q[2] = X[6];
				q[3] = X[3];
				UTILS_dcm2q(C_O2B, qO2B); // quaternions denotes the rotation to transform from OCF to body
				for (unsigned int i = 0; i < 4; i++) {
					quat_real_mat[k][i] = q[i];
					qReal_Body_mat[k][i] = qO2B[i];
					quat_est_mat[k][i] = adcs_hk.output_qest[i];
				}
				k++;
			}
			

			// Convert output for use by iMTQ
			LINALG_matrix_matXvec(3, 3, axismap_torquer, 3, com_dip_out,
						  com_dip_imtq);


			//------------------------------

			//HOUSE KEEPING DATA
			HOUSEKEEP_sync();
			HOUSEKEEP_retrieve(&adcs_hk);

			house_time[k-1][0] = adcs_hk.timestamp.day;
			house_time[k-1][1] = adcs_hk.timestamp.month;
			house_time[k-1][2] = adcs_hk.timestamp.year;
			house_time[k-1][3] = adcs_hk.timestamp.hours;
			house_time[k-1][4] = adcs_hk.timestamp.minutes;
			house_time[k-1][5] = adcs_hk.timestamp.seconds;

			house_outq[k-1][0] = adcs_hk.output_qest[0];
			house_outq[k-1][1] = adcs_hk.output_qest[1];
			house_outq[k-1][2] = adcs_hk.output_qest[2];
			house_outq[k-1][3] = adcs_hk.output_qest[3];

			house_inmode[k-1] = adcs_hk.input_mode;
			house_indt[k-1] = adcs_hk.input_dt;
			house_calceclipse[k-1] = adcs_hk.calc_eclipse;
			house_errorcovsum[k-1] = adcs_hk.error_covariance_sum;

			for (unsigned int i = 0; i < 3; i++) {
				house_inmag[k-1][i] = adcs_hk.input_magbrf[i]*1E9;
				house_insol[k-1][i] = adcs_hk.input_solbrf[i];
				house_calcpos[k-1][i] = adcs_hk.calc_poskm[i];
				house_calcvel[k-1][i] = adcs_hk.calc_velkms[i];
				house_calcsol[k-1][i] = adcs_hk.calc_soleci[i];
				house_calcmag[k-1][i] = adcs_hk.calc_mageci[i]*1E9;
				house_outvel[k-1][i] = adcs_hk.output_ohmest[i];
				house_outmagbias[k-1][i] = adcs_hk.output_magbiasest[i];
			}
			MOCK_GetCurrentTime(&now); // update time
		}

		//////////////////////////////////////////////////Euler Integrator///////////////////////////////////////////////////////
		///////////////////////////////////////////////////RK4-integrator////////////////////////////////////////////////////////
		/* Runge-Kutta 4 integrator: EOM fnc produces the slopes required by the algorithm, attitude dependant variables
		   are calculated even at the mid step*/

		//---------------------------------------------------------------------
		LINALG_vector_cross_product(3, com_dip_out, 3, mag_in, T_real);

		// Check if actuation is required at this sub-cycle number and set flag value
		if( torque_actuation_counter == 0 && torque_actuation_counter < torque_actuation_desired_number){
			torque_actuation_flag = 1;
			torque_actuation_counter ++;
		}
		else if(torque_actuation_counter < torque_actuation_desired_number){
			torque_actuation_counter ++;
		}
		else{
			torque_actuation_flag = 0;
			torque_actuation_counter ++;
		}

		assert(torque_actuation_counter <= mod);
		assert(torque_actuation_flag == 1 || torque_actuation_flag == 0);
		
		// since disturbance torque is constant during one second
		if (loopCounter % mod == 0) {
			Rho = ExpAtmModel(rOrbCurrent, R_earth);	// Calculating current density with the exponential model
			/*-Add disturbance torques-*/
			DistTorq(activation_flag,rnd_dipole,inertiatensor,muRc,cD,cM,sigS,sigIR,C_O2B,rOrbCurrent,vOrbCurrent,
				BcurrentOCF,uSunCurrentOCF,*Eclipse_flag,jD,Rho,dnface,darea,drface,dskewomearth,tDT1,TD,TGG,TRD,TS);
			// printf("%6.17f \n", Dens);
		}
		for (unsigned int i=0;i<3;i++)
		{
			T_real_int[i] = T_real[i] * torque_actuation_flag + tDT1[i];
		}
		if (rnd_dipole)
		{
			rnd_dipole=false; // falsifies the randomness of the dipole for the mid steps
			falsifier=true;
		}
		/*-End dist. torques-*/
		//____________NEW VARIABLES for ISIS propagators _____________//
		double ohm_new[3];
		double q_new[4];
		double qF_out[4][4];
		double h_internal[3]; // cubesat's total internal angular momentum
		MatrixInverse(inertiatensor,InverseInertia); // inverse of the inertia is calculated
		for (int i = 0; i < 3; i++) {
			h_internal[i] = 0.0;
			for (int j = 0; j < 3; j++) {
				h_internal[i] += inertiatensor[i][j] * ohm[j];
			}
		}

		ECI2OCF(rOrbCurrent, vOrbCurrent, dcmE2O);
		// Propagate angular velocity
		UTILS_propagate_angular_rate(ohm, T_real_int, h_internal, inertiatensor,InverseInertia, stepsize, ohm_new);

		// Propagate quaternion
		UTILS_propagate_quaternion(quaternion4, ohm, T_real_int, inertiatensor, InverseInertia, h_internal, stepsize, q_new, qF_out);
		
		// Update values for the next iteration
		for (int i = 0; i < 3; i++) {
			ohm[i] = ohm_new[i];
		}
		for (int i = 0; i < 4; i++) {
			quaternion4[i] = q_new[i];
		}
		
		X[0] = ohm[0];
		X[1] = ohm[1];
		X[2] = ohm[2];
		X[3] = quaternion4[3]; //  fourth is the real 
		X[4] = quaternion4[0];
		X[5] = quaternion4[1];
		X[6] = quaternion4[2];

		/////////////////////////////////////////////////////////////End of Integration////////////////////////////////////////////////////////////

		double normQuat = sqrt(pow(X[3], 2) + pow(X[4], 2) +
					   pow(X[5], 2) + pow(X[6], 2));
		for (unsigned int i = 3; i < 7; i++) {
			X[i] = X[i] / normQuat;
		}

		if (*runningTime > endTime) {
			printf("--> Break while loop\n");
			break;
		}

		// DON'T USE
		//Calculate the time step from the measured time
		/*if (loopCounter % mod == 0 && selector_check(imtq)) {
			mingw_gettimeofday(&et, NULL);
			t_temp1 =  (et.tv_sec - st.tv_sec) + (et.tv_usec - st.tv_usec)/1000000.0;
			stepsize = (t_temp1) / mod;
			mingw_gettimeofday(&st, NULL);

			//printf("looptime = %f\n", t_temp1);
		}*/

		loopCounter++;
		*runningTime = *runningTime + stepsize;

		if (loopCounter == 10000 * L) {
			printf("Error report on Cycle %i \n", loopCounter);
			printf("runningTime = %f\n", *runningTime);
			L++;
		}
		
	}
	// "END OF MAIN LOOP" <<<<<<<<<<<<<<<<<<<<<<<<<< =<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	time_t endd_while = time(NULL);

	printf("End of Simulation Cycle, number of ADCS service calls = %li \n",
		   adcs_hk.EKF7a_callNr);
	printf("Time: %f \n", runningTime[0]);

	//Save variables to TEXT files: ///////////////////////////////////////////////////

	// Write running time to text file (MY ADDITION)
	FILE *run_time;
	run_time=fopen("..\\data\\output\\run_time.txt", "w");
		fprintf(run_time, "%f \n", runningTime[0]);
	fclose(run_time);
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	//HOUSEKEEP
	FILE *housefile;
	housefile = fopen("..\\data\\output\\housekeep.txt", "w");
	/* write text into the file stream*/

	for (unsigned int i = 0; i < k; i++) {
		fprintf(housefile, "--------------------------------------\n\
	ADCS housekeeping frame:              |\n\
	--------------------------------------\n\
	Time-stamp: %i-%i-%i-%i-%i-%i\n\n\
	The magnetic input vector:\n\n\
	%6.3f, %6.3f, %6.3f\n\n\
	The input vector to the sun:\n\n\
	%6.3f, %6.3f, %6.3f\n\n\
	The requested ADCS mode: %i\n\
	The given timestamp definition: %f\n\n\
	The calculated ECI position: \n\n\
	%5.8f, %5.8f, %5.8f\n\n\
	The calculated ECI velocity: \n\n\
	%5.8f, %5.8f, %5.8f\n\n\
	The calculated ECI vector to the sun: \n\n\
	%6.15f, %6.15f, %6.15f\n\n\
	The calculated ECI magnetic field vector: \n\n\
	%6.3f, %6.3f, %6.3f\n\n\
	The caclulated eclipse boolean: %i\n\n\
	The estimated quaternion: \n\n\
	%6.15f, %6.15f, %6.15f, %6.15f\n\n\
	The estimated angular rate: \n\n\
	%6.15f, %6.15f, %6.15f\n\n\
	The sum of the error covariance matrix: %6.15f\n\n\
	The estimated magnetic bias:\n\n\
	%6.3f, %6.3f, %6.3f\n\n",
								house_time[i][2], house_time[i][1], house_time[i][0],
								house_time[i][3], house_time[i][4], house_time[i][5],
								house_inmag[i][0], house_inmag[i][1], house_inmag[i][2],
								house_insol[i][0], house_insol[i][1], house_insol[i][2],
								house_inmode[i], house_indt[i],
								house_calcpos[i][0], house_calcpos[i][1], house_calcpos[i][2],
								house_calcvel[i][0], house_calcvel[i][1], house_calcvel[i][2],
								house_calcsol[i][0], house_calcsol[i][1], house_calcsol[i][2],
								house_calcmag[i][0], house_calcmag[i][1], house_calcmag[i][2],
								house_calceclipse[i],
								house_outq[i][0], house_outq[i][1], house_outq[i][2], house_outq[i][3],
								house_outvel[i][0], house_outvel[i][1], house_outvel[i][2],
								house_errorcovsum[i],
								house_outmagbias[i][0], house_outmagbias[i][1], house_outmagbias[i][2]);
		}
		fclose(housefile);

	// Write adcs modes switching times to text file (MY ADDITION)
	FILE *switch_time_vec;
	switch_time_vec=fopen("..\\data\\output\\switch_modes.txt", "w");
	for (unsigned int i = 0; i < modes_number; i++) {
		fprintf(switch_time_vec, "%6.3f \n",
			switch_modes[i]);
	}
	fclose(switch_time_vec);

	// Write TLE UPDATE times to text file (MY ADDITION)
	FILE *TLE_time_vec;
	TLE_time_vec=fopen("..\\data\\output\\switch_TLE.txt", "w");
	for (unsigned int i = 0; i < TLE_number; i++) {
		fprintf(TLE_time_vec, "%6.3f \n",
			switch_TLE[i]);
	}
	fclose(TLE_time_vec);

	// Write adcs parameters switching times to text file (MY ADDITION)
	FILE *ADCS_param_time_vec;
	ADCS_param_time_vec=fopen("..\\data\\output\\switch_ADCS_param.txt", "w");
	for (unsigned int i = 0; i < param_sets; i++) {
		fprintf(ADCS_param_time_vec, "%6.3f \n",
			switch_param[i]);
	}
	fclose(ADCS_param_time_vec);

	// Write ECI position to text file (MY ADDITION)
	FILE *r_ECI;
	r_ECI=fopen("..\\data\\output\\r_ECI.txt", "w");
	for (unsigned int i = 0; i < k; i++) {
		fprintf(r_ECI, "%5.8f, %5.8f, %5.8f \n",
			house_calcpos[i][0], house_calcpos[i][1], house_calcpos[i][2]);
	}
	fclose(r_ECI);

	// Write date 2 to text file (MY ADDITION)
	FILE *date2;
	date2=fopen("..\\data\\output\\date2.txt", "w");
	for (unsigned int i = 0; i < k; i++) {
		fprintf(date2,"%i-%i-%i-%i-%i-%i \n",
			date[i][0], date[i][1], date[i][2], date[i][3], date[i][4], date[i][5]);
	}
	fclose(date2);

	// Write datetime to text file (MY ADDITION)
	FILE *date_time;
	date_time=fopen("..\\data\\output\\datetime.txt", "w");
	for (unsigned int i = 0; i < k; i++) {
		fprintf(date_time, "%i-%i-%i-%i-%i-%i \n",
			house_time[i][2], house_time[i][1], house_time[i][0],
								house_time[i][3], house_time[i][4], house_time[i][5]);
	}
	fclose(date_time);

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	//Write estimated quaternion to text file
	FILE *quatE;
	quatE = fopen("..\\data\\output\\quatE.txt", "w");
	/* write text into the file stream*/

	for (unsigned int i = 0; i < k; i++) {
		fprintf(quatE, "%6.15f%s%6.15f%s%6.15f%s%6.15f \n",
			quat_est_mat[i][0], ",", quat_est_mat[i][1], ",",
			quat_est_mat[i][2], ",", quat_est_mat[i][3]);
	}
	fclose(quatE);

	//Write Omega Estimation to text file
	FILE *Om_est;
	Om_est = fopen("..\\data\\output\\OmEst.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Om_est, "%6.15f%s%6.15f%s%6.15f \n", Omega_mat[i][0],
			",", Omega_mat[i][1], ",", Omega_mat[i][2]);
	}
	fclose(Om_est);

	//Write Magnetic Field to text file
	FILE *B_file;
	B_file = fopen("..\\data\\output\\BE.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(B_file, "%6.15f%s%6.15f%s%6.15f \n", B_mat[i][0], ",",
			B_mat[i][1], ",", B_mat[i][2]);
	}
	fclose(B_file);

	//Write Sun Vector to text file
	FILE *uSun_file;
	uSun_file = fopen("..\\data\\output\\uSunE.txt", "w");

	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(uSun_file, "%6.15f%s%6.15f%s%6.15f \n", uSun_mat[i][0],
			",", uSun_mat[i][1], ",", uSun_mat[i][2]);
	}
	fclose(uSun_file);

	//Write Computed ISIS Sun Vector to text file
	FILE *uSun_ISIS_file;
	uSun_ISIS_file = fopen("..\\data\\output\\uSunE_ISIS.txt", "w");

	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(uSun_ISIS_file, "%6.15f%s%6.15f%s%6.15f \n", uSun_ISIS[i][0],
			",", uSun_ISIS[i][1], ",", uSun_ISIS[i][2]);
	}
	fclose(uSun_ISIS_file);

	// Write Eclipse flag from Matlab to text file
	FILE *Eclipse_Matlab;
	Eclipse_Matlab = fopen("..\\data\\output\\Eclipse_Matlab.txt", "w");

	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Eclipse_Matlab, "%6.15f \n", Eclipse_flag_output[i]);
	}
	fclose(Eclipse_Matlab);

	// // Write Eclipse flag from ISIS to text file
	// FILE *Eclipse_ISIS;
	// Eclipse_ISIS = fopen("..\\data\\output\\Eclipse_ISIS.txt", "w");

	// /* write text into the file stream*/
	// for (unsigned int i = 0; i < k; i++) {
	// 	fprintf(Eclipse_ISIS, "%6.15f \n", eclipseflag[i]);
	// }
	// fclose(Eclipse_ISIS);

	// Write Albedo_Threshold to text file
	FILE *Albedo_Threshold;
	Albedo_Threshold = fopen("..\\data\\output\\Albedo_Threshold.txt", "w");

	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Albedo_Threshold, "%6.15f \n", largest_element[i]);
	}
	fclose(Albedo_Threshold);

	//Write Comanded Dipole to text file
	FILE *Dcon_file;
	Dcon_file = fopen("..\\data\\output\\Dipcon.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Dcon_file, "%6.15f%s%6.15f%s%6.15f \n", Dip_mat[i][0],
			",", Dip_mat[i][1], ",", Dip_mat[i][2]);
	}
	fclose(Dcon_file);

	//Write Control Torque to text file
	FILE *Tcon_file;
	Tcon_file = fopen("..\\data\\output\\Tcon.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Tcon_file, "%6.15f%s%6.15f%s%6.15f \n", Tcon_mat[i][0],
			",", Tcon_mat[i][1], ",", Tcon_mat[i][2]);
	}
	fclose(Tcon_file);

	//Write real quaternion to text file
	FILE *quatR;
	quatR = fopen("..\\data\\output\\quatR_ECI.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(quatR, "%6.15f%s%6.15f%s%6.15f%s%6.15f \n",
			quat_real_mat[i][0], ",", quat_real_mat[i][1], ",",
			quat_real_mat[i][2], ",", quat_real_mat[i][3]);
	}
	fclose(quatR);

	FILE *quatR_B;
	quatR_B = fopen("..\\data\\output\\quatR_B.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(quatR_B, "%6.15f%s%6.15f%s%6.15f%s%6.15f \n",
			qReal_Body_mat[i][0], ",", qReal_Body_mat[i][1], ",",
			qReal_Body_mat[i][2], ",", qReal_Body_mat[i][3]);
	}
	fclose(quatR_B);

	// Write simulated angular velocity to text file
	FILE *Om_sym;
	Om_sym = fopen("..\\data\\output\\Om_sym.txt", "w");
	/* write text into the file stream*/
	for (unsigned int i = 0; i < k; i++) {
		fprintf(Om_sym, "%6.15f%s%6.15f%s%6.15f \n",
			Omega_sym_mat[i][0], ",", Omega_sym_mat[i][1], ",",
			Omega_sym_mat[i][2]);
	}
	fclose(Om_sym);

	time_t enddd = time(NULL); /* Stop timer */
	printf("Simulation duration: %f seconds \n",
		   (double)(enddd -
			strtt)); /* Calculate and display elapsed time */
	//printf("While-loop measured time: %f\n", (((et.tv_sec - st.tv_sec) * 1000000) + (et.tv_usec - st.tv_usec)) / 1000000.0);
	printf("While-loop duration: %f seconds \n",
		   (double)(endd_while -
			strt_while)); /* Calculate and display elapsed time */

	//  Release allocated memory by malloc functions- MY ADDITION

	free(adcs_param);
    free(switch_param);
    free(adcs_modes);
    free(switch_modes);
	free(switch_TLE);
    free(tlesets);
	return 0;
}

