#ifndef H_BMG_MPI
#define H_BMG_MPI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw_mpi.h>
#include "bmgFlags.h"

//typedef float real;
typedef double real;
//#define MPI_real MPI_FLOAT
#define MPI_real MPI_DOUBLE
#define ABORT \
	MPI_Barrier(MPI_COMM_WORLD);	\
	MPI_Abort(MPI_COMM_WORLD, -3)

/* define system size */
#define NX	512
#define NY	512
#define NXY (NX*NY)
#define DIM 2
#define DE 3
/* potential energy landscape */
#define MODES 20
#define XMODES_SIX 6  // only 3 independent slip modes; the last three are negative of the first three
#define XMODES_TWO 2  // only 1 independent slip modes; the last three are negative of the first three
#define XMODES_TEN 12  // only 5 independent slip modes; the last three are negative of the first three


/* temperature */
#define T__K 300.
#if defined(HIGH_T_Vit1)
#define Tg__K 623.
#elif defined(HIGH_T_PdNiP)
#define Tg__K 587.5
#endif
#define TRIAL_FREQUENCY 1e13
#define E_gl 80.7
#define E_sc 30.0
#define k_E	39.0E-3
#define DeltaE	24.0
#define nu_gl 0.360
#define nu_sc 0.447

/* relaxation */
#ifdef RELX_BIMOLE_FLAG
#define QC_BIMOLE 1.85
//#define ETA_MAX_BIMOLE 0.223
#define ETA_EQI_BIMOLE 0.1
#define ETA_INT_BIMOLE 1.0
#endif

/* Lu et al. Acta Mater. 51 (2003) 3429 */
#ifdef RELX_FVT_FLAG
#define FVT_ALPHA 1.72E6
#define FVT_BETA 0.85
#define FVT_T0 426.3
#define FVT_VISCT_N (4.0E-5 * exp((16.5*FVT_T0)/(T__K-FVT_T0)))
#endif


/* units */
#define EV_IN_J 1.60217733e-19
#define J_IN_EV (1/EV_IN_J)
#define EPSILON 1e-17
#define BOLZ_IN_J__K 1.380658e-23
#define BOLZ BOLZ_IN_J__K
#define kT (T__K*BOLZ*J_IN_EV)
#define PI 3.14159265359
#define TWO_PI (2*PI)


/* voxel geometry */
#ifdef XTAL_VARIANTS
#define RECT_GRID_RATIO 1.73205081
#define LX 1.3
#define LY (1.3*RECT_GRID_RATIO)
#else
#define LX 0.7
#define LY 0.7
//#define LX 1.0
//#define LY 1.0
#endif
#define LZ 1.0
#define VOXEL0 (LX*LY*LZ)


#define space_loop for(i=0;i<lnx;i++)for(j=0;j<NY;j++)
#define GAMMA0 0.14
#ifdef STZ_DILATATION
#define DILATE0 (0.15*GAMMA0)
#endif
#define Sqr(x)	((x)*(x))

/* elastic medium */
extern real E, nu, lambda, mu, alpha, beta;	/* elastic moduli */
#ifdef BMG_GLOBALS_ONCE
real E = 0.0;
real nu = 0.0;
real mu = 0.0;
real lambda = 0.0;
real alpha = 0.0;
real beta = 0.0;
#endif

/* energetic parameters */
extern real Q0amp;
extern real Qmax;
extern real DW;
extern real XtalQ0amp;
extern real XtalShear;
extern real *Work[MODES];
extern real *Q0[MODES];
extern real *Q1[MODES];
extern real *Q2[MODES];
extern real *pQ1[MODES];
extern real *pQ2[MODES];
#ifdef MICRO_MEMORY
extern real *PreQ0[MODES];
#endif
extern real *Q[MODES];
extern real *tlast;
extern int softening_FLAG;	
#ifdef BMG_GLOBALS_ONCE
real Q0amp = 0.0;
real Qmax = 10.0;
real DW = 0.0;
real XtalQ0amp = 0.0;
real XtalShear = 0.0;
real *Work[MODES] = {NULL};
real *Q0[MODES] = {NULL};
real *Q1[MODES] = {NULL};
real *Q2[MODES] = {NULL};
real *pQ1[MODES] = {NULL};
real *pQ2[MODES] = {NULL};
#ifdef MICRO_MEMORY
real *PreQ0[MODES] = {NULL};
#endif
real *Q[MODES] = {NULL};
real *tlast = NULL;
int softening_FLAG = 0;
#endif

void ActEnergy(int step);
void SofteningFactor(void);
real AvgPermSoftening();
void ElasticEnergy(void);

/* Softening parameters */
extern real PermSoft, PermSoftmro, TempSoft, TempSoftmro, PermSofteningCap;
extern real *PermSoftening, *TempSoftening, *Softening;
extern int *MRO;   // MRO type 0-glass, 1-SixFold, 2-TwoFold
extern real *MROthe0;   // MRO shear direction (for one mode)
extern real xQ0amp;   // Activation energy for xtal-like MRO
extern real xMROFrac; // Volume fraction of xtal-like MRO
extern int InitMRO; // Initialization condition for MRO
extern real tTemp;
extern real tPerm;
extern real tResi;
#ifdef STRAIN_RATE_SENSITIVITY
extern real StrainRateSensitivity;
#endif
#ifdef BMG_GLOBALS_ONCE
real PermSoft = 0, PermSoftmro = 0, TempSoft=0, TempSoftmro = 0, PermSofteningCap=0;
real *PermSoftening = NULL, *TempSoftening = NULL, *Softening = NULL;
int *MRO = 0;   // MRO type
real *MROthe0 = 0;   // MRO shear direction (for one mode)
real xQ0amp = 0.;;   // Activation energy for xtal-like MRO
real xMROFrac = 0;; // Volume fraction of xtal-like MRO
int InitMRO = 0; // Initialization condition for MRO
real tTemp = 0.0;
real tPerm = 0.0;
real tResi = 0.0;
#ifdef STRAIN_RATE_SENSITIVITY
real StrainRateSensitivity = 0.;
#endif
#endif

/* local heating */
#define Cp 420.0	/* J/kg-K; Heat capacity */
#define rho 6125.0	/* kg/m^3; density */
#define thermal_diffusivity 3E-6	/* m^2/s */
#define thermal_coeff1	(1/Cp/rho)	/* K*m^3/J */
extern real *Tlocal;
#ifdef BMG_GLOBALS_ONCE
real *Tlocal = NULL;
#endif
void HeatConducting(real);

extern real strainrate;
extern real strainincrement;
extern real totime;
extern real dt;
extern real RandomRatio;
extern int printperiod;
extern real EndStrain;
extern real elEnergy;
extern real RelxTime;
#ifdef BMG_GLOBALS_ONCE
real strainrate = 0.0;
real RandomRatio = 0.0;
real strainincrement = 0.0;
real totime = 0;
real dt = 0.0;
int printperiod = 0;
real EndStrain = 0;
real elEnergy = 0.;
real RelxTime = 0.;
#endif

/* Initialize the simulation from the stream "in" */
void ReadInput(FILE * in);
void Init(void);
void FFT_MPI_Init(void);
void Grad(void);
void Energy_Init(void);
void Alloc(void);
void Free(void);
void ParameterPrint(void);
void RecordPrep(void);
void RecordFree(void);
void DispField(real *XR[2]);
void DispPrep(void);
void DispFree(void);
void WriteScalarMPI(real *field, char *s, int step);
void WriteScalarMPI_INT(int *field, char *s, int step);
void WriteMROMPI_INT(int *field, char *s, int step);


/* SS */
extern real *DEps0[MODES][DE];
#ifdef MICRO_MEMORY
extern real *PreDEps0[MODES][DE];
extern int	*PreSelectMode;
#endif
extern real *Eps0[DE];
extern real *Sigma[DE];
extern real *TotEps[DE];
extern real *RdStress[DE];
extern char ResStressDataFile[80];
extern real ResiRatio;
extern real Stress[3];
extern real EpsAvg[3];
extern real Epsyield;
extern real hisstress;
extern real *EqStrain;
extern real *EqStress;
extern real *strainhist;
extern real *stresshist;
extern real *elEnergyhist;
extern real *Relx_elEnergyhist;
extern real *XR[2];
#ifdef BMG_GLOBALS_ONCE
real *DEps0[MODES][DE] = {{NULL}};
#ifdef MICRO_MEMORY
real *PreDEps0[MODES][DE] = {{NULL}};
int *PreSelectMode = NULL;
#endif
real *Eps0[DE] = {NULL};
real *Sigma[DE] = {NULL};
real *TotEps[DE] = {0};
real *RdStress[DE] = {NULL};
char ResStressDataFile[80] = {0};
real ResiRatio = 0.;
real Stress[3] = {0};
real EpsAvg[3] = {0};
real Epsyield = 0.;
real hisstress = 0.;
real *EqStrain = NULL;
real *EqStress = NULL;
real *strainhist = NULL;
real *stresshist = NULL;
real *elEnergyhist = NULL;
real *Relx_elEnergyhist = NULL;
real *XR[2] = {0};
#endif

#ifdef INHOM_ELAST
extern real *Lambda, *Mu;
extern real *Tau0o[DE], *Tau1o[DE], *Tauo[DE];
extern real *epsolast[DE];
#ifdef BMG_GLOBALS_ONCE
real *Lambda=0, *Mu=0;
real *Tau0o[DE]={0};
real *Tau1o[DE]={0};
real *Tauo[DE]={0};
real *epsolast[DE]={0};
#endif

void HomogeneousStressSolver(real lambda, real mu, real **eps0o, real *epsavg,
		real **uo, real **epso, real **sigmao, real *stress, int uflag);
void InhomogeneousStressSolver(real *Lambdao, real *Muo, 
		real **eps0o, real *epsavg,
		real lambda, real mu,
		real **uo, real **epso, real **sigmao, real *stress, int uflag);
#endif

void DStrain(void);
void RdStrain(void);
void HomStressSolver(void);
void IsotropicMul(real lambda, real mu, fftw_complex *strain[], fftw_complex *stress[]);
void RIsotropicMul(real lambda, real mu, real *strain[]);
void athermal_plasticity(void);
void KMC(void);
void isotropicinv(real lambda, real mu, real *stress, real *strain);
void VoxelShearTransform(int, int);
void PreRelaxation(real RelxTime);
real relx_KMC(void);
void relx_athermal_plasticity(void);
void CheckStressSelfBalance(void);
void PrintResidualStress(void);
void ProjectResidualStress(char*,int);

/* Random generator parameters */
extern const gsl_rng_type *RT;
extern gsl_rng *r;
extern int seed;
#ifdef BMG_GLOBALS_ONCE
const gsl_rng_type *RT;
gsl_rng *r;
int seed = 0;
#endif

void Rand_Init(void);

/* counter for SOC analysis */
extern int *countSOC[2];
extern real SOCelast;
extern int sizeSOC;
extern real *elastSOC[2];
extern int sizeElaSOC;
extern int avlanCounter;
extern real *wtimeCounter;
extern real wtimeStart;
extern real wtimeEnd;
#ifdef BMG_GLOBALS_ONCE
int *countSOC[2] = {0};
real SOCelast = 0.;
int sizeSOC = 0;
real *elastSOC[2] = {0};
int sizeElaSOC = 0;
int avlanCounter = 0;
real *wtimeCounter = 0;
real wtimeStart = 0.;
real wtimeEnd = 0.;
#endif
void UpdateSOC(void);

/* step counter */
extern int step;
extern int astep;
extern int kstep;
extern int estep;
extern int flag, aflag;
#ifdef BMG_GLOBALS_ONCE
int step = 0;
int astep = 0;
int kstep = 0;
int estep = 0;
int flag = 0;
int aflag = 0;
#endif

/* record actual activaion barriers */
extern FILE *fq;
#ifdef BMG_GLOBALS_ONCE
FILE *fq;
#endif

void PrintFlags(void);

#endif /* end bmg_mpi.h */
