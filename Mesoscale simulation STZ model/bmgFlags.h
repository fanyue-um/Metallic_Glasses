#ifndef H_BMGFLAGS
#define H_BMGFLAGS

/**************************************
  The set-up of simulation nature
  *************************************/
//#define COMPRESSION
#define TENSION
//#define SHEAR

/***************************************
  The applied residual stress can be
  either const or relaxed out. In addtion
  the system can be "pre-annealed".
 ***************************************/
//#define RESIDUAL_STRESS_CONST
//#define RESIDUAL_STRESS_RELX
//#define PRE_RELAXATION

/***************************************
  Boundary condition applied
  **************************************/
//#define STRESS_FREE_BOUNDARY

/***************************************
  Direction softening schemes.
  Except for SF_EVEN, all the other
  schemes have to be used in combination with
  the flag "MICRO_MEMORY"
 ***************************************/
//#define SF_LINEAR_FWD
//#define SF_LINEAR_BWD
//#define SF_EVEN
#define SF_SQUARE_FWD
//#define SF_Cubic_FWD
//#define SF_GAUSS_FWD

/****************************************
  Different systems 
 ****************************************/
#define LOW_T	/* general cases for low-T inhomogeneous */
//#define HIGH_T_Vit1
//#define HIGH_T_PdNiP
//#define RELX_BIMOLE_FLAG
//#define RELX_FVT_FLAG
//#define AVG_PERMSOFTENING

/***************************************
  Micro memory effect
  Right now only record the previous
  transformation modes and selected mode.
  No further effect on following transistions
 ***************************************/
#define MICRO_MEMORY

/***************************************
  Strain-rate sensitivity effect
  The LOCAL_RATE_V2 requires the recording
  of previous transformation
  **************************************/
//#define STRAIN_RATE_SENSITIVITY
//#define LOCAL_RATE_V2

/***************************************
  Inhomogeneous elasticity
 ***************************************/
#define INHOM_ELAST


/***************************************
  Dilatation effect
 ***************************************/
//#define STZ_DILATATION

/***************************************
  Subtract error to ensure stress boundary condition
 ***************************************/
//#define NO_SUBTRACT_ERROR


/***************************************
  Other functions and usages
 ***************************************/
//#define HEAT_CONDUCT
//#define XTAL_VARIANTS


/***************************************
  MRO related
 ***************************************/
#define NO_SF_XTAL_MRO
//#define FCC_SYMMETRY
//#define TWO_FOLD_SYMMETRY
//#define FULL_TENSOR_RECORD



#endif

