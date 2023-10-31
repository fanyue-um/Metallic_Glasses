/* my_fft.h */
#ifndef H_MY_FFT
#define H_MY_FFT

#include <fftw_mpi.h>
#include "bmg.h"

/* MPI parameters */
extern int rank, totnodes;
extern int lnx, lxs, lnyt, lyst, lsize;
extern MPI_Status status;
#ifdef BMG_GLOBALS_ONCE
int rank = 0, totnodes = 0;
int lnx = 0, lxs = 0, lnyt = 0, lyst = 0, lsize = 0;
MPI_Status status;
#endif

/* k-space grid */
extern real *g[DIM];
extern real *gg[3];
extern real *g2;
#ifdef BMG_GLOBALS_ONCE
real *g[DIM] = {NULL};
real *gg[3] = {NULL};
real *g2 = 0;
#endif

/* fft parameters */
extern fftwnd_mpi_plan plan, iplan;
extern fftw_complex *fft_data, *fft_work;
extern fftw_complex **EPSnSIGMA;
extern fftw_complex **Fiso;
extern fftw_complex *Giso;
#ifdef BMG_GLOBALS_ONCE
fftwnd_mpi_plan plan, iplan;
fftw_complex *fft_data = NULL, *fft_work = NULL;
fftw_complex **EPSnSIGMA = 0;
fftw_complex **Fiso = 0;
fftw_complex *Giso = 0;
#endif


#endif /* end my_fft.h */
