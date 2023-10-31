/* init.c */
#include "bmg.h"
#include "my_fft.h"
#include <fftw_mpi.h>

void Init(void){

	/* for 2D */
#ifdef XTAL_VARIANTS
	/* material parameters from single crystal Cu */
	nu = 0.3;
	mu = 40;
	E = mu*(2+2*nu);
//	nu = nu * 2;
	/* 2D isotropic elastic medium */
//	mu = mu - (4E-3)*T__K;
//	E = 114.36 - (6.12E-2)*T__K;
#elif defined(LOW_T)
//	mu = 37.15 - 2.2090/(exp(252.9/T__K)-1);	/* Vashni expression */
//	E = 100.48 - 5.7992/(exp(259.3/T__K)-1);	/* Vashni expression */
//	nu = 0.3531 + 0.0041/(exp(209.1/T__K)-1);	/* Vashni expression */

//	Q0amp = 5;
//	Q0amp = Q0amp*(1+T__K/10000.0);	/* inverted Meyer-Neldel rule for low T */
	//Q0amp = Q0amp*(1-T__K/(Tg__K+100));	/* Meyer-Neldel rule for high T */
	nu *=2;
	mu = E/(2+2*nu);
#elif defined(HIGH_T_Vit1)
	E = (E_gl + E_sc)/2. - (E_gl - E_sc)*tanh((T__K-Tg__K)/DeltaE)/2.
		- k_E*(T__K-Tg__K);
	nu = (nu_gl + nu_sc)/2. - (nu_gl - nu_sc)*tanh((T__K-Tg__K)/DeltaE)/2.;
	mu = E / (2 + 2 * nu);

//	Q0amp = (Q_gl + Q_sc)/2. - (Q_gl - Q_sc)*tanh((T__K - Tg__K)/DeltaQ)/2.;
	Q0amp = 2.2;
#elif defined(HIGH_T_PdNiP)
	E = 96;
	nu = 0.36;
	mu = E / (2 + 2*nu);

	Q0amp = 3.0;
#endif

	lambda = E * nu / (1 - nu * nu);
	alpha = (1 + nu) / 2;
	beta = nu;

	dt = strainincrement/strainrate;
	tTemp = 1/(TRIAL_FREQUENCY * exp(-0.37/kT));
//	tPerm = 1/(TRIAL_FREQUENCY * exp(-1.7/kT));

//	tTemp = 1/(TRIAL_FREQUENCY * exp(-1.45/kT));
//	tPerm = 1/(TRIAL_FREQUENCY * exp(-1.85/kT));
	tResi = 1/(TRIAL_FREQUENCY * exp(-1.80/kT));
#if defined(RELX_FVT_FLAG)
	Q0amp *= (1-exp(-pow(FVT_ALPHA/strainrate/FVT_VISCT_N, FVT_BETA)));
#endif

	/* Damage trackers */
	PermSoft = 13;
	PermSoftmro = 20;
	TempSoft = 39;
	TempSoftmro = 60;
//	PermSoft = 20.;
//	TempSoft = 60.;
//	PermSoft = 5;
//	TempSoft = 15;
	PermSofteningCap = -log(0.6);

#ifdef STRAIN_RATE_SENSITIVITY
	StrainRateSensitivity = 0.1316;
	//StrainRateSensitivity = 0.5;
#endif
	return;
}/* end Init() */

void FFT_MPI_Init(void){
  int i;

	plan = fftw2d_mpi_create_plan(MPI_COMM_WORLD, NX, NY,
				FFTW_FORWARD, FFTW_ESTIMATE);
	iplan = fftw2d_mpi_create_plan(MPI_COMM_WORLD, NX, NY,
				FFTW_BACKWARD, FFTW_ESTIMATE);

	/* slab decomposition */
	fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);

	/* malloc data and work vectors */
	fft_data = (fftw_complex *)fftw_malloc(lsize * sizeof(fftw_complex));
	if(!fft_data) printf("malloc failed for fft_data\n");
	fft_work = (fftw_complex *)fftw_malloc(lsize * sizeof(fftw_complex));
	if(!fft_work) printf("malloc failed for fft_work\n");
	EPSnSIGMA = (fftw_complex **)fftw_malloc(3*sizeof(fftw_complex*));
	if(!EPSnSIGMA) printf("malloc failed for EPSnSIGMA\n");
	for(i=0;i<3;i++){
		EPSnSIGMA[i] = (fftw_complex *)fftw_malloc(lsize*sizeof(fftw_complex));
		if(!EPSnSIGMA[i]) printf("malloc failed for EPSnSIGMA[%d]\n",i);
	}
	Fiso = (fftw_complex **)fftw_malloc(2*sizeof(fftw_complex*));
	if(!Fiso) printf("malloc failed for Fiso\n");
	for(i=0;i<2;i++){
		Fiso[i] = (fftw_complex *)fftw_malloc(lsize*sizeof(fftw_complex));
		if(!Fiso[i]) printf("malloc failed for Fiso[%d]\n",i);
	}
	Giso = (fftw_complex *)fftw_malloc(lsize*sizeof(fftw_complex));
	if(!Giso) printf("malloc failed for Giso\n");


	return;
}/* end FFT_MPI_Init() */

void Grad(void){
	int i, j, ir, idx;
	real s1, s2;

	for(i=0; i<DIM; i++) g[i] = (real *)malloc(lnx*NY * sizeof(real));
	for(i=0; i<3; i++) gg[i] = (real *)malloc(lnx*NY * sizeof(real));
	g2 = (real *)malloc(lnx*NY * sizeof(real));

	space_loop{
		ir = i+lxs;
		idx = i*NY + j;
#ifdef XTAL_VARIANTS
		s1 = (real)(j) * (real)(1.0/NY); s2 = (real)(ir) * (real)(1.0/(NX*RECT_GRID_RATIO));
#else
		s1 = (real)(j) * (real)(1.0/NY); s2 = (real)(ir) * (real)(1.0/NX);
#endif
		g[0][idx] = (s1 - round(s1)) * TWO_PI;
		g[1][idx] = (s2 - round(s2)) * TWO_PI;
		g2[idx] = sqrt(g[0][idx]*g[0][idx] + g[1][idx]*g[1][idx]);
		if((ir+j) == 0) g2[idx] = 1.0;
		g[0][idx] /= g2[idx];
		g[1][idx] /= g2[idx];
		gg[0][idx] = g[0][idx] * g[0][idx];
		gg[1][idx] = g[1][idx] * g[1][idx];
		gg[2][idx] = g[0][idx] * g[1][idx];
	}

	return;
}/* end Grad() */

void Rand_Init(void){
	seed = seed + rank*lsize*4;
	gsl_rng_default_seed = seed;

	RT = gsl_rng_default;
	r = gsl_rng_alloc(RT);

	return;
}/* end Rand_Init() */

void Energy_Init(void){
	int i, j, m, idx;

	
//		space_loop{
//			idx = i*NY + j;
/*      if(MRO[idx]==1){  // six-fold
        if(m<XMODES_SIX){
			    Q0[m][idx] = xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
        }else{
          Q0[m][idx] = INFINITY;
        }
      }
      else if(MRO[idx]==2){  // two-fold
        if(m<XMODES_TWO){
			    Q0[m][idx] = xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
        }else{
          Q0[m][idx] = INFINITY;
        }
      }
      else if(MRO[idx]==3){  // ten-fold
        if(m<XMODES_TEN){
			    Q0[m][idx] = xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
        }else{
          Q0[m][idx] = INFINITY;
        }
      }
      else{
 */   	 
	for(m=0; m<MODES; m++){
		space_loop{
			idx = i*NY + j;
	//   			if(gsl_rng_uniform(r)>=0.3){
					Q0[m][idx]=gsl_ran_gaussian(r, 0.25)+Q0amp;
			//		if(Q0[m][idx]<0){
			//			Q0[m][idx]=0;
			//		}
					Q1[m][idx]=Q0[m][idx];
	//				pQ1[m][idx]=1;
	//			}
	/*			else{
					Q0[m][idx]=gsl_ran_gaussian(r, 0.25)+Q0amp;
					Q1[m][idx]=Q0[m][idx];
					pQ1[m][idx]=2;
				}
	*/	}
	}
		//	  Q0[m][idx] = Q0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
   //   }
//		}
//	}


	return;
}/* end Energy_Init() */
double gsl_ran_gaussianr (const gsl_rng * r, const double sigma){
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x = -1 + 2 * gsl_rng_uniform_pos (r);
      y = -1 + 2 * gsl_rng_uniform_pos (r);

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

void Alloc(void){
	int i, j, m, de, idx;

	for(m=0; m<MODES; m++){
		Q0[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Q0[m]){printf("malloc failed for Q0\n");}
		Q1[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Q1[m]){printf("malloc failed for Q1\n");}
		Q2[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Q2[m]){printf("malloc failed for Q2\n");}
		pQ1[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!pQ1[m]){printf("malloc failed for pQ1\n");}
		pQ2[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!pQ2[m]){printf("malloc failed for pQ2\n");}

  // Allocate disp field array
	DispPrep();

#ifdef MICRO_MEMORY
		PreQ0[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!PreQ0[m]){printf("malloc failed for PreQ0\n");}
#endif

		Q[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Q[m]){printf("malloc failed for Q\n");}
		Work[m] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Work[m]){printf("malloc failed for Work\n");}

		for(de=0; de<DE; de++){
			DEps0[m][de] = (real *)malloc(lnx*NY * sizeof(real));
			if(!DEps0[de]){printf("malloc failed for DEps0\n");}
#ifdef MICRO_MEMORY
			PreDEps0[m][de] = (real *)malloc(lnx*NY * sizeof(real));
			if(!PreDEps0[de]){printf("malloc failed for DEps0\n");}
#endif
		}
	}

	for(de=0; de<DE; de++){
		Eps0[de] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Eps0[de]){printf("malloc failed for Eps0\n");}

		Sigma[de] = (real *)malloc(lnx*NY * sizeof(real));
		if(!Sigma[de]){printf("malloc failed for Sigma\n");}

		TotEps[de] = (real *)malloc(lnx*NY * sizeof(real));
		if(!TotEps[de]){printf("malloc failed for TotEps\n");}

		RdStress[de] = (real *)malloc(lnx*NY * sizeof(real));
		if(!RdStress[de]){printf("malloc failed for RdStress\n");}
	}

	tlast = (real *)malloc(lnx*NY * sizeof(real));
	if(!tlast){printf("malloc failed for tlast\n");}

#ifdef MICRO_MEMORY
	PreSelectMode = (int *)malloc(lnx*NY * sizeof(int));
	if(!tlast){printf("malloc failed for PreSelectMode\n");}
#endif


	/* initialize Sigma, Eps0, RdStress */
	for(de=0; de<DE; de++){
		space_loop{
			m = i*NY + j;
			Eps0[de][m] = 0.;
			Sigma[de][m] = 0.;
      TotEps[de][m] = 0.;
			RdStress[de][m] = 0.;
		}
	}

	/* initialize PreDEps0 */
#ifdef MICRO_MEMORY
	for(m=0; m<MODES; m++){
		for(de=0; de<DE; de++){
			space_loop{
				idx = i*NY + j;
				PreDEps0[m][de][idx] = 0.;
			}
		}
	}
#endif

	/* initialize tlast */
	space_loop{
		m = i*NY + j;
		tlast[m] = -100*tTemp;
#ifdef MICRO_MEMORY
		PreSelectMode[m] = 0;
#endif
	}


	PermSoftening = (real *)malloc(lnx*NY * sizeof(real));
	if(!PermSoftening){printf("malloc failed for PermSofteing\n");}
	/* initialize PermSoftening */
	space_loop{
		m = i*NY + j;
		PermSoftening[m] = 0;
	}


	TempSoftening = (real *)malloc(lnx*NY * sizeof(real));
	if(!TempSoftening){printf("malloc failed for TempSoftening\n");}
	/* initialize TempSoftening */
	space_loop{
		m = i*NY + j;
		TempSoftening[m] = 0;
	}

	
	Softening = (real *)malloc(lnx*NY * sizeof(real));
	if(!Softening){printf("malloc failed for Softening\n");}
	/* initialize Softening */
	space_loop{
		m = i*NY + j;
		Softening[m] = 0;
	}

	MRO = (int *)malloc(lnx*NY * sizeof(int));
  if(!MRO){printf("malloc failed for MRO\n");}
	MROthe0 = (real *)malloc(lnx*NY * sizeof(real));
  if(!MROthe0){printf("malloc failed for MROthe0\n");}
  /* initialize MRO*/
  if(InitMRO==0){
	int n;
	n=0;
//	FILE *f3;
//	f3 = fopen("particle.out", "w");
    	  space_loop{
    	    m = i*NY+j;
//    	if(rank==3){
	  //  	m=32767;
	  //  	m=32512;
	    real tmp = gsl_rng_uniform(r);
    	    if(tmp<=xMROFrac){ // xtal-like
    		n++;
	//	int d,the0;
/*		d=10;
	//	d=20*(0.75+0.5*gsl_rng_uniform(r));
//		fprintf(f3, "%d\t%d\t%d\n", i, j, d);
        	the0 = gsl_rng_uniform(r)*PI-PI/2.0;  // [-PI/2, PI/2]
		int i1,j1;
		for(i1=0; i1<d; i1++){
			for(j1=0; j1<d; j1++){
				if((i1-d/2)*(i1-d/2)+(j1-d/2)*(j1-d/2)<d*d/4 && (m-(i1)*NY-j1)>=0){
				MRO[m-(i1)*NY-j1] = 2;	
        MROthe0[m-(i1)*NY-j1] = the0;  // [-PI/2, PI/2]
    //    MROthe0[m-(i1)*NY-j1] = PI/1.0;  // [-PI/2, PI/2]
				}
			}
		}
  */	       MRO[m] = 1; // six-fold
/*	    if(m>1){
		int i1,j1;
		for(i1=0; i1<30; i1++){
			for(j1=0; j1<30; j1++){
				MRO[m-i1*NY-j1] = 1;
			}
		}
	}
*/	      //MRO[m] = 2; // two-fold
	      //MRO[m] = 3; // ten-fold
        MROthe0[m] = gsl_rng_uniform(r)*PI-PI/2.0;  // [-PI/2, PI/2]
     //     MROthe0[m] = PI/2.0;  // [-PI/2, PI/2]
	    }
	    else{ // normal glass
	      MRO[m] = 0;
        MROthe0[m] = 0.;
	    }
	  }
  }
    else if(InitMRO==1){
    FILE *fp;
    int ii,jj;
    fp = fopen("MROstructure.ms","r");
    // empty reading to go to the corresponding slabbed region
    int EmptySteps = rank*NXY/totnodes;
    char buffer[80] = {0};
    for(int i=0;i<EmptySteps;++i){
      fgets(buffer,80,fp);
    }
    space_loop{
	     m = i*NY+j;
       fgets(buffer,80,fp); 
       sscanf(buffer,"%d %d %d %lf",&ii,&jj,&MRO[m],&MROthe0[m]);
       idx = (ii-1-lxs)*NY+jj-1;
       if(idx!=m){
	        perror("Error in reading MROstructure.ms file (inconsistent index system)!!\n");
		      MPI_Abort(MPI_COMM_WORLD, -1);
       }
    }

#if 1
    space_loop{
	     m = i*NY+j;
       if(MRO[m]==3){ // change to a mix of ten-fold and two-fold
	       real tmp = gsl_rng_uniform(r);
         if(tmp<0.5) MRO[m] = 2;
       }
       //if(MRO[m]==3){ // change ten-fold to two-fold
        // MRO[m] = 2;
       //}
    }
#endif

  }
  else{
	  perror("Error: wrong option for MRO initialization\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
  }
  WriteMROMPI_INT(MRO, "Initial_MROstructure",0);

#ifdef INHOM_ELAST
	Lambda = (real *)malloc(lsize * sizeof(real));
	Mu = (real *)malloc(lsize * sizeof(real));
  for(de=0;de<DE;de++){
		Tau0o[de] = (real *)malloc(lsize * sizeof(real));
		if(!Tau0o[de]){printf("malloc failed for Tau0o\n");}

		Tau1o[de] = (real *)malloc(lsize * sizeof(real));
		if(!Tau1o[de]){printf("malloc failed for Tau1o\n");}

		Tauo[de] = (real *)malloc(lsize * sizeof(real));
		if(!Tauo[de]){printf("malloc failed for Tauo\n");}

		epsolast[de] = (real *)malloc(lsize * sizeof(real));
		if(!epsolast[de]){printf("malloc failed for epsolast\n");}
  }

  // Inhomogeneous elastic medium
  space_loop{
		m = i*NY + j;
    //real tmp = 0.8 + (!MRO[m])*0.2;
    //real tmp = 1.0 + (int)(MRO[m]>0)*0.2;
    //real tmp = 1.0 + (int)(MRO[m]>0)*0.5;
    real tmp = 1.0; 
    Lambda[m] = tmp*lambda;
    Mu[m] = tmp*mu;
  }

#endif


	Tlocal = (real *)malloc(lnx*NY * sizeof(real));
	if(!Tlocal){printf("malloc failed for Tlocal\n");}
	/* initialize local temperature */
	space_loop{
		m = i*NY + j;
		Tlocal[m] = T__K;
	}

	return;
}/* end Alloc() */


void Free(void){
	int i, m, de;

	fftwnd_mpi_destroy_plan(plan);
	fftwnd_mpi_destroy_plan(iplan);
	fftw_free(fft_data);
	fftw_free(fft_work);
	for(i=0;i<3;i++)
		fftw_free(EPSnSIGMA[i]);
	free(EPSnSIGMA);
	for(i=0;i<2;i++)
		fftw_free(Fiso[i]);
	free(Fiso);
	free(Giso);

	gsl_rng_free(r);

	DispFree();
	
	for(m=0; m<MODES; m++){
		free(Q0[m]);
		free(Q1[m]);
#ifdef MICRO_MEMORY
		free(PreQ0[m]);
#endif
		free(Q[m]);
		for(de=0; de<DE; de++){
			free(DEps0[m][de]);
#ifdef MICRO_MEMORY
			free(PreDEps0[m][de]);
#endif
		}
	}

	for(de=0; de<DE; de++){
		free(Eps0[de]);
		free(Sigma[de]);
		free(TotEps[de]);
		free(RdStress[de]);
	}

	for(de=0; de<DIM; de++) free(g[de]);
	for(de=0; de<3; de++) free(gg[de]);
	free(g2);

	free(tlast);
#ifdef MICRO_MEMORY
	free(PreSelectMode);
#endif

#ifdef INHOM_ELAST
  free(Lambda);
  free(Mu);
  for(de=0;de<DE;de++){
    free(Tau0o[de]);
    free(Tau1o[de]);
    free(Tauo[de]);
    free(epsolast[de]);
  }
#endif

	free(PermSoftening);
	free(TempSoftening);
	free(Softening);
  free(MRO);
	free(Tlocal);

	free(strainhist);
	free(stresshist);
	free(elEnergyhist);

	for(de=0; de<2; de++){
		free(countSOC[de]);
		free(elastSOC[de]);
	}

	return;
}/* end Free() */

void ActEnergy(int step){
	int i, j, m, idx, mm, temp1;
	real DWork;
	real DotProdt, norm1, norm2, DireModifier;
//	printf("act1");
//				printf("EpsAvg2=%.6e\n", EpsAvg[2]);
//				printf("Stress0=%.6e\n", Stress[0]);
	for(m=0; m<MODES; m++)
		space_loop{
			idx = i*NY + j;

#ifdef MICRO_MEMORY
			mm = PreSelectMode[idx];
#endif

#if defined(RESIDUAL_STRESS_CONST)
			DWork = ((Sigma[0][idx]+RdStress[0][idx]) * DEps0[m][0][idx] + 
					(Sigma[1][idx]+RdStress[1][idx]) * DEps0[m][1][idx] +
					2 * (Sigma[2][idx]+RdStress[2][idx]) * DEps0[m][2][idx]) * VOXEL0 * 1e-18 * J_IN_EV;
#elif defined(RESIDUAL_STRESS_RELX)
			DWork = ((Sigma[0][idx]+RdStress[0][idx]*exp(-totime/tResi)) * DEps0[m][0][idx] +
					(Sigma[1][idx]+RdStress[1][idx]*exp(-totime/tResi)) * DEps0[m][1][idx] +
					2 * (Sigma[2][idx]+RdStress[2][idx]*exp(-totime/tResi)) * DEps0[m][2][idx]) * VOXEL0 * 1e-18 * J_IN_EV;
#else
			DWork = (Sigma[0][idx] * DEps0[m][0][idx] + Sigma[1][idx] * DEps0[m][1][idx] +
					2 * Sigma[2][idx] * DEps0[m][2][idx]) * VOXEL0 * 1e-18 * J_IN_EV;
#endif

			/* different softening scheme */
#if defined(SF_LINEAR_FWD)
			DotProdt = PreDEps0[mm][0][idx] * DEps0[m][0][idx] + PreDEps0[mm][1][idx] * DEps0[m][1][idx] +
					2 * PreDEps0[mm][2][idx] * DEps0[m][2][idx];
			norm1 = sqrt(DEps0[m][0][idx] * DEps0[m][0][idx] + DEps0[m][1][idx] * DEps0[m][1][idx] +
					2 * DEps0[m][2][idx] * DEps0[m][2][idx]);
			norm2 = sqrt(PreDEps0[mm][0][idx] * PreDEps0[mm][0][idx] + PreDEps0[mm][1][idx] * PreDEps0[mm][1][idx] +
					2 * PreDEps0[mm][2][idx] * PreDEps0[mm][2][idx]);

			DireModifier = DotProdt / (norm1 + EPSILON) / (norm2 + EPSILON);
			DireModifier = (DireModifier + 1)/2;	/* forward preferred */
#elif defined(SF_LINEAR_BWD)
			DotProdt = PreDEps0[mm][0][idx] * DEps0[m][0][idx] + PreDEps0[mm][1][idx] * DEps0[m][1][idx] +
					2 * PreDEps0[mm][2][idx] * DEps0[m][2][idx];
			norm1 = sqrt(DEps0[m][0][idx] * DEps0[m][0][idx] + DEps0[m][1][idx] * DEps0[m][1][idx] +
					2 * DEps0[m][2][idx] * DEps0[m][2][idx]);
			norm2 = sqrt(PreDEps0[mm][0][idx] * PreDEps0[mm][0][idx] + PreDEps0[mm][1][idx] * PreDEps0[mm][1][idx] +
					2 * PreDEps0[mm][2][idx] * PreDEps0[mm][2][idx]);

			DireModifier = DotProdt / (norm1 + EPSILON) / (norm2 + EPSILON);
			DireModifier = (-DireModifier + 1)/2;	/* backward prefered */
#elif defined(SF_EVEN)
			DireModifier = 0.5;						/* even prefered */
#elif defined(SF_SQUARE_FWD)
			DotProdt = PreDEps0[mm][0][idx] * DEps0[m][0][idx] + PreDEps0[mm][1][idx] * DEps0[m][1][idx] +
					2 * PreDEps0[mm][2][idx] * DEps0[m][2][idx];
			norm1 = sqrt(DEps0[m][0][idx] * DEps0[m][0][idx] + DEps0[m][1][idx] * DEps0[m][1][idx] +
					2 * DEps0[m][2][idx] * DEps0[m][2][idx]);
			norm2 = sqrt(PreDEps0[mm][0][idx] * PreDEps0[mm][0][idx] + PreDEps0[mm][1][idx] * PreDEps0[mm][1][idx] +
					2 * PreDEps0[mm][2][idx] * PreDEps0[mm][2][idx]);

			DireModifier = DotProdt / (norm1 + EPSILON) / (norm2 + EPSILON);
			DireModifier = (1 + DireModifier)*(1 + DireModifier)/4;		/* SquareModified; forward */
#elif defined(SF_CUBIC_FWD)
			DotProdt = PreDEps0[mm][0][idx] * DEps0[m][0][idx] + PreDEps0[mm][1][idx] * DEps0[m][1][idx] +
					2 * PreDEps0[mm][2][idx] * DEps0[m][2][idx];
			norm1 = sqrt(DEps0[m][0][idx] * DEps0[m][0][idx] + DEps0[m][1][idx] * DEps0[m][1][idx] +
					2 * DEps0[m][2][idx] * DEps0[m][2][idx]);
			norm2 = sqrt(PreDEps0[mm][0][idx] * PreDEps0[mm][0][idx] + PreDEps0[mm][1][idx] * PreDEps0[mm][1][idx] +
					2 * PreDEps0[mm][2][idx] * PreDEps0[mm][2][idx]);

			DireModifier = DotProdt / (norm1 + EPSILON) / (norm2 + EPSILON);
			DireModifier = (DireModifier + 1)*(DireModifier + 1)*(DireModifier + 1)/8;	/* CubicModified forward */
#elif defined(SF_GAUSS_FWD)
			DotProdt = PreDEps0[mm][0][idx] * DEps0[m][0][idx] + PreDEps0[mm][1][idx] * DEps0[m][1][idx] +
					2 * PreDEps0[mm][2][idx] * DEps0[m][2][idx];
			norm1 = sqrt(DEps0[m][0][idx] * DEps0[m][0][idx] + DEps0[m][1][idx] * DEps0[m][1][idx] +
					2 * DEps0[m][2][idx] * DEps0[m][2][idx]);
			norm2 = sqrt(PreDEps0[mm][0][idx] * PreDEps0[mm][0][idx] + PreDEps0[mm][1][idx] * PreDEps0[mm][1][idx] +
					2 * PreDEps0[mm][2][idx] * PreDEps0[mm][2][idx]);

			DireModifier = DotProdt / (norm1 + EPSILON) / (norm2 + EPSILON);
			DireModifier = exp(-DireModifier*DireModifier*9/2);		/* GaussModified forward */
#endif

      /* Consider if the xtal-like MRO should experience softening */
#ifdef NO_SF_XTAL_MRO
      if(MRO[idx]==1){  // six-fold
         DireModifier *= 0.1;
         Softening[idx] *=0.1;
      }
      else if(MRO[idx]==2){ // two-fold
         DireModifier *= 0.01;
         Softening[idx] *=0.01;
      }
      else if(MRO[idx]==3){ // ten-fold
         DireModifier *= 0.9;
         Softening[idx] *= 0.9;
      }
      else if(MRO[idx]==0){ // glass
        /* do nothing */
      }
      else{ // invalid phase
	printf("%d\n",idx);
	printf("%d\n",MRO[idx]);
	      perror("Error 2 in MRO type (only {0,1,2,3} are allowed)!!\n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
      }
#endif

#ifdef STRAIN_RATE_SENSITIVITY
			real tmpEps;
			real tmpRate;
#ifdef LOCAL_RATE_V2
			tmpEps = sqrt(PreDEps0[mm][0][idx]*PreDEps0[mm][0][idx]+PreDEps0[mm][1][idx]*PreDEps0[mm][1][idx]+
					2*PreDEps0[mm][2][idx]*PreDEps0[mm][2][idx]);
#else
			tmpEps = sqrt(Eps0[0][idx]*Eps0[0][idx]+Eps0[1][idx]*Eps0[1][idx]+
					2*Eps0[2][idx]*Eps0[2][idx]);
#endif
			tmpRate = 1.0+tmpEps/(totime-tlast[idx]+EPSILON)/strainrate;
			Q[m][idx] = Q0[m][idx]*exp(-DireModifier*Softening[idx])*pow(tmpRate,
					StrainRateSensitivity)- DWork/2;
#else
		
			Q[m][idx] = Q0[m][idx]*exp(-DireModifier*Softening[idx]) - DWork/2;
			Work[m][idx] = DWork/2;
#endif
		}
		
	return;
}/* end ActEnergy() */

real AvgPermSoftening(void){
	real sum_loc, avg_all;
	int i, j, idx;

	sum_loc = 0.0;
	avg_all = 0.0;
	space_loop{
		idx = i*NY + j;
		sum_loc += PermSoftening[idx];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&sum_loc, &avg_all, 1, MPI_DOUBLE,
			MPI_SUM, MPI_COMM_WORLD);

	avg_all /= (real)NXY;
	
	return avg_all;
}/* end AvgPermSoftening() */

void SofteningFactor(void){
	int i, j, idx;
	real PermAvg;

	space_loop{
		idx = i*NY + j;
#if defined(RELX_BIMOLE_FLAG)
		//Softening[idx] = ETA_MAX_BIMOLE / (1.+(ETA_MAX_BIMOLE/ETA_INT_BIMOLE-1.)*exp(-totime/tPerm));
		//Softening[idx] = ETA_EQI_BIMOLE / (1.+(ETA_EQI_BIMOLE/ETA_INT_BIMOLE-1.)*exp(-totime/tPerm));
		Softening[idx] = ETA_EQI_BIMOLE / (1.+(ETA_EQI_BIMOLE/(PermSoftening[idx]+EPSILON)-1.)*exp((tlast[idx]-totime)/tPerm));
#elif defined(RELX_FVT_FLAG)
		//Softening[idx] = totime*PermSoftening[idx]/(1-exp(-pow(FVT_ALPHA*(totime-tlast[idx])/FVT_VISCT_N, FVT_BETA)))
		//		+ TempSoftening[idx]*exp((tlast[idx]-totime)/tTemp);
		Softening[idx] = PermSoftening[idx]*exp((tlast[idx]-totime)/tPerm)	+ TempSoftening[idx]*exp((tlast[idx]-totime)/tTemp);
#elif defined(AVG_PERMSOFTENING)
		PermAvg = AvgPermSoftening();
		Softening[idx] = PermAvg + TempSoftening[idx] * exp((tlast[idx] - totime)/tTemp);
		//Softening[idx] = PermAvg + TempSoftening[idx] * exp(( - totime)/tTemp);
		//Softening[idx] = PermAvg*exp((tlast[idx]-totime)/tTemp);
#else
		Softening[idx] = PermSoftening[idx] + TempSoftening[idx] * exp((tlast[idx] - totime)/tTemp);
		//Softening[idx] = PermSoftening[idx]*exp((tlast[idx]-totime)/tPerm) + TempSoftening[idx] * exp((tlast[idx] - totime)/tTemp);
		//Softening[idx] = TempSoftening[idx] * exp((tlast[idx] - totime)/tTemp);
//		Softening[idx] = PermSoftening[idx]*tanh((totime)/tPerm) + TempSoftening[idx] * exp((tlast[idx] - totime)/tTemp);
#endif
	}
	
	return;
}/* end SofteningFactor() */

void PreRelaxation(real TimeLimit){
	int i,j,m,de;
	int idx;
	int step_relx = 0;
	real time_relx = 0.;
	real *time_relx_hist = 0;
	real *Relx_elEnergyhist = 0;
	real *Relx_strainhist[3] = {0};
	real *Relx_stresshist[3] = {0};
	real elEnergyAll = 0.;
	FILE *fp;
	char fname[20] = {0};
	int tmp_astep, tmp_kstep;

	tmp_astep = 0;
	tmp_kstep = 0;

	while(1){
		step_relx++;
		elEnergyAll = 0.;

		/* for relaxation, no softening! */
		/* Softening factor */
//		if(softening_FLAG==1) SofteningFactor();

		/* Activation Energy */
		ActEnergy(step_relx);
		
		/* Sorting activation energy */
		flag = 0;	aflag = 0;
		for(m=0; m<MODES; m++)
			space_loop{
				idx = i*NY + j;
				if(Q[m][idx] < 0) flag++;
			}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&flag, &aflag, 1, MPI_INT,
				MPI_SUM, MPI_COMM_WORLD);

		
		/* evolution */
		if(aflag >0){
			/* athermal plasticity */
			MPI_Barrier(MPI_COMM_WORLD);
			relx_athermal_plasticity();
			tmp_astep++;
		}
		else{
			/* thermal */
			MPI_Barrier(MPI_COMM_WORLD);
			time_relx += fabs(relx_KMC());
			tmp_kstep++;
		}

		/* store mean stress and strain */
		for(i=0; i<3; i++){
			Relx_strainhist[i] = (real *)realloc((void *)Relx_strainhist[i], step_relx*sizeof(real));
			Relx_stresshist[i] = (real *)realloc((void *)Relx_stresshist[i], step_relx*sizeof(real));
			Relx_strainhist[i][step_relx-1] = EpsAvg[i];
			Relx_stresshist[i][step_relx-1] = Stress[i];
		}

		/* calculate elastic energy */
		ElasticEnergy();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&elEnergy, &elEnergyAll, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		Relx_elEnergyhist = (real *)realloc((void *)Relx_elEnergyhist, step_relx*sizeof(real));
		Relx_elEnergyhist[step_relx-1] = elEnergyAll/NXY;
	
//		/* check stress-self balance */
//		Relx_Divehist = (real *)realloc((void *)Relx_Divehist, step_relx*sizeof(real));
//		Relx_Divehist[step_relx-1] = StressSelfBalance();

		/* clock advances */
		time_relx_hist = (real *)realloc((void *)time_relx_hist, step_relx*sizeof(real));
		time_relx_hist[step_relx-1] = time_relx;

		if(time_relx>TimeLimit){
			if(rank==0)
				printf("Pre-relaxing for %.3es time is completed!(Expecting %.3es)\n", time_relx, TimeLimit);
			break;
		}
		else {
			if(rank==0){
				printf("Relaxation step = %d (time = %.3es)\tElastic Energy = %.3e\nastep=%d\tkstep=%d\n\n", step_relx, time_relx, Relx_elEnergyhist[step_relx-1], tmp_astep, tmp_kstep);
				/* record elastic energy */
				sprintf(fname, "RelaxationElasticEnergy.out");
				fp = fopen(fname, "w");
				for(i=0; i<step_relx; i++){
					fprintf(fp, "%d\t%.6e\t%.6e\n", i, time_relx_hist[i], Relx_elEnergyhist[i]);
				}
				fclose(fp);

				
				sprintf(fname, "RelaxationStress.out");
				fp = fopen(fname, "w");
				for(i=0; i<step_relx; i++){
					fprintf(fp, "%d\t%.6e\t%.6e\t%.6e\t%.6e\n", i, time_relx_hist[i], Relx_stresshist[0][i],
							Relx_stresshist[1][i], Relx_stresshist[2][i]);
				}
				fclose(fp);

				sprintf(fname, "RelaxationEpsAvg.out");
				fp = fopen(fname, "w");
				for(i=0; i<step_relx; i++){
					fprintf(fp, "%d\t%.6e\t%.6e\t%.6e\t%.6e\n", i, time_relx_hist[i], Relx_strainhist[0][i],
							Relx_strainhist[1][i], Relx_strainhist[2][i]);
				}
				fclose(fp);
			}
		}
	}
	free(Relx_elEnergyhist);
	free(time_relx_hist);
	for(i=0; i<3; i++){
		free(Relx_strainhist[i]);
		free(Relx_stresshist[i]);
	}
	
	return;
}/* end PreRelxation() */

void ProjectResidualStress(char *fname, int nulty)
{





	return;
}/*end ProjectResidualStress()*/
