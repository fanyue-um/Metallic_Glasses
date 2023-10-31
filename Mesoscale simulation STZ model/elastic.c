/* elastic.c */
#include "bmg.h"
#include "my_fft.h"
#include <fftw_mpi.h>

void DStrain(void){
	int i, j, m, idx;

#if defined(XTAL_VARIANTS)
	for(m=0; m<MODES; m++){
		space_loop{
			idx = i*NY + j;
			switch(m){
				case 0:
					DEps0[m][0][idx] = 0;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				case 1:
					DEps0[m][0][idx] = 0;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = -XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				case 2:
					DEps0[m][0][idx] = 0.5*sqrt(3)*XtalShear;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = -0.5*XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				case 3:
					DEps0[m][0][idx] = -0.5*sqrt(3)*XtalShear;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = 0.5*XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				case 4:
					DEps0[m][0][idx] = -0.5*sqrt(3)*XtalShear;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = -0.5*XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				case 5:
					DEps0[m][0][idx] = 0.5*sqrt(3)*XtalShear;
					DEps0[m][1][idx] = -DEps0[m][0][idx];
					DEps0[m][2][idx] = 0.5*XtalShear;
					Q0[m][idx] = XtalQ0amp;
					break;
				default:
					DEps0[m][0][idx] = 0;
					DEps0[m][1][idx] = 0;
					DEps0[m][2][idx] = 0;
					Q0[m][idx] = INFINITY;
					break;
			}
		}
	}
#else
	for(m=0; m<MODES; m++){
		space_loop{	
			idx = i*NY + j;
      if(MRO[idx]==1){  // six-fold
        if(m==0){
			    //DEps0[m][0][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    //DEps0[m][1][idx] = -DEps0[m][0][idx];
			    //DEps0[m][2][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[m][0][idx] = 0.5*GAMMA0 * cos(MROthe0[idx]);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = 0.5*GAMMA0 * sin(MROthe0[idx]);
        }else if(m<XMODES_SIX){
          real tmp_theta = m*2.*PI/XMODES_SIX*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[m][0][idx] = DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = -DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
        }else{
			    DEps0[m][0][idx] = 0.;
			    DEps0[m][1][idx] = 0.;
			    DEps0[m][2][idx] = 0.;
        }
      }
      else if(MRO[idx]==2){  // two-fold
        if(m==0){
			    //DEps0[m][0][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    //DEps0[m][1][idx] = -DEps0[m][0][idx];
			    //DEps0[m][2][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[m][0][idx] = 0.8*GAMMA0 * cos(MROthe0[idx]);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = 0.8*GAMMA0 * sin(MROthe0[idx]);
        }else if(m<XMODES_TWO){
          real tmp_theta = m*2.*PI/XMODES_TWO*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[m][0][idx] = DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = -DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
      }else{
			    DEps0[m][0][idx] = 0.;
			    DEps0[m][1][idx] = 0.;
			    DEps0[m][2][idx] = 0.;
        }
      }
      else if(MRO[idx]==3){  // ten-fold
        if(m==0){
			    //DEps0[m][0][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    //DEps0[m][1][idx] = -DEps0[m][0][idx];
			    //DEps0[m][2][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[m][0][idx] = 0.5*GAMMA0 * cos(MROthe0[idx]);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = 0.5*GAMMA0 * sin(MROthe0[idx]);
        }else if(m<XMODES_TEN){
          real tmp_theta = m*2.*PI/XMODES_TEN*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[m][0][idx] = DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[m][1][idx] = -DEps0[m][0][idx];
			    DEps0[m][2][idx] = -DEps0[0][0][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][idx]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
        }else{
			    DEps0[m][0][idx] = 0.;
			    DEps0[m][1][idx] = 0.;
			    DEps0[m][2][idx] = 0.;
        }
      }
      else if(MRO[idx]==0){ // glass
			  DEps0[m][0][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			  DEps0[m][1][idx] = -DEps0[m][0][idx];
			  DEps0[m][2][idx] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
      }
      else{ // invalid phase
	printf("%d\n",idx);
	printf("%d\n",MRO[idx]);
	      perror("Error 1 in MRO type (only {0,1,2,3} are allowed)!!\n");
		    MPI_Abort(MPI_COMM_WORLD, -1);
      }
		}
	}
#endif
	
	return;
}/* end DStrain() */

void RdStrain(void){
	int i, j, de, idx;
	real avg[3] = {0.0};
	real avgAll[3] = {0.};

	space_loop{
		idx = i*NY + j;
		Eps0[0][idx] = ResiRatio*0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
		Eps0[1][idx] = ResiRatio*0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
		Eps0[2][idx] = ResiRatio*0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
		avg[0] += Eps0[0][idx];
		avg[1] += Eps0[1][idx];
		avg[2] += Eps0[2][idx];
	}

	for(de=0; de<DE; de++){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&avg[de], &avgAll[de], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	}

	for(de=0; de<3; de++){
		avgAll[de] /= NXY;
		space_loop{
			idx = i*NY + j;
			Eps0[de][idx] -= avgAll[de];
		}
	}

	return;
}/* end RdStrain() */

void HomStressSolver(void){
	int i, j, de, idx;

	real *TrueStrain[DE] = {NULL};
	fftw_complex F[2];
	fftw_complex G;
	fftw_complex *tmp[DE] = {NULL};
	real StressAll[3] = {0.};
	fftw_complex K0strain[3]={0.};
	fftw_complex K0trace = {0.};
	
	for(de=0; de<DE; de++){
		tmp[de] = (fftw_complex *)fftw_malloc(lsize * sizeof(fftw_complex));
		if(!(tmp[de])) perror("failed to allocate tmp\n");
		TrueStrain[de] = (real *)malloc(lnx*NY * sizeof(real));
		if(!(TrueStrain[de])) perror("failed to allocate TrueStrain\n");
	}


	for(de=0; de<DE; de++){

		/* establish fft_data */
		space_loop{
			idx = i*NY + j;
			fft_data[idx].re = Eps0[de][idx];
			fft_data[idx].im = 0.0;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		/* forward FFT for strain */
		fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);

		/* store fft results into tmp */
		space_loop{
			idx = i*NY + j;
			if(isnan(fft_data[idx].re)||isnan(fft_data[idx].im)){
				printf("\n\nPE#%d-PYZ: NAN happens!\n\n",rank);
				MPI_Abort(MPI_COMM_WORLD,-4);
			}
			tmp[de][idx].re = fft_data[idx].re;
			tmp[de][idx].im = fft_data[idx].im;
		}
	}
#ifdef STRESS_FREE_BOUNDARY
	if(rank==0){
		for(de=0; de<DE; de++){
			K0strain[de].re = tmp[de][0].re;
			K0strain[de].im = tmp[de][0].im;
		}
	}
#endif

	/*strain-free stress */
	IsotropicMul(lambda, mu, tmp, tmp);

	/* F and G */
	space_loop{
		idx = i*NY + j;
		F[0].re = tmp[0][idx].re * g[0][idx] + tmp[2][idx].re * g[1][idx];
		F[0].im = tmp[0][idx].im * g[0][idx] + tmp[2][idx].im * g[1][idx];
		F[1].re = tmp[2][idx].re * g[0][idx] + tmp[1][idx].re * g[1][idx];
		F[1].im = tmp[2][idx].im * g[0][idx] + tmp[1][idx].im * g[1][idx];
		G.re = g[0][idx] * F[0].re + g[1][idx] * F[1].re;
		G.im = g[0][idx] * F[0].im + g[1][idx] * F[1].im;

		tmp[0][idx].re = (F[0].re * g[0][idx] - alpha * G.re * gg[0][idx]) / mu;
		tmp[0][idx].im = (F[0].im * g[0][idx] - alpha * G.im * gg[0][idx]) / mu;
		tmp[1][idx].re = (F[1].re * g[1][idx] - alpha * G.re * gg[1][idx]) / mu;
		tmp[1][idx].im = (F[1].im * g[1][idx] - alpha * G.im * gg[1][idx]) / mu;
		tmp[2][idx].re = (F[0].re * g[1][idx] + F[1].re * g[0][idx] - 2*alpha*G.re*gg[2][idx])/2/mu;
		tmp[2][idx].im = (F[0].im * g[1][idx] + F[1].im * g[0][idx] - 2*alpha*G.im*gg[2][idx])/2/mu;
	}

	/* The above calculated total strain in k-space need to satisfy 
	   the applied boundary condition, i.e. stress free or strain free.
	   for strain free, nothing special is needed here cause the total
	   strain in k-space at k=0 is zero according to above calculation.
	   If we want stree-free boundaries, the supercell has to change its
	   shape in order to relax the transformation-induced stress. */
#ifdef STRESS_FREE_BOUNDARY
	if(rank==0){
		for(i=0; i<2; i++){
			K0trace.re += tmp[i][0].re - K0strain[i].re;
			K0trace.im += tmp[i][0].im - K0strain[i].im;
		}
		K0trace.re = -lambda*K0trace.re/2.0/mu;
		K0trace.im = -lambda*K0trace.im/2.0/mu;
		/* The x-axis  is strain-controlled and the other
		   directions are stress-free */
	//	tmp[0][0].re = K0trace.re + K0strain[0].re;
	//	tmp[0][0].im = K0trace.im + K0strain[0].im;
		tmp[1][0].re = K0trace.re + K0strain[1].re;
		tmp[1][0].im = K0trace.im + K0strain[1].im;
		tmp[2][0].re = K0strain[2].re;
		tmp[2][0].im = K0strain[2].im;
	}
#endif

	/* backward FFT */
	for(de=0; de<DE; de++){
		/* establish fft_data */
		space_loop{
			idx = i*NY + j;
			fft_data[idx].re = tmp[de][idx].re;
			fft_data[idx].im = tmp[de][idx].im;
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		/* backward FFT for strain */
		fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);

		space_loop{
			idx = i*NY + j;
			if(isnan(fft_data[idx].re)||isnan(fft_data[idx].im)){
				printf("\n\nPE#%d-PYZ: NAN happens!\n\n",rank);
				MPI_Abort(MPI_COMM_WORLD,-4);
			}
			TrueStrain[de][idx] = fft_data[idx].re/(real)NXY + EpsAvg[de];
			TrueStrain[de][idx] -= Eps0[de][idx];
		}
	}

	/* compute stress field */
	RIsotropicMul(lambda, mu, TrueStrain);

	/* macro stress */
	for(de=0; de<DE; de++){
		Stress[de] = 0.0;
		space_loop{
			idx = i*NY + j;
			Stress[de] += Sigma[de][idx];
		}
	}

	for(de=0; de<DE; de++){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&Stress[de], &StressAll[de], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		Stress[de] = StressAll[de]/NXY;
	}
	
	for(de=0; de<DE; de++){
		fftw_free(tmp[de]);
		free(TrueStrain[de]);
	}

	return;
}/* end HomStressSolver() */

void IsotropicMul(real lambda, real mu, fftw_complex *strain[], fftw_complex *stress[]){
	int i, j, idx;
	fftw_complex tra;

	space_loop{
		idx = i*NY + j;
		tra.re = lambda * (strain[0][idx].re + strain[1][idx].re);
		tra.im = lambda * (strain[0][idx].im + strain[1][idx].im);
		stress[0][idx].re = tra.re + 2*mu*strain[0][idx].re;
		stress[0][idx].im = tra.im + 2*mu*strain[0][idx].im;
		stress[1][idx].re = tra.re + 2*mu*strain[1][idx].re;
		stress[1][idx].im = tra.im + 2*mu*strain[1][idx].im;
		stress[2][idx].re = 2*mu*strain[2][idx].re;
		stress[2][idx].im = 2*mu*strain[2][idx].im;
	}

	return;
}/* end IsotropicMul() */

void RIsotropicMul(real lambda, real mu, real *strain[]){
	int i, j, idx;
	real tra;

	space_loop{
		idx = i*NY + j;
		tra = lambda * (strain[0][idx] + strain[1][idx]);
		Sigma[0][idx] = tra + 2*mu*strain[0][idx];
		Sigma[1][idx] = tra + 2*mu*strain[1][idx];
		Sigma[2][idx] = 2 * mu * strain[2][idx];
	}

	return;
}/* end RIsotropicMul */

void VoxelShearTransform(int mm, int IDX){
	int i, j, m, idx1, ii, jj;
	real DeltaSoftening;
	real DeltaHeat;

#if defined(MICRO_MEMORY)
	/* Microscopic memory: system should remeber the transition
	   paths of at least the most recent generation */
	for(m=0; m<MODES; m++){
		for(i=0; i<DE; i++){
			PreDEps0[m][i][IDX] = DEps0[m][i][IDX];
			PreQ0[m][IDX] = Q0[m][IDX];
		}
	}
#endif

	for(i=0; i<DE; i++){
		Eps0[i][IDX] += DEps0[mm][i][IDX];
	}

	DeltaHeat = (Sigma[0][IDX] * DEps0[mm][0][IDX] + Sigma[1][IDX] * DEps0[mm][1][IDX] +
			2 * Sigma[2][IDX] * DEps0[mm][2][IDX]) * 1E9;	/* J/m^3 */
	Tlocal[IDX] += fabs(DeltaHeat) * thermal_coeff1;

	DeltaSoftening = DEps0[mm][0][IDX]*DEps0[mm][0][IDX] + 
			DEps0[mm][1][IDX]*DEps0[mm][1][IDX] + 2*DEps0[mm][2][IDX]*DEps0[mm][2][IDX];
	if(MRO[IDX]==1)
	PermSoftening[IDX] += DeltaSoftening*PermSoftmro;
	else
	PermSoftening[IDX] += DeltaSoftening*PermSoft;
	if(PermSoftening[IDX] > PermSofteningCap) PermSoftening[IDX] = PermSofteningCap;

	if(MRO[IDX]==1)
	TempSoftening[IDX] = DeltaSoftening*TempSoftmro;
	else
	TempSoftening[IDX] = DeltaSoftening*TempSoft;
	tlast[IDX] = totime;
#if defined(XTAL_VARIANTS)
	return;
//#elif defined(MICRO_MEMORY)
//	for(i=0; i<MODES-1; i++){
//		Q0[i][IDX] = Q0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
//		DEps0[i][0][IDX] = 0.5 * GAMMA0 * gsl_ran_gaussian(r, 1.0);
//		DEps0[i][1][IDX] = -DEps0[i][0][IDX];
//		DEps0[i][2][IDX] = 0.5 * GAMMA0 * gsl_ran_gaussian(r, 1.0);
//#ifdef STZ_DILATATION
//		real tmpDilate;
//		tmpDilate = DILATE0*gsl_ran_gaussian(r,1.0);
//		DEps0[i][0][IDX] += fabs(tmpDilate/2.0);
//		DEps0[i][1][IDX] += fabs(tmpDilate/2.0);
//#endif
//	}
//	/* The reverse path leads to previous state.
//	 Also notice that the barrier should be somehow
//	 altered, otherwise matrix-back-froce will always
//	 transform it back and forth repeatedly */
//	Q0[MODES-1][IDX] = Q0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
//	DEps0[MODES-1][0][IDX] = -PreDEps0[mm][0][IDX];
//	DEps0[MODES-1][1][IDX] = -PreDEps0[mm][1][IDX];
//	DEps0[MODES-1][2][IDX] = -PreDEps0[mm][2][IDX];
//
//	//Q0[MODES-1][IDX] = PreQ0[mm][IDX];
//	//DEps0[MODES-1][0][IDX] = PreDEps0[mm][0][IDX];
//	//DEps0[MODES-1][1][IDX] = PreDEps0[mm][1][IDX];
//	//DEps0[MODES-1][2][IDX] = PreDEps0[mm][2][IDX];
//#ifdef STZ_DILATATION
//	real tmpDilate;
//	tmpDilate = DILATE0*gsl_ran_gaussian(r,1.0);
//	DEps0[MODES-1][0][IDX] += fabs(tmpDilate/2.0);
//	DEps0[MODES-1][1][IDX] += fabs(tmpDilate/2.0);
//#endif

#else
		//		printf("EpsAvg0=%.6e\n", EpsAvg[0]);
		//		printf("Stress0=%.6e\n", Stress[0]);
		//		printf("stresshist=%.6e\n", stresshist[N_print-1]);
		//		printf("strainhist=%.6e\n", strainhist[N_print-1]);
/*avg Eps0*/
	real avg[3] = {0.0};
	space_loop{
		idx1 = i*NY + j;
		avg[0] += Eps0[0][idx1];
		avg[1] += Eps0[1][idx1];
		avg[2] += Eps0[2][idx1];
	}
	avg[0] /= NXY;
	avg[1] /= NXY;
	avg[2] /= NXY;

	for(i=0; i<MODES; i++){
    if(MRO[IDX]==1){  // six-fold
		    Q0[i][IDX]= xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
#if 0   // assume static
        if(i==0){
			    DEps0[i][0][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
        }else if(i<XMODES_SIX){
          real tmp_theta = i*2.*PI/XMODES_SIX*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[i][0][IDX] = DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = -DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
        }else{
          /* do nothing */
        }
#endif
    }
    else if(MRO[IDX]==2){  // two-fold
		    Q0[i][IDX]= xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
#if 0   // assume static
        if(i==0){
			    DEps0[i][0][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
        }else if(i<XMODES_TWO){
          real tmp_theta = i*2.*PI/XMODES_TWO*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[i][0][IDX] = DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = -DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
        }else{
          /* do nothing */
        }
#endif
    }
    else if(MRO[IDX]==3){  // ten-fold
		    Q0[i][IDX]= xQ0amp * (1-RandomRatio + RandomRatio*(gsl_rng_uniform(r)));
#if 0   // assume static
        if(i==0){
			    DEps0[i][0][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = 0.5*GAMMA0 * gsl_ran_gaussian(r, 1.0);
        }else if(i<XMODES_TEN){
          real tmp_theta = i*2.*PI/XMODES_TEN*(1.+0.2*(2*gsl_rng_uniform(r)-1.));
			    DEps0[i][0][IDX] = DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta);
			    DEps0[i][1][IDX] = -DEps0[i][0][IDX];
			    DEps0[i][2][IDX] = -DEps0[0][0][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*sin(2*tmp_theta)
            +DEps0[0][2][IDX]*(1.+0.2*(2*gsl_rng_uniform(r)-1.))*cos(2*tmp_theta);
        }else{
          /* do nothing */
        }
#endif
    }
    else if(MRO[IDX]==0){ // glass 
		  DEps0[i][0][IDX] = 0.5 * GAMMA0 * gsl_ran_gaussian(r, 1.0);
		  DEps0[i][1][IDX] = -DEps0[i][0][IDX];
		  DEps0[i][2][IDX] = 0.5 * GAMMA0 * gsl_ran_gaussian(r, 1.0);
    }
    else{
	    perror("Error: invalid MRO type in VoxelTransform()!!\n");
		  MPI_Abort(MPI_COMM_WORLD, -1);
    }
	}
	return;
#endif
}/* end VoxelShearTransform() */

void isotropicinv(real lambda, real mu, real *stress, real *strain){
	real tra;
	tra = -lambda/2/mu/(2*lambda + 2*mu) * (stress[0] + stress[1]);

	strain[0] = tra + stress[0]/2/mu;
	strain[1] = tra + stress[1]/2/mu;
	strain[2] = stress[2]/2/mu;

	return;
}/* end isotropicinv() */

void ElasticEnergy(void){
	real tmpElStrain[DE] = {0.};
	real tmpSigma[DE] = {0.};
	int i, j, ii, idx;

	elEnergy = 0.;
	space_loop{
		idx = i*NY + j;
		for(ii=0; ii<DE; ii++)
			tmpSigma[ii] = Sigma[ii][idx];
		isotropicinv(lambda, mu, tmpSigma, tmpElStrain);
		elEnergy += 0.5*(tmpSigma[0]*tmpElStrain[0] + tmpSigma[1]*tmpElStrain[1] + 2*tmpSigma[2]*tmpElStrain[2]);
	}
	return;
}

void HeatConducting(real detaTime){
	int i, j, idx;

	space_loop{
		idx = i*NY + j;
		fft_data[idx].re = Tlocal[idx];
		fft_data[idx].im = 0.0;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	/* forward FFT for strain */
	fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);

	space_loop{
		idx = i*NY + j;
		fft_data[idx].re *= exp(-g2[idx] * detaTime);
		fft_data[idx].im *= exp(-g2[idx] * detaTime);
		if(isnan(fft_data[idx].re)||isnan(fft_data[idx].im)){
			printf("\n\nPE#%d-PYZ: NAN happens!\n\n",rank);
			MPI_Abort(MPI_COMM_WORLD,-4);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	/* backward FFT for strain */
	fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
	
	space_loop{
		idx = i*NY + j;
		Tlocal[idx] = fft_data[idx].re/NXY;
	}

	return;
}/* end HeatConducting() */

void CheckStressSelfBalance(void){
	int i, j, de, idx;

	fftw_complex *SIGMA[DE] = {NULL};
	real residual[2]={0.0};
	real residualAll[2]={0.0};


	/* self-balancing the resulted stress */
	/* k-space stress */
	
	for(de=0; de<DE; de++){
		SIGMA[de] = (fftw_complex *)fftw_malloc(lsize * sizeof(fftw_complex));
		if(!(SIGMA[de])) perror("failed to allocate tmp\n");
	}


	for(de=0; de<DE; de++){

		/* establish fft_data */
		space_loop{
			idx = i*NY + j;
			fft_data[idx].re = Sigma[de][idx];
			fft_data[idx].im = 0.0;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		/* forward FFT to get k-space stress */
		fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
		space_loop{
			idx = i*NY + j;
			if(isnan(fft_data[idx].re)||isnan(fft_data[idx].im)){
				printf("\n\nPE#%d-PYZ: NAN happens!\n\n",rank);
				MPI_Abort(MPI_COMM_WORLD,-4);
			}
			SIGMA[de][idx].re = fft_data[idx].re;
			SIGMA[de][idx].im = fft_data[idx].im;
		}
	}

	residual[0] = 0.0; residual[1] = 0.0;
	space_loop{
		idx = i*NY + j;

		residual[0] += SIGMA[0][idx].re*g[0][idx]+SIGMA[2][idx].re*g[1][idx];
		residual[1] += SIGMA[2][idx].re*g[0][idx]+SIGMA[1][idx].re*g[1][idx];
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for(i=0; i<2; i++){
		MPI_Allreduce(&residual[i], &residualAll[i], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
	//	residualAll[i] /= (real)NXY;
	}

	if(rank == 0)
		printf("Stress-self balance check in k-space: Residual = [%e, %e]\n", residualAll[0], residualAll[1]);


	for(de=0; de<DE; de++){
		fftw_free(SIGMA[de]);
	}

	return;
}/* end CheckStressSelfBalance() */

#ifdef INHOM_ELAST
static void IsotropicMul_C(real lambda, real mu, fftw_complex **strain, fftw_complex **stress)
{
	fftw_complex tra;
  int i,j,idx;

	space_loop{
    idx = i*NY+j;
		tra.re = lambda * (strain[0][idx].re + strain[1][idx].re);
		tra.im = lambda * (strain[0][idx].im + strain[1][idx].im);
		stress[0][idx].re = tra.re + 2*mu*strain[0][idx].re;
		stress[0][idx].im = tra.im + 2*mu*strain[0][idx].im;
		stress[1][idx].re = tra.re + 2*mu*strain[1][idx].re;
		stress[1][idx].im = tra.im + 2*mu*strain[1][idx].im;
		stress[2][idx].re = 2*mu*strain[2][idx].re;
		stress[2][idx].im = 2*mu*strain[2][idx].im;
	}

	return;
}/* end IsotropicMul_C() */

static void HadamardProdnSum_C(fftw_complex *m1, real *m2,
		fftw_complex *n1, real *n2, fftw_complex *out)
{
  int i,j,idx;

	space_loop{
    idx = i*NY+j;
		out[idx].re = m1[idx].re*m2[idx] + n1[idx].re*n2[idx];
		out[idx].im = m1[idx].im*m2[idx] + n1[idx].im*n2[idx];
	}

	return;
}/*end HadamardProdnSum_C()*/

static void HadamardProd2nSSum_C(fftw_complex *m1, real *m2,
		fftw_complex *m3, real *m4, real s0,
		fftw_complex *n1, real *n2, real s, fftw_complex *out)
{
  int i,j,idx;

	space_loop{
    idx = i*NY+j;
		out[idx].re = s0*(m1[idx].re*m2[idx] + m3[idx].re*m4[idx]) +
			s*(n1[idx].re*n2[idx]);
		out[idx].im = s0*(m1[idx].im*m2[idx] + m3[idx].im*m4[idx]) +
			s*(n1[idx].im*n2[idx]);
	}

	return;
}/*end HadamardProd2nSSum_C()*/

static void IsotropicMul_R2(real lambda, real mu, real **totstraino,
		real **strain0o, real **stresso)
{
  int i,j,pIDX;
	real tra;
	space_loop{
    pIDX = i*NY+j;
		tra = lambda * ((totstraino[0][pIDX]-strain0o[0][pIDX])+
				(totstraino[1][pIDX]-strain0o[1][pIDX]));
		stresso[0][pIDX] = tra + 2.0*mu*(totstraino[0][pIDX]-strain0o[0][pIDX]);
		stresso[1][pIDX] = tra + 2.0*mu*(totstraino[1][pIDX]-strain0o[1][pIDX]);
		stresso[2][pIDX] = 2.0*mu*(totstraino[2][pIDX]-strain0o[2][pIDX]);
	}

	return;
}/*end IsotropicMul_R2()*/

static void HadamardProdnSSum_C(fftw_complex *m1, real *m2, real s0,
		fftw_complex *n1, real *n2, real s, fftw_complex *out)
{
  int i,j,idx;
	space_loop{
    idx = i*NY+j;
		out[idx].re = s0*(m1[idx].re*m2[idx]) + s*(n1[idx].re*n2[idx]);
		out[idx].im = s0*(m1[idx].im*m2[idx]) + s*(n1[idx].im*n2[idx]);
	}

	return;
}/*end HadamardProdnSSum_C()*/

static void AvgStress(real **so, real *stress)
{
	real sum, glb_sum;
	int de;
  int i,j,pIDX;
	for(de=0;de<3;de++){
		sum=0.0;
		glb_sum=0.0;
		space_loop{
      pIDX = i*NY+j;
			sum += so[de][pIDX];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&sum, &glb_sum, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		stress[de] = glb_sum/NXY;
	}

	return;
}/*end AvgStress()*/

static void HadamardSiDivd(fftw_complex *m1, real *m2, real s, fftw_complex *out)
{
  int i,j,idx;
	space_loop{
    idx = i*NY+j;
		out[idx].re = m1[idx].im*s/m2[idx];
		out[idx].im = -1.0*m1[idx].re*s/m2[idx];
	}

	return;
}/*end HadamardSiDivd()*/

static void HadamardSumnSProd(fftw_complex *m1, real s, fftw_complex *n1, real *n2,
		fftw_complex *out)
{
  int i,j,idx;
	space_loop{
    idx = i*NY+j;
		out[idx].re = m1[idx].re + s*(n1[idx].re*n2[idx]);
		out[idx].im = m1[idx].im + s*(n1[idx].im*n2[idx]);
	}

	return;
}/*end HadamardSumnSProd()*/

void HomogeneousStressSolver(real lambda, real mu, real **eps0o, real *epsavg,
		real **uo, real **epso, real **sigmao, real *stress, int uflag)
{
	/*****************
	  the last uflag indicates whether the displacement field will be calculated or not
	  *****************/
	int i,j,idx,de;
	real alpha = (lambda+mu)/(lambda+2.0*mu);
	for(de=0;de<DE;de++){
    space_loop{
      idx = i*NY+j;
      fft_data[idx].re = eps0o[de][idx];
      fft_data[idx].im = 0.;
    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
    space_loop{
      idx = i*NY+j;
      EPSnSIGMA[de][idx].re = fft_data[idx].re;
      EPSnSIGMA[de][idx].im = fft_data[idx].im;
    }
	}

	IsotropicMul_C(lambda, mu, EPSnSIGMA, EPSnSIGMA);

	HadamardProdnSum_C(EPSnSIGMA[0],g[0],EPSnSIGMA[2],g[1],Fiso[0]);
	HadamardProdnSum_C(EPSnSIGMA[2],g[0],EPSnSIGMA[1],g[1],Fiso[1]);
	HadamardProdnSum_C(Fiso[0],g[0],Fiso[1],g[1],Giso);
	HadamardProdnSSum_C(Fiso[0], g[0], 1.0/mu, Giso, gg[0], -1.0*alpha/mu, EPSnSIGMA[0]);
	HadamardProdnSSum_C(Fiso[1], g[1], 1.0/mu, Giso, gg[1], -1.0*alpha/mu, EPSnSIGMA[1]);
	HadamardProd2nSSum_C(Fiso[0], g[1], Fiso[1], g[0], 0.5/mu, Giso, gg[2], -1.0*alpha/mu, EPSnSIGMA[2]);

	if(uflag==0){	// only calculate strain
		for(de=0;de<DE;de++){
      space_loop{
        idx = i*NY+j;
        fft_data[idx].re = EPSnSIGMA[de][idx].re;
        fft_data[idx].im = EPSnSIGMA[de][idx].im;
      }
	    MPI_Barrier(MPI_COMM_WORLD);
	    fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
      space_loop{
        idx = i*NY+j;
        epso[de][idx] = fft_data[idx].re/NXY + epsavg[de];
      }
		}
		IsotropicMul_R2(lambda, mu, epso, eps0o, sigmao);
		AvgStress(sigmao,stress);

	}
	else if(uflag==1){	// calculate both strain and displacement
		for(de=0;de<3;de++){
      space_loop{
        idx = i*NY+j;
        fft_data[idx].re = EPSnSIGMA[de][idx].re;
        fft_data[idx].im = EPSnSIGMA[de][idx].im;
      }
	    MPI_Barrier(MPI_COMM_WORLD);
	    fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
      space_loop{
        idx = i*NY+j;
        epso[de][idx] = fft_data[idx].re/NXY + epsavg[de];
      }
		}
		IsotropicMul_R2(lambda, mu, epso, eps0o, sigmao);
		AvgStress(sigmao,stress);

		for(de=0;de<2;de++){
      // EPSnSIGMA[2] is used as a tmp variable
			HadamardSumnSProd(Fiso[de], -1.0*alpha, Giso, g[de], EPSnSIGMA[2]);
			HadamardSiDivd(EPSnSIGMA[2],g2,1.0/mu,EPSnSIGMA[de]);

      space_loop{
        idx = i*NY+j;
        fft_data[idx].re = EPSnSIGMA[de][idx].re;
        fft_data[idx].im = EPSnSIGMA[de][idx].im;
      }
	    MPI_Barrier(MPI_COMM_WORLD);
	    fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
      space_loop{
        idx = i*NY+j;
        uo[de][idx] = fft_data[idx].re/NXY;
      }
		}

		/***********************
		* supercell position *
		 ***********************/
		space_loop{
			idx = i*NY + j;
			uo[0][idx] += (i+lxs)*(1+epsavg[0])*LX + j*epsavg[2]*LY;
			uo[1][idx] += (i+lxs)*epsavg[2]*LX + j*(1+epsavg[1])*LY;
		}
	}
	else{	// only calculate displacement
		for(de=0;de<2;de++){
      // EPSnSIGMA[2] is used as a tmp variable
			HadamardSumnSProd(Fiso[de], -1.0*alpha, Giso, g[de], EPSnSIGMA[2]);
			HadamardSiDivd(EPSnSIGMA[2],g2,1.0/mu,EPSnSIGMA[de]);

      space_loop{
        idx = i*NY+j;
        fft_data[idx].re = EPSnSIGMA[de][idx].re;
        fft_data[idx].im = EPSnSIGMA[de][idx].im;
      }
	    MPI_Barrier(MPI_COMM_WORLD);
	    fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
      space_loop{
        idx = i*NY+j;
        uo[de][idx] = fft_data[idx].re/NXY;
      }
		}

		/***********************
		* supercell position *
		 ***********************/
		space_loop{
			idx = i*NY + j;
			uo[0][idx] += (i+lxs)*(1+epsavg[0])*LX + j*epsavg[2]*LY;
			uo[1][idx] += (i+lxs)*epsavg[2]*LX + j*(1+epsavg[1])*LY;
		}
	}

	return;
}/*end HomogeneousStressSolver()*/

static void IsotropicMul_RF(real *lambdao, real *muo, real **straino,
		real **stresso)
{
  int i,j,pIDX;
	real tra;

	space_loop{
    pIDX = i*NY+j;
		tra = lambdao[pIDX] * (straino[0][pIDX] + straino[1][pIDX]);
		stresso[0][pIDX] = tra + 2*muo[pIDX]*straino[0][pIDX];
		stresso[1][pIDX] = tra + 2*muo[pIDX]*straino[1][pIDX];
		stresso[2][pIDX] = 2*muo[pIDX]*straino[2][pIDX];
	}

	return;
}/*end IsotropicMul_RF()*/

static void IsotropicMul_RF_D(real *lambdao, real *muo,
		real lambda, real mu,
		real **straino, real **stresso)
{
	real tra;
  int i,j,pIDX;

	space_loop{
    pIDX = i*NY+j;
		tra = (lambda-lambdao[pIDX]) * (straino[0][pIDX] + straino[1][pIDX]);
		stresso[0][pIDX] = tra + 2*(mu-muo[pIDX])*straino[0][pIDX];
		stresso[1][pIDX] = tra + 2*(mu-muo[pIDX])*straino[1][pIDX];
		stresso[2][pIDX] = 2*(mu-muo[pIDX])*straino[2][pIDX];
	}

	return;
}/*end IsotropicMul_RF_D()*/

static void IsotropicInv(real lambda, real mu, real **stresso,
		real **straino)
{
  int i, j, pIDX;
	real tra;

	space_loop{
    pIDX = i*NY+j;
		tra = -1.0*lambda/2.0/mu/(2.0*lambda+2.0*mu)*(stresso[0][pIDX]+stresso[1][pIDX]);
		straino[0][pIDX] = tra + stresso[0][pIDX]/2.0/mu;
		straino[1][pIDX] = tra + stresso[1][pIDX]/2.0/mu;
		straino[2][pIDX] = stresso[2][pIDX]/2.0/mu;
	}

	return;
}/*end IsotropicInv()*/

static real AvgChange(real **strainlast, real **strain) 
{
	int i,j,pIDX,de;
	real avg = 0.0;
	real glb_avg = 0.0;
	for(de=0;de<3;de++){
		space_loop{
      pIDX = i*NY+j;
			avg += pow((strain[de][pIDX]-strainlast[de][pIDX]),2.0);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&avg, &glb_avg, 1, MPI_real,
			MPI_SUM, MPI_COMM_WORLD);
	glb_avg /= NXY;
	return(sqrt(glb_avg));

}/*end AvgChange()*/

void InhomogeneousStressSolver(real *Lambdao, real *Muo, 
		real **eps0o, real *epsavg,
		real lambda, real mu,
		real **uo, real **epso, real **sigmao, real *stress, int uflag)
{
	/*************************************************************************
	  Lambda-o and Mu-o define the elastically inhomogeneous system and
	  eps0-o is the associated transformation strain field.
	  lambda and mu define the reference homogeneous system.
	  The last row of arguments are the final solutions
	*************************************************************************/

	int i,j,pIDX,de;
	real convg, normepsavg;

	/* Need an initial guess for the total strain field */
	for(de=0;de<DE;de++){
		space_loop{
      pIDX = i*NY+j;
			epso[de][pIDX] = eps0o[de][pIDX]+epsavg[de];
		}
	}
	IsotropicMul_RF(Lambdao,Muo,eps0o,Tau0o);
	IsotropicMul_RF_D(Lambdao,Muo,lambda,mu,epso,Tau1o);
	for(de=0;de<DE;de++){
		space_loop{
      pIDX = i*NY+j;
			Tauo[de][pIDX] = Tau0o[de][pIDX] + Tau1o[de][pIDX];
		}
	}
	/* epso is used here to store the transformation
	   strain fields of reference system */
	IsotropicInv(lambda,mu,Tauo,epso);

	int step=0;
	normepsavg = sqrt(epsavg[0]*epsavg[0]+epsavg[1]*epsavg[1]+epsavg[2]*epsavg[2]);
	while(1){
    step++;

		/* record the reference system */
		for(de=0;de<3;de++){
			space_loop{
        pIDX = i*NY+j;
				epsolast[de][pIDX] = epso[de][pIDX];
			}
		}

		/* solve the reference homogeneous system */
		HomogeneousStressSolver(lambda,mu,epsolast,epsavg,
				uo, epso, sigmao, stress, 0); // epso stores the total strain
		IsotropicMul_RF_D(Lambdao,Muo,lambda,mu,epso,Tau1o);
		/* total stress
		   Note the first part is fixed by the inhom system itself */
		for(de=0;de<3;de++){
			space_loop{
        pIDX = i*NY+j;
				Tauo[de][pIDX] = Tau0o[de][pIDX]+Tau1o[de][pIDX];
			}
		}
		/* epso is used here to store the transformation
			strain fields of reference system */
		IsotropicInv(lambda,mu,Tauo,epso);

    real mix=0.5;
		for(de=0;de<DE;de++){
			space_loop{
        pIDX = i*NY+j;
				epso[de][pIDX] = mix*epso[de][pIDX]+(1-mix)*epsolast[de][pIDX];
			}
		}
    

		convg = AvgChange(epsolast,epso);
		
		if((convg <((1E-3)*normepsavg))||normepsavg<1E-4){
			if(uflag!=0){
				/* update the displacement field */
		    HomogeneousStressSolver(lambda,mu,epsolast,epsavg,
				    uo, epso, sigmao, stress, uflag); 
			}
			break;
		}

    if(step>100){
      break;
    }
	}

	return;
}/* end InhomogeneousStressSolver() */

#endif
