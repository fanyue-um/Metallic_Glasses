#include "bmg.h"
#include "my_fft.h"

void ReadInput(FILE *in){
	int matchek = 0;

	matchek =  fscanf(in, "Moduli(E, nu):\n"
					"%lf %lf\n\n",
					&E, &nu);
	if(matchek != 2){
		perror("ReadInput: Could not match Moduli\n");
	}

	matchek = fscanf(in, "Stress-free energy barrier:\n"
					"%lf\n\n", &Q0amp);
	if(matchek !=1){
		perror("ReadInput: Could not match Stress-free energy barrier\n");
	}

	matchek = fscanf(in, "Applied strain rate:\n"
					"%lf\n\n", &strainrate);
	if(matchek !=1){
		perror("ReadInput: Could not match Applied strain rate\n");
	}

	matchek = fscanf(in, "Strain increment:\n"
					"%lf\n\n", &strainincrement);
	if(matchek != 1){
		perror("ReadInput: Could not match Strain increment\n");
	}

	matchek = fscanf(in, "Random seed:\n%d\n\n", &seed);
	if(matchek != 1){
		perror("ReadInput: Could not match Random seed\n");
	}

	matchek = fscanf(in, "Random ratio of fluctuation in energy barrier:\n"
					"%lf\n\n", &RandomRatio);
	if(matchek != 1){
		perror("ReadInput: Could not match Random ratio of fluctuation in energy barrier\n");
	}

	matchek = fscanf(in, "Print period:\n%d\n\n", &printperiod);
	if(matchek != 1){
		perror("ReadInput: Could not match Print period\n");
	}

	matchek = fscanf(in, "End Strain:\n%lf\n\n", &EndStrain);
	if(matchek != 1){
		perror("ReadInput: Could not match EndStrain\n");
	}

	matchek = fscanf(in, "Softening Flag:\n%d\n\n", &softening_FLAG);
	if(matchek != 1){
		perror("ReadInput: Could not match Softening Flag\n");
	}

	matchek = fscanf(in, "Relaxation Time:\n%lf\n\n", &RelxTime);
	if(matchek != 1){
		perror("ReadInput: Could not match Relaxation Time\n");
	}

	matchek = fscanf(in, "Characteristic Xtal shear:\n%lf\n\n", &XtalShear);
	if(matchek != 1){
		perror("ReadInput: Could not match Characteristic xtal shear\n");
	}

	matchek = fscanf(in, "Barrier for Xtal shear:\n%lf\n\n", &XtalQ0amp);
	if(matchek != 1){
		perror("ReadInput: Could not match Barrier for xtal shear\n");
	}

	matchek = fscanf(in, "Ratio for generating residual stress:\n%lf\n\n", &ResiRatio);
	if(matchek != 1){
		perror("ReadInput: Could not match Residual stress Ratio\n");
	}

	matchek = fscanf(in, "Residual stress data file:\n%s\n\n", &ResStressDataFile);
	if(matchek != 1){
		perror("ReadInput: Could not match Residual stress data file\n");
	}

  matchek = fscanf(in, "Volume fraction of crystal-like MRO:\n%lf\n\n", &xMROFrac);
  if(matchek !=1){
    perror("ReadInput: Could not match volume fraction of crystal-like MRO\n");
  }

  matchek = fscanf(in, "Initialization condition for MRO:\n%d\n\n", &InitMRO);
  if(matchek !=1){
    perror("ReadInput: Could not match initialization condition for MRO\n");
  }

	matchek = fscanf(in, "Stress-free energy barrier for crystal-like MRO:\n"
					"%lf\n\n", &xQ0amp);
	if(matchek !=1){
		perror("ReadInput: Could not match Stress-free energy barrier for crystal-like MRO\n");
	}

	return;

}/* end ReadInput() */

void ParameterPrint(void){
	printf("\nSimulation Geometry:\n");
	printf("System size:\t[%d x %d]\n", NX, NY);
	printf("Voxel size:\t[%.3fnm x %.3fnm]\n", LX, LY);

	printf("\nMaterial Properties at T=%f\n", T__K);
	printf("Activation energy barrier amplituide: \t%.3f\n", Q0amp);
	printf("Young's modulus:\t%.3f\n", E);
	printf("Shear modulus:\t%.3f\n", mu);
	printf("Poisson's ratio:\t%.3f\n", nu);
	
	printf("\nPEL Parameters:\n");
	printf("# of STZ modes:\t%d\n", MODES);
	printf("Stress-free energy barrier prefactor:\t%.3f\n", Q0amp);

	printf("\nSoftening Parameters:\n");
	printf("Permanent softening prefactor:\t%.3f\n", PermSoft);
	printf("Temporary softening prefactor:\t%.3f\n", TempSoft);
	printf("PermSoftening cap:\t%.3f\n", PermSofteningCap);

	printf("\nSimulation Environment:\n");
	printf("Applied strain rate:\t%.3e\n", strainrate);
	printf("Strain increment in recording data:\t%.3e\n", strainincrement);
	printf("End Strain:\t%.3f\n", EndStrain);
	if(softening_FLAG==0)
		printf("Softening model: OFF\n");
	else
		printf("Softening model: ON\n");
	printf("Relaxation time:\t%.3e\n", RelxTime);
	printf("Ratio for generating residual stress:\t%.3f\n", ResiRatio);

	return;
}/* end ParameterPrint() */

void DispField(real *XR[2]){
	int i, j, idx, de;
	fftw_complex F[2];
	fftw_complex G;
	fftw_complex *tmp[DE] = {0};

	for(de=0; de<DE; de++){
		tmp[de] = (fftw_complex *)fftw_malloc(lsize * sizeof(fftw_complex));
		if(!(tmp[de])) perror("failed to allocate tmp\n");
	}

	/***********************
	* supercell position *
	 ***********************/
	space_loop{
		idx = i*NY + j;
		XR[0][idx] = (i+lxs)*(1+EpsAvg[0])*LX + j*EpsAvg[2]*LY;
		XR[1][idx] = (i+lxs)*EpsAvg[2]*LX + j*(1+EpsAvg[1])*LY;
	}

	/***********************
	 * local displacement *
	 ***********************/
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

//		tmp[0][idx].re = (F[0].re * g[0][idx] - alpha * G.re * gg[0][idx]) / mu;
//		tmp[0][idx].im = (F[0].im * g[0][idx] - alpha * G.im * gg[0][idx]) / mu;
//		tmp[1][idx].re = (F[1].re * g[1][idx] - alpha * G.re * gg[1][idx]) / mu;
//		tmp[1][idx].im = (F[1].im * g[1][idx] - alpha * G.im * gg[1][idx]) / mu;
//		tmp[2][idx].re = (F[0].re * g[1][idx] + F[1].re * g[0][idx] - 2*alpha*G.re*gg[2][idx])/2/mu;
//		tmp[2][idx].im = (F[0].im * g[1][idx] + F[1].im * g[0][idx] - 2*alpha*G.im*gg[2][idx])/2/mu;

		tmp[0][idx].re = (F[0].im - alpha*G.im*g[0][idx])/mu/g2[idx];
		tmp[0][idx].im = (-F[0].re + alpha*G.re*g[0][idx])/mu/g2[idx];
		tmp[1][idx].re = (F[1].im - alpha*G.im*g[1][idx])/mu/g2[idx];
		tmp[1][idx].im = (-F[1].re + alpha*G.re*g[1][idx])/mu/g2[idx];
	}

	/* bwd FFT for total strain */
	for(de=0; de<2; de++){
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
			/***********************
			* final position *
			***********************/
			XR[de][idx] += fft_data[idx].re/NXY;
		}
	}
	for(de=0; de<DE; de++){
		fftw_free(tmp[de]);
	}

	return;
}/* end DispField() */

void DispPrep(void){
	int i, j;
	int idx;
	for(i=0; i<2; i++)
		XR[i] = (real*)malloc(lsize * sizeof(real));

	space_loop{
		idx = i*NY + j;
		XR[0][idx] = 0;
		XR[1][idx] = 0;
	}

	return;
}/* end DispPrep() */

void RecordPrep(void){
	int i, j;
	int idx;
	EqStrain = (real*)malloc(lsize * sizeof(real));
	EqStress = (real*)malloc(lsize * sizeof(real));

	space_loop{
		idx = i*NY + j;
		EqStrain[idx] = sqrt(Eps0[0][idx]*Eps0[0][idx] + Eps0[2][idx]*Eps0[2][idx]);
		EqStress[idx] = sqrt(Sigma[0][idx]*Sigma[0][idx] + Sigma[2][idx]*Sigma[2][idx]);
	}

	return;
}/* end RecordPrep() */

void DispFree(void){
	free(XR[0]);
	free(XR[1]);

	return;
}/* end DispFree() */

void RecordFree(void){
	free(EqStrain);
	free(EqStress);
	
	return;
}/* end RecordFree() */

void WriteScalarMPI(real *field, char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;

	sprintf(fname, "%s_S%03d.iout", s ,step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(real), MPI_real,
			MPI_real, "native", MPI_INFO_NULL);
	MPI_File_write(fp, field, lsize, MPI_real, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteScalarMPI()*/

void WriteScalarMPI_INT(int *field, char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;

	sprintf(fname, "%s_S%03d.iout", s ,step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, field, lsize, MPI_INT, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteScalarMPI_INT()*/

void WriteMROMPI_INT(int *field, char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
  int i,j,m;
  int *tmp;

	tmp = (int *)malloc(lnx*NY * sizeof(int));

  space_loop{
	     m = i*NY+j;
       //tmp[m] = field[m]>0;
       tmp[m] = field[m];
  }

	sprintf(fname, "%s_S%03d.iout", s ,step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmp, lsize, MPI_INT, &status);
	MPI_File_close(&fp);

  free(tmp);

	return;
}/*end WriteMROMPI_INT()*/

void PrintResidualStress(void){
	char fname[20];
	FILE *f1;
	int i, j, idx, de;

	for(de=0; de<3; de++){
		sprintf(fname, "RdStress%02d_P%02d.out", de+1, rank);
		f1=fopen(fname, "w");
		if(f1 == NULL){
			perror("Error opening file\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		for(i=0; i<lnx; i++){
			for(j=0; j<NY; j++){
				idx = i*NY + j;
				fprintf(f1, "%lf ", RdStress[de][idx]);
			}
			fprintf(f1, "\n");
		}
		fclose(f1);
	}

	return;
}/* end PrintResidualStress() */

void PrintFlags(void){

	printf("Active flags in bmgFlags.h:\n");
#ifdef HEAT_CONDUCT
	printf("#define HEAT_CONDUCT\n");
#endif

#ifdef XTAL_VARIANTS
	printf("#define XTAL_VARIANTS\n");
#endif

#ifdef COMPRESSION
	printf("#define COMPRESSION\n");
#endif

#ifdef TENSION
	printf("#define TENSION\n");
#endif

#ifdef SHEAR
	printf("#define SHEAR\n");
#endif

#ifdef PRE_RELAXATION
	printf("#define PRE_RELAXATION\n");
#endif

#ifdef RESIDUAL_STRESS
	printf("#define RESIDUAL_STRESS\n");
#endif

#ifdef UNIFORM_DEFECTS
	printf("#define UNIFORM_DEFECTS\n");
#endif

#ifdef LOW_T
	printf("#define LOW_T\n");
#endif

#ifdef HIGH_T_Vit1
	printf("#define HIGH_T_Vit1\n");
#endif

#ifdef RELX_BIMOLE_FLAG
	printf("#define RELX_BIMOLE_FLAG\n");
#endif

#ifdef RELX_FVT
	printf("#define RELX_FVT_FLAG\n");
#endif

#ifdef AVG_PERMSOFTENING
	printf("#define AVG_PERMSOFTENING\n");
#endif

#ifdef STRAIN_RATE_SENSITIVITY
	printf("#define STRAIN_RATE_SENSITIVITY\n");
#endif

	return;
}/* end PrintFlags() */

