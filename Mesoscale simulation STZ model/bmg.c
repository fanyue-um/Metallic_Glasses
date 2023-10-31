#define BMG_GLOBALS_ONCE
#include "bmg.h"
#include "my_fft.h"
#undef BMG_GLOBALS_ONCE

FILE *f1;


int main(int argc, char **argv){
	int m, i, j, idx, de, m1, ii, jj, idx1, yd = 0;
	int N_print = 1;
	char fname[20];
	real rectimeLast = 0.;
	real elEnergyAll;
	int na;
	na=0;
//	printf("test0");

	strainhist = (real *)realloc((void *)strainhist, N_print*sizeof(real));
	stresshist = (real *)realloc((void *)stresshist, N_print*sizeof(real));
	elEnergyhist = (real *)realloc((void *)elEnergyhist, N_print*sizeof(real));
	strainhist[0] = stresshist[0] = 0.0;
	/* SOC analysis counters
	   countSOC[0] - avalanche size;
	   countSOC[1] - # of corresponding event; */
	sizeSOC = 1;
	countSOC[0] = (int *)realloc((void *)countSOC[0], sizeSOC*sizeof(int)); countSOC[0][sizeSOC-1] = 1;
	countSOC[1] = (int *)realloc((void *)countSOC[1], sizeSOC*sizeof(int)); countSOC[1][sizeSOC-1] = 0;

	/* Input file reading */
	f1 = fopen("input", "r");
	ReadInput(f1);
	printf("Input file read\n");
	fclose(f1);

	
	/* MPI initialized */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank==0) PrintFlags();
	/* Initialize global variables */
	Init();
	if(rank==0) printf("global variables initialized\n");
				printf("iE=%.3f\nimu=%.3f\ninu=%.3f\nilambda=%.3f\n", E, mu, nu, lambda);

  /* global file to record actual Q */
  if(rank==0) fq = fopen("actQ.out","w");

	/* Initialize MPI FFT */
	FFT_MPI_Init();
	if(rank==0) printf("MPI_FFT initialized\n");
	//printf("PE#%d: %d-%d\n", rank, lxs, lxs+lnx);

	/* Initialize random number generator */
	Rand_Init();
	if(rank==0) printf("Random generator initialized\n");

	/* alloc and initialize Eps0, DEps0, Sigma, etc */
	Alloc();

	/* Initialize k-space */
	Grad();
	if(rank==0) printf("k-space initialized\n");


	/* Fluctuation in stress-free energy barrier */
	Energy_Init();
	if(rank==0) printf("Energy fluctuation initialized\n");


	/* Generate SFTSs */
	DStrain();

	/* Simulation starts */
	if(rank==0){
		printf("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		printf("Simulation begins\n");
		printf("--------------------------------------\n\n");
		ParameterPrint();
	}


/* Gnerate residual stersses */
#if defined(RESIDUAL_STRESS_CONST)||defined(RESIDUAL_STRESS_RELX)
	/* generate distributed "pre-transformation" strain */
  int nullity=0;
	RdStrain();
#ifdef INHOM_ELAST
	//HomogeneousStressSolver(lambda, mu, Eps0, EpsAvg,
	//		XR, TotEps, Sigma, Stress, 0);
	InhomogeneousStressSolver(Lambda, Mu,
      Eps0, EpsAvg,
      lambda, mu,
			XR, TotEps, Sigma, Stress, 0);
#else
	HomStressSolver();
#endif
	ProjectResidualStress(ResStressDataFile, nullity);
	printf("pyz:%s\n",ResStressDataFile);
	exit(3);


	for(de=0; de<DE; de++){
		space_loop{
			idx = i*NY + j;
			Eps0[de][idx] = 0.0;
			Sigma[de][idx] = 0.0;
		}
	}

	PrintResidualStress();
#endif


	EpsAvg[0] = 0.;	EpsAvg[1] = 0.;	EpsAvg[2] = 0.;
#ifdef INHOM_ELAST
	//HomogeneousStressSolver(lambda, mu, Eps0, EpsAvg,
	//		XR, TotEps, Sigma, Stress, 0);
	InhomogeneousStressSolver(Lambda, Mu,
      Eps0, EpsAvg,
      lambda, mu,
			XR, TotEps, Sigma, Stress, 0);
#else
	HomStressSolver();
#endif

	if(rank == 0){
		printf("PE#%d:step=%d\ttime=%lf\tstrain=[%.3e, %.3e, %.3e]\n"
				"stress=[%.3e, %.3e, %.3e]\n\n",
				rank, step, totime, EpsAvg[0], EpsAvg[1], EpsAvg[2],
				Stress[0], Stress[1], Stress[2]);
	}
				printf("E=%.3f\nmu=%.3f\nnu=%.3f\nlambda=%.3f\n", E, mu, nu, lambda);

	MPI_Barrier(MPI_COMM_WORLD);
	while(1){
		step++;
		/* Softening factor */
		if(softening_FLAG==1) SofteningFactor();

		/* Activation Energy */
		ActEnergy(step);
		
		/* Sorting activation energy */
		flag = 0;	aflag = 0;
		for(m=0; m<MODES; m++)
			space_loop{
				idx = i*NY + j;
				//if(Q[m][idx] < 0 || fabs(Q[m][idx])<0.33) flag++;
				if(Q[m][idx] < 0) flag++;
			}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&flag, &aflag, 1, MPI_INT,
				MPI_SUM, MPI_COMM_WORLD);

		
		/* evolution */
		if(aflag >0){
			/* athermal plasticity */
			MPI_Barrier(MPI_COMM_WORLD);
			athermal_plasticity();
		}
		else{
			UpdateSOC();
			/* thermal OR pure elasticity */
			MPI_Barrier(MPI_COMM_WORLD);
			KMC();
		}
		
		/* calculate elastic energy */
		ElasticEnergy();
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allreduce(&elEnergy, &elEnergyAll, 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);

		/* print & record */
#if defined(TENSION)
		if(fabs(EpsAvg[0] - strainhist[N_print-1]) >=
	   			(printperiod * strainincrement * 0.99)){
#elif defined(SHEAR)
		if(fabs(EpsAvg[2] - strainhist[N_print-1]) >=
	   			(printperiod * strainincrement * 0.99)){
#endif
			N_print++;
			strainhist = (real *)realloc((void *)strainhist, N_print*sizeof(real));
			stresshist = (real *)realloc((void *)stresshist, N_print*sizeof(real));
			elEnergyhist = (real *)realloc((void *)elEnergyhist, N_print*sizeof(real));
	
#if defined(TENSION)
			strainhist[N_print-1] = EpsAvg[0];
			stresshist[N_print-1] = Stress[0];
#elif defined(COMPRESSION)
			strainhist[N_print-1] = EpsAvg[0];
			stresshist[N_print-1] = Stress[0];
#elif defined(SHEAR)
			strainhist[N_print-1] = EpsAvg[2];
			stresshist[N_print-1] = Stress[2];
#endif
/*avg Eps0*/
	real avg[3] = {0.0};
	space_loop{
		idx = i*NY + j;
		avg[0] += Eps0[0][idx];
		avg[1] += Eps0[1][idx];
		avg[2] += Eps0[2][idx];
	}
	avg[0] /= NXY;
	avg[1] /= NXY;
	avg[2] /= NXY;

			RecordPrep();
 	for(m1=0; m1<MODES; m1++){
		for(ii=0;ii<lnx;ii++)
			for(jj=0;jj<NY;jj++){
			idx1 = ii*NY + jj;

		  	Q0[m1][idx1]=Q1[m1][idx1];
			}
	}

			elEnergyhist[N_print-1] = elEnergyAll;
			if(rank == 0){
				printf("PE#%d:step=%d\ttime=%lf\tstrain=[%.3e, %.3e, %.3e]\n"
					"stress=[%.3e, %.3e, %.3e]\n",
					rank, step, totime, EpsAvg[0], EpsAvg[1], EpsAvg[2],
					Stress[0], Stress[1], Stress[2]);
				printf("step=%d\nastep=%d\nkstep=%d\nestep=%d\nN_print=%d\n\n", step, astep, kstep, estep, N_print);
			}

			/* strain & stress record */
		//	RecordPrep();

			/* Von Mises Strain */
			sprintf(fname, "Strain%03d_P%02d.out", N_print-1, rank);
			f1 = fopen(fname, "w");
			if(f1 == NULL){
				perror("Error opening file\n");
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			for(i=0; i<lnx; i++){
				for(j=0; j<NY; j++){
					idx = i*NY + j;
					fprintf(f1, "%lf ", EqStrain[idx]);
				}
				fprintf(f1, "\n");
			}
			fclose(f1);
			
			/* Von Mises Stress */
			sprintf(fname, "Stress%03d_P%02d.out", N_print-1, rank);
			f1 = fopen(fname, "w");
			if(f1 == NULL){
				perror("Error opening file\n");
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			for(i=0; i<lnx; i++){
				for(j=0; j<NY; j++){
					idx = i*NY + j;
					fprintf(f1, "%lf ", EqStress[idx]);
				}
				fprintf(f1, "\n");
			}
			fclose(f1);
			
#ifdef FULL_TENSOR_RECORD
			/* Strain and stress */
      for(i=0;i<3;i++){
			  sprintf(fname, "Strain_%02d", i+1);
        WriteScalarMPI(Eps0[i],fname,N_print-1);
			  sprintf(fname, "Stress_%02d", i+1);
        WriteScalarMPI(Sigma[i],fname,N_print-1);
      }
#endif


#ifdef HEAT_CONDUCT

#endif

			
				
			RecordFree();

			/* displacement record */
#ifdef INHOM_ELAST
	//HomogeneousStressSolver(lambda, mu, Eps0, EpsAvg,
	//		XR, TotEps, Sigma, Stress, 2);
	InhomogeneousStressSolver(Lambda, Mu,
      Eps0, EpsAvg,
      lambda, mu,
			XR, TotEps, Sigma, Stress, 2);
#else
			DispField(XR);
#endif

			/* X1 */
			sprintf(fname, "X1R%03d_P%02d.out", N_print-1, rank);
			f1 = fopen(fname, "w");
			if(f1 == NULL){
				perror("Error opening file\n");
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			for(i=0; i<lnx; i++){
				for(j=0; j<NY; j++){
					idx = i*NY + j;
					fprintf(f1, "%lf ", XR[0][idx]);
				}
				fprintf(f1, "\n");
			}
			fclose(f1);

			/* X2 */
			sprintf(fname, "X2R%03d_P%02d.out", N_print-1, rank);
			f1 = fopen(fname, "w");
			if(f1 == NULL){
				perror("Error opening file\n");
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			for(i=0; i<lnx; i++){
				for(j=0; j<NY; j++){
					idx = i*NY + j;
					fprintf(f1, "%lf ", XR[1][idx]);
				}
				fprintf(f1, "\n");
			}
			fclose(f1);



			/* stress-strain curve */
			if(rank == 0){
				f1 = fopen("SS.out", "w");
				for(i=0; i<N_print; i++){
					fprintf(f1, "%.6e\t%.6e\n", strainhist[i], stresshist[i]);
				}
				fclose(f1);
		
				f1 = fopen("SOC.out", "w");
				for(i=0; i<sizeSOC; i++){
					fprintf(f1, "%d\t%d\n", countSOC[0][i], countSOC[1][i]);
				}
				fclose(f1);
		
				f1 = fopen("ElastSOC.out", "w");
				for(i=0; i<sizeElaSOC; i++){
					fprintf(f1, "%.6e\t%.6e\n", elastSOC[0][i], elastSOC[1][i]);
				}
				fclose(f1);
		
				f1= fopen("ElasticEnergy.out", "w");
				for(i=0; i<N_print; i++){
					fprintf(f1, "%.6e\t%.6e\n", strainhist[i], elEnergyhist[i]);
				}
				fclose(f1);
					
				f1 = fopen("EA.out","w");
				fprintf(f1, "%.6e\t%.6e\t%.6e\n", EpsAvg[0], EpsAvg[1], EpsAvg[2]);
				fclose(f1);
		
				printf("step=%d\nastep=%d\nkstep=%d\nestep=%d\n", step, astep, kstep, estep);
			}
			/* EpsAvg */
			if(rank==0){
				f1 = fopen("EpsAvg.out","a");
				fprintf(f1, "%lf %lf %lf\n", EpsAvg[0], EpsAvg[1], EpsAvg[2]);
				fclose(f1);
			}
			/* avg */
			if(rank==0){
				f1 = fopen("avg.out","a");
				fprintf(f1, "%lf %lf %lf\n", avg[0], avg[1], avg[2]);
				fclose(f1);
			}


		}

#if defined(TENSION)
		if(fabs(EpsAvg[0]) > EndStrain){
			if(rank == 0) printf("PE#%d: macro strain exceeds %lf\n", rank, EndStrain);
			break;
		}

		if(Stress[0] < -0.5){
			if(rank == 0){
				printf("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				printf("Unstable system! Solution is divergent (occilating)\n");
				printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			}
			break;
		}
#elif defined(SHEAR)
		if(fabs(EpsAvg[2]) > EndStrain){
			if(rank == 0) printf("PE#%d: macro strain exceeds %lf\n", rank, EndStrain);
			break;
		}

		if(fabs(Stress[2]) < EPSILON){
			if(rank == 0){
				printf("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				printf("Unstable system! Solution is divergent (occilating)\n");
				printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			}
			break;
		}
#endif
	}

	/* data writing */

	/* stress-strain curve */
	if(rank == 0){
		f1 = fopen("SS.out", "w");
		for(i=0; i<N_print; i++){
			fprintf(f1, "%.6e\t%.6e\n", strainhist[i], stresshist[i]);
		}
		fclose(f1);

		f1 = fopen("SOC.out", "w");
		for(i=0; i<sizeSOC; i++){
			fprintf(f1, "%d\t%d\n", countSOC[0][i], countSOC[1][i]);
		}
		fclose(f1);

		f1 = fopen("ElastSOC.out", "w");
		for(i=0; i<sizeElaSOC; i++){
			fprintf(f1, "%.6e\t%.6e\n", elastSOC[0][i], elastSOC[1][i]);
		}
		fclose(f1);

		f1= fopen("ElasticEnergy.out", "w");
		for(i=0; i<N_print; i++){
			fprintf(f1, "%.6e\t%.6e\n", strainhist[i], elEnergyhist[i]);
		}
		fclose(f1);
			
		f1 = fopen("EA.out","w");
		fprintf(f1, "%.6e\t%.6e\t%.6e\n", EpsAvg[0], EpsAvg[1], EpsAvg[2]);
		fclose(f1);

		printf("step=%d\nastep=%d\nkstep=%d\nestep=%d\n", step, astep, kstep, estep);
	}
//	/* plastic strain */
//	sprintf(fname, "E1_P%d.out", rank);
//	f1=fopen(fname, "w");
//	space_loop{
//		idx = i*NY + j;
//		fprintf(f1, "%.6e", Eps0[0][idx]);
//		if(j==(NY-1))
//			fprintf(f1, "\n");
//		else
//			fprintf(f1, "\t");
//	}
//	fclose(f1);
//	
//	sprintf(fname, "E3_P%d.out", rank);
//	f1=fopen(fname, "w");
//	space_loop{
//		idx = i*NY + j;
//		fprintf(f1, "%.6e", Eps0[2][idx]);
//		if(j==(NY-1))
//			fprintf(f1, "\n");
//		else
//			fprintf(f1, "\t");
//	}
//	fclose(f1);


	MPI_Barrier(MPI_COMM_WORLD);
	printf("\nWrap up PE#%d\n", rank);

  if(rank==0) fclose(fq);

	if(rank==0){
		printf("\n\n--------------------------------------\n");
		printf("\nProgram finished\n");
		printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	}

	/* Free memory */
	Free();

	/* MPI finished */
	MPI_Finalize();

	return 0;
}
