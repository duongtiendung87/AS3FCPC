#include "mylib.h"
int PERCENT_SELECT = 20;
int PERCENT_WRONG = 30;
int PERCENT_REMOVE = 30;

int Bool[100];
int permu[100];
int *clust = NULL;
double caval = 0;
double alpha_expert = 1.0;

int **must_link = NULL;
int **cannot_link = NULL;

// extern double data[N][D];
// extern int label[N];
// void outPutSFCM(char *filename, double **v, double **u, int n, int iter){
//     int i,j,k;
//     FILE *f;
//     f = fopen(filename,"w");
//     if(!f){
//          printf("Can't open file %s!",filename);
//          scanf("%d",&i);
//          exit(0);
//     }
//    fprintf(f,"N=%d\n",n);
//	fprintf(f,"R=%d\n",D);
//	fprintf(f,"C=%d\n",C);
//	fprintf(f,"U:\n");
//	for(i = 0;i<n;i++){
//		for(j = 0;j<C;j++){
//			fprintf(f,"%7.6lf ",u[i][j]);
//		}
//		fprintf(f,"\n");
////		printf("i = %d\n",i);
//	}
//	fprintf(f,"V:\n");
//    for(i = 0;i<C;i++){
//        for(j = 0;j<D;j++){
//            fprintf(f,"%7.6lf ",v[i][j]);
//        }
//        fprintf(f,"\n");
//    }
//    fprintf(f,"Leap= %d\n",iter);
//    fclose(f);
//}
void detectBoundaries()
{
	if (queried == NULL)
		queried = (int *)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)
	{
		double max1 = 0, max2 = 0;
		for (int j = 0; j < C; j++)
		{
			if (U[i][j] > max1)
			{
				max2 = max1;
				max1 = U[i][j];
			}
			else if (U[i][j] > max2)
			{
				max2 = U[i][j];
			}
		}
		double diff = max1 - max2;
		queried[i] = (diff < 0.2);
	}
}
void adjustMembership()
{
	for (int i = 0; i < N; i++)
	{
		if (queried[i])
		{
			double max1 = 0, max2 = 0;
			int idx1 = -1, idx2 = -1;
			for (int j = 0; j < C; j++)
			{
				if (U[i][j] > max1)
				{
					max2 = max1;
					idx2 = idx1;
					max1 = U[i][j];
					idx1 = j;
				}
				else if (U[i][j] > max2)
				{
					max2 = U[i][j];
					idx2 = j;
				}
			}
			double epsilon = 0.01;
			if ((max1 - max2) > epsilon)
			{
				double mid = (U[i][idx1] + U[i][idx2]) / 2.0;
				U[i][idx1] = mid - epsilon / 2.0;
				U[i][idx2] = mid + epsilon / 2.0;
			}
		}
	}
}
void applyConstraints()
{
	// Allocate must_link and cannot_link
	must_link = (int **)malloc(N * sizeof(int *));
	cannot_link = (int **)malloc(N * sizeof(int *));
	for (int i = 0; i < N; i++)
	{
		must_link[i] = (int *)calloc(N, sizeof(int));
		cannot_link[i] = (int *)calloc(N, sizeof(int));
	}
	for (int i = 0; i < N; i++)
	{
		if (queried[i])
		{
			label[i] = CLASS[i];
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (label[i] != -1 && label[j] != -1)
			{
				if (label[i] == label[j])
					must_link[i][j] = 1;
				else
					cannot_link[i][j] = 1;
			}
		}
	}
}

void xuat()
{
	double max, tmp;
	int i, j, k, count = 0;

	if (clust == NULL)
	{
		for (i = 0; i < N; i++)
		{
			max = 0;
			k = 0;
			for (j = 0; j < C; j++)
			{
				if (max < U[i][j])
				{
					max = U[i][j];
					k = j;
				}
			}
			if (cluster[i] == permu[k])
			{
				count++;
			}
		}
		tmp = (double)count / N;
		if (tmp > caval)
			caval = tmp;
	}
	else
	{
		for (i = 0; i < N; i++)
		{
			max = 0;
			k = 0;
			for (j = 0; j < C; j++)
			{
				if (max < U[i][j])
				{
					max = U[i][j];
					k = j;
				}
			}
			if (clust[i] == permu[k])
			{
				count++;
			}
		}
		tmp = (double)count / NLA;
		if (tmp > caval)
			caval = tmp;
	}
}

void trypermu(int k)
{
	int i;
	for (i = 0; i < C; i++)
	{
		if (!Bool[i])
		{
			permu[k] = i;
			Bool[i] = 1;
			if (k == C - 1)
				xuat();
			else
				trypermu(k + 1);
			Bool[i] = 0;
		}
	}
}

double validity_kuhn()
{
	int i, j, k;
	int count = 0, num = 0;
	for (i = 0; i < C; i++)
	{
		Bool[i] = 0;
		permu[i] = i;
	}
	caval = 0;
	trypermu(0);
	return caval;
}

double checkCentroid(double **vold, double **vnew)
{
	int i, j;
	double tg;
	for (j = 0; j < C; j++)
	{
		tg = 0;
		for (i = 0; i < D; i++)
		{
			tg += (double)pow(vnew[j][i] - vold[j][i], 2);
		}
		if (tg >= EPS)
			return tg;
	}
	return -1;
}

void ramdomLabel()
{
	int ita = 0, i, j, k, l, ij;
	int thres, tg, *markcla;
	if (clust == NULL)
		clust = (int *)malloc(N * sizeof(int));

	NLA = (int)(PERCENT_SELECT * N / 100);
	thres = (int)ceil(NLA / C);

	Xla = (double **)malloc(NLA * sizeof(double *));
	Ula = (double **)malloc(NLA * sizeof(double *));
	CLASSLA = (int *)malloc(NLA * sizeof(int));
	for (i = 0; i < NLA; i++)
	{
		Xla[i] = (double *)malloc(D * sizeof(double));
		Ula[i] = (double *)malloc(C * sizeof(double));
	}

	for (i = 0; i < NLA; i++)
	{
		clust[i] = -1;
	}

	// Random select label
	// printf("Begin random select label\n");
	j = 0;
	while (ita < NLA)
	{
		k = rand() % CLASS[j];
		l = 0;
		for (i = 0; i < N; i++)
		{
			if (cluster[i] == j && l < k && clust[i] == -1)
			{
				for (ij = 0; ij < D; ij++)
				{
					Xla[ita][ij] = X[i][ij];
				}
				CLASSLA[ita] = j;
				clust[i] = j;
				ita++;
				l++;
				break;
			}
		}
		j++;
		if (j == C)
			j = 0;
	}

	// Make wrong label
	// printf("Begin make wrong label\n");
	if (PERCENT_WRONG > 0)
	{
		tg = (int)(PERCENT_WRONG * NLA / 100);
		markcla = (int *)malloc(NLA * sizeof(int));
		for (i = 0; i < NLA; i++)
		{
			markcla[i] = 0;
		}
		ij = 0;
		while (ij < tg)
		{
			k = rand() % tg;
			j = 0;
			for (i = 0; i < NLA; i++)
			{
				if (markcla[i] == 0)
				{
					if (j == k)
					{
						if (C == 2)
							CLASSLA[i] = (CLASSLA[i] + 1) % 2;
						else
							CLASSLA[i] = ((rand() % (C - 1) + 1) + CLASSLA[i]) % C;
						markcla[i] = 1;
					}
					else
					{
						j++;
					}
				}
			}
			ij++;
		}
		free(markcla);
	}
}

double FCM(double prepareTime, char *fileout, int *iter)
{
	int i, j, k, l, step = 0, Nrm, numrm;
	double isNext = 1, lambda = 1;
	double tg1 = 0, tgg, tg2 = 0, tg3, runtime = 0;
	double **vold;
	char s[50];
	clock_t t;
	t = clock();
	// find minmax
	double *min, *max;
	double **U1;
	int *markX;
	ramdomLabel();
	Nrm = round(N * PERCENT_SELECT * PERCENT_REMOVE / 10000);

	int *nall = (int *)malloc(NLA * sizeof(int));
	int *nsame = (int *)malloc(NLA * sizeof(int));
	int *ndiff = (int *)malloc(NLA * sizeof(int));
	double radius;

	// printf("finish random label!\n");

	min = (double *)malloc(D * sizeof(double));
	max = (double *)malloc(D * sizeof(double));
	vold = (double **)malloc(C * sizeof(double *));
	U1 = (double **)malloc(N * sizeof(double *));
	markX = (int *)malloc(NLA * sizeof(int));
	for (i = 0; i < N; i++)
	{
		U1[i] = (double *)malloc(D * sizeof(double));
		U1[i][0] = -1;
	}

	// printf("Calculate min and max X\n");

	for (j = 0; j < D; j++)
	{
		min[j] = Xla[0][j];
		max[j] = Xla[0][j];
		for (i = 1; i < NLA; i++)
		{
			if (min[j] > Xla[i][j])
			{
				min[j] = Xla[i][j];
			}
			if (max[j] < Xla[i][j])
			{
				max[j] = Xla[i][j];
			}
		}
	}

	tg1 = -1;
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			tg2 = 0;
			for (l = 0; l < D; l++)
			{
				tg2 += pow(X[i][l] - X[j][l], 2);
			}
			if (tg1 == -1 || tg1 > tg2)
			{
				tg1 = tg2;
			}
		}
	}
	EPS = tg1;
	if (EPS < 0.001)
		EPS = 0.001;
	// printf("eps = %10.5lf",EPS);

	// fill the parameter for label
	tg1 = -1;
	tg2 = -1;
	for (i = 0; i < N; i++)
	{
		for (k = i + 1; k < N; k++)
		{
			tgg = 0;
			for (j = 0; j < D; j++)
			{
				tgg += pow(X[i][j] - X[k][j], 2);
			}
			tgg = sqrt(tgg);
			if (tg1 == -1 || tg1 > tgg)
			{
				tg1 = tgg;
			}
			if (tg2 < tgg)
			{
				tg2 = tgg;
			}
		}
	}
	radius = (tg1 + tg2) * alpha_expert / 2;
	for (i = 0; i < NLA; i++)
	{
		// getting label in radius
		nall[i] = 0;
		nsame[i] = 0;
		ndiff[i] = 0;
		for (k = 0; k < N; k++)
		{
			tgg = 0;
			for (j = 0; j < D; j++)
			{
				tgg += pow(Xla[i][j] - X[k][j], 2);
			}
			if (tgg > 0 && tgg < radius)
			{
				nall[i]++;
			}
		}

		for (k = 0; k < NLA; k++)
		{
			if (i != k)
			{
				tgg = 0;
				for (j = 0; j < D; j++)
				{
					tgg += pow(Xla[i][j] - Xla[k][j], 2);
				}
				if (tgg > 0 && tgg < radius)
				{
					if (CLASSLA[i] == CLASSLA[k])
					{
						nsame[i]++;
					}
					else
					{
						ndiff[i]++;
					}
				}
			}
		}
		nall[i] = nall[i] - nsame[i] - ndiff[i] - 1;
		nall[i] = (int)round(nall[i] * PERCENT_SELECT / 100);
	}

	// randomV
	for (i = 0; i < C; i++)
	{
		vold[i] = (double *)malloc(D * sizeof(double));
		for (j = 0; j < D; j++)
		{
			Vla[i][j] = min[j] + (max[j] - min[j]) * (double)(rand() % 100) / 100;
			vold[i][j] = Vla[i][j];
		}
	}

	// outPutSFCM("demo/init.txt", Vla, Ula, NLA, 0);

	// FCM for label
	do
	{
		// calculate U[i][j]
		for (j = 0; j < C; j++)
		{
			for (i = 0; i < NLA; i++)
			{
				tg2 = 0;
				for (l = 0; l < C; l++)
				{
					tg2 += pow(calcX_subtract_V2(Xla[i], Vla[j]) / calcX_subtract_V2(Xla[i], Vla[l]), 2 / (M - 1));
				}
				Ula[i][j] = 1 / tg2;
			}
		}
		// if(step == 0)	outPutSFCM("demo/cal_u_first.txt", Vla, Ula, NLA, 0);
		// defuzzified label FCM
		for (i = 0; i < NLA; i++)
		{
			tg1 = Ula[i][0];
			tg2 = 0;
			for (j = 1; j < C; j++)
			{
				if (tg1 < Ula[i][j])
				{
					tg1 = Ula[i][j];
					tg2 = j;
				}
			}
			markX[i] = tg2;
		}
		// reduce membership of wrong label
		for (i = 0; i < NLA; i++)
		{
			tg1 = 0;
			tg2 = 0;
			for (j = 0; j < NLA; j++)
			{
				if (i != j && markX[i] == markX[j])
				{
					if (CLASSLA[i] == CLASSLA[j])
					{
						tg1++;
					}
					if (CLASSLA[i] != CLASSLA[j])
					{
						tg2++;
					}
				}
			}
			if (tg2 != 0)
			{
				if (tg1 < (float)(tg1 + tg2) / C)
				{
					tgg = Ula[i][markX[i]] / 2;
					Ula[i][markX[i]] = tgg;
					for (j = 0; j < C; j++)
					{
						if (j != markX[i])
						{
							Ula[i][j] += tgg / (C - 1);
						}
					}
				}
			}
		}

		// calculate V[i][j]
		for (j = 0; j < C; j++)
		{
			for (i = 0; i < D; i++)
			{
				tg1 = 0;
				tg2 = 0;
				for (k = 0; k < NLA; k++)
				{
					tg1 += ((nall[k] + nsame[k]) / (ndiff[k] + 1)) * pow(Ula[k][j], M) * Xla[k][i];
					tg2 += ((nall[k] + nsame[k]) / (ndiff[k] + 1)) * pow(Ula[k][j], M);
				}
				if (tg2 != 0)
					Vla[j][i] = tg1 / tg2;
			}
		}

		// if(step == 0)	outPutSFCM("demo/reduce_u.txt", Vla, Ula, NLA, 0);

		isNext = checkCentroid(vold, Vla);
		if (isNext != -1 && step < MAXSTEPS)
		{
			// Saving old V
			for (j = 0; j < C; j++)
			{
				for (i = 0; i < D; i++)
				{
					vold[j][i] = Vla[j][i];
				}
			}
			// if(step % 10 == 0) printf("FCM first step:   iter %d ok with eps = %10.8lf\n",step,isNext);
			step++;
			// printf("step label=%d\n",step);
		}
		runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();
	} while (isNext != -1);
	printf("FCM finished with %d iter\n", step);
	// outPutSFCM("demo/finish_FCM.txt", Vla, Ula, NLA, step);

	EPS = 0.001;

	// defuzzified label FCM
	for (i = 0; i < NLA; i++)
	{
		tg1 = Ula[i][0];
		tg2 = 0;
		for (j = 1; j < C; j++)
		{
			if (tg1 < Ula[i][j])
			{
				tg1 = Ula[i][j];
				tg2 = j;
			}
		}
		markX[i] = tg2;
	}
	// reduce membership of wrong label
	for (i = 0; i < NLA; i++)
	{
		tg1 = 0;
		tg2 = 0;
		for (j = 0; j < NLA; j++)
		{
			if (i != j && markX[i] == markX[j])
			{
				if (CLASSLA[i] == CLASSLA[j])
				{
					tg1++;
				}
				if (CLASSLA[i] != CLASSLA[j])
				{
					tg2++;
				}
			}
		}
		if (tg2 != 0)
		{
			if (tg1 < (float)(tg1 + tg2) / C)
			{
				tgg = Ula[i][markX[i]] / 2;
				Ula[i][markX[i]] = tgg;
				for (j = 0; j < C; j++)
				{
					if (j != markX[i])
					{
						Ula[i][j] += tgg / (C - 1);
					}
				}
			}
		}
	}
	detectBoundaries();
	adjustMembership();
	applyConstraints();
	//	// Mark wrong label
	//	for(i = 0;i<NLA;i++){
	//    	tg1 = 0;
	//    	tg2 = 0;
	//    	for(j = 0;j<NLA;j++){
	//    		if(i != j && markX[i] == markX[j]){
	//				if(CLASSLA[i] == CLASSLA[j]){
	//					tg1++;
	//				}
	//				if(CLASSLA[i] != CLASSLA[j]){
	//					tg2++;
	//				}
	//			}
	//		}
	//		if(tg2 != 0){
	//			if(tg1 < (float)(tg1 + tg2)/C){
	//				markX[i] = -2;
	//			}
	//		}
	//	}

	//	double **listW = (double**)malloc(NLA*sizeof(double*));
	//	for(i = 0;i<NLA;i++){
	//		listW[i] = (double*)malloc(2*sizeof(double));
	//		listW[i][0] = i;
	//		listW[i][1] = (double)(nall[i] + nsame[i]) / (ndiff[i] + 1);
	//	}

	//	for(i = 0;i<NLA-1;i++){
	//		for(j = i+1;j<NLA;j++){
	//			if(listW[i][1] < listW[j][1]){
	//				tgg = listW[i][1];
	//				listW[i][1] = listW[j][1];
	//				listW[j][1] = tgg;
	//
	//				tgg = listW[i][0];
	//				listW[i][0] = listW[j][0];
	//				listW[j][0] = tgg;
	//			}
	//		}
	//	}
	//	numrm = 0;
	//	for(i = 0;i<NLA && numrm < Nrm;i++){
	//		if(markX[(int)listW[i][0]] == -2){
	//			numrm++;
	//			U1[(int)listW[i][0]][0] = -2;
	//		}
	//	}
	//
	// printf("calculate U1 for all\n");
	// calculate U1 for all
	l = 0;
	for (i = 0; i < N; i++)
	{
		if (clust[i] >= 0 && U1[i][0] != -2)
		{
			for (j = 0; j < C; j++)
			{
				U1[i][j] = Ula[l][j];
			}
			l++;
		}
		else
		{
			for (j = 0; j < C; j++)
			{
				tg2 = 0;
				for (l = 0; l < C; l++)
				{
					tg2 += pow(calcX_subtract_V2(X[i], Vla[j]) / calcX_subtract_V2(X[i], Vla[l]), 2 / (M - 1));
				}
				U1[i][j] = 1 / tg2;
			}
		}
	}
	// outPutSFCM("demo/cal_u_init.txt", Vla, U1, N, 0);
	//  find min max
	for (j = 0; j < D; j++)
	{
		min[j] = X[0][j];
		max[j] = X[0][j];
		for (i = 1; i < N; i++)
		{
			if (min[j] > X[i][j])
			{
				min[j] = X[i][j];
			}
			if (max[j] < X[i][j])
			{
				max[j] = X[i][j];
			}
		}
	}

	// random V
	for (i = 0; i < C; i++)
	{
		for (j = 0; j < D; j++)
		{
			V[i][j] = min[j] + (max[j] - min[j]) * (double)(rand() % 100) / 100;
			vold[i][j] = V[i][j];
		}
	}

	printf("Semisupervise clustering\n");
	// Semisupervise clustering
	step = 0;
	isNext = 1;
	do
	{
		// calculate U[i][j]
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < C; j++)
			{
				tg2 = 0;
				tg3 = calcX_subtract_V2(X[i], V[j]);
				for (l = 0; l < C; l++)
				{
					if (l == j)
						tg2 += 1;
					else
					{
						tg1 = calcX_subtract_V2(X[i], V[l]);
						if (tg1 != 0)
						{
							tg1 = pow(tg3 / tg1, 2 / (M - 1));
						}
						else
						{
							tg1 = 0;
						}
						tg2 += tg1;
					}
				}
				U[i][j] = (U1[i][j] * lambda + 1 / tg2) / (1 + lambda);
				if (U[i][j] >= 1.0)
				{
					printf("U = %5.4lf > 1 U1=%5.4lf tg2=%5.4lf\n", U[i][j], U1[i][j], tg2);
					// getch();
				}
			}
		}

		// calculate V[i][j]
		for (j = 0; j < C; j++)
		{
			for (i = 0; i < D; i++)
			{
				tg1 = 0;
				tg2 = 0;
				for (k = 0; k < N; k++)
				{
					tg1 += (pow(U[k][j], M) + lambda * pow(U[k][j] - U1[k][j], M)) * X[k][i];
					tg2 += (pow(U[k][j], M) + lambda * pow(U[k][j] - U1[k][j], M));
				}
				V[j][i] = tg1 / tg2;
			}
		}

		// if(step == 0) outPutSFCM("demo/ts3fcm_first.txt", V, U, N, 0);

		// printf("check centroid\n");
		isNext = checkCentroid(vold, V);
		if (isNext != -1 && step < MAXSTEPS)
		{
			// Saving old V
			for (j = 0; j < C; j++)
			{
				for (i = 0; i < D; i++)
				{
					vold[j][i] = V[j][i];
				}
			}
			// if(step % 10 == 0) printf("TS3FCM:   iter %d ok with eps = %10.8lf\n",step,isNext);
			step++;
		}
		else
		{
			isNext = -1;
		}
		runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();
	} while (isNext != -1);
	// if(step == 0)
	// outPutSFCM("demo/sfcm_final.txt", V, U, N, 0);

	// printf("SFCM:   Finish with iter %d\n",step);

	printf("AS3FCPC\n");
	// Semisupervise clustering
	step = 0;
	isNext = 1;
	do
	{
		// calculate U[i][j]
		for (i = 0; i < N; i++)
		{

			for (k = 0; k < C; k++)
			{
				// --- Distance squared ---
				double dik2 = calcX_subtract_V2(X[i], V[k]); // d_{ik}^2
				if (dik2 == 0)
					dik2 = 1e-6; // avoid division by zero

				// --- Numerator terms ---
				double term1 = 0.0, term2 = 0.0, term3 = 0.0;
				double sum_dik2_inv = 0.0;

				// Term 1: 2 * U1[i][k] * dik^2
				term1 = 2.0 * U1[i][k] * dik2;

				// Subtract: alpha * sum_{(i,j) in M} sum_{l!=k} U[j][l]
				for (int mj = 0; mj < num_must_link[i]; mj++)
				{
					int j = must_link[i][mj];
					for (int l = 0; l < C; l++)
					{
						if (l != k)
						{
							term1 -= alpha * U[j][l];
						}
					}
				}

				// Subtract: alpha * sum_{(i,j) in C} U[j][k]
				for (int cj = 0; cj < num_cannot_link[i]; cj++)
				{
					int j = cannot_link[i][cj];
					term1 -= alpha * U[j][k];
				}

				// --- Term 2: sum over all clusters ---
				double numerator2 = 0.0;
				sum_dik2_inv = 0.0;
				for (int kk = 0; kk < C; kk++)
				{
					double temp = 2.0 * U1[i][kk];
					for (int mj = 0; mj < num_must_link[i]; mj++)
					{
						int j = must_link[i][mj];
						for (int l = 0; l < C; l++)
						{
							if (l != kk)
							{
								temp -= alpha * U[j][l] / dik2;
							}
						}
					}
					numerator2 += temp;
					sum_dik2_inv += 1.0 / dik2;
				}
				if (sum_dik2_inv != 0.0)
					term2 = numerator2 / sum_dik2_inv;
				else
					term2 = 0.0;

				// --- Term 3: cannot-link penalty ---
				double numerator3 = 0.0;
				for (int kk = 0; kk < C; kk++)
				{
					for (int cj = 0; cj < num_cannot_link[i]; cj++)
					{
						int j = cannot_link[i][cj];
						numerator3 += alpha * U[j][kk] / dik2;
					}
				}
				if (sum_dik2_inv != 0.0)
					term3 = numerator3 / sum_dik2_inv;
				else
					term3 = 0.0;

				// --- Combine all terms ---
				double numerator = term1 - term2 + term3 + 4.0;
				double denominator = 4.0 * dik2;

				U[i][k] = numerator / denominator;

				// Clamp U[i][k] to [0, 1]
				if (U[i][k] < 0.0)
					U[i][k] = 0.0;
				if (U[i][k] > 1.0)
					U[i][k] = 1.0;
			}
		}

		// calculate V[i][j]
		for (j = 0; j < C; j++)
		{
			for (i = 0; i < D; i++)
			{
				tg1 = 0;
				tg2 = 0;
				for (k = 0; k < N; k++)
				{
					tg1 += (pow(U[k][j], M) + lambda * pow(U[k][j] - U1[k][j], M)) * X[k][i];
					tg2 += (pow(U[k][j], M) + lambda * pow(U[k][j] - U1[k][j], M));
				}
				V[j][i] = tg1 / tg2;
			}
		}

		// if(step == 0) outPutSFCM("demo/ts3fcm_first.txt", V, U, N, 0);

		// printf("check centroid\n");
		isNext = checkCentroid(vold, V);
		if (isNext != -1 && step < MAXSTEPS)
		{
			// Saving old V
			for (j = 0; j < C; j++)
			{
				for (i = 0; i < D; i++)
				{
					vold[j][i] = V[j][i];
				}
			}
			// if(step % 10 == 0) printf("TS3FCM:   iter %d ok with eps = %10.8lf\n",step,isNext);
			step++;
		}
		else
		{
			isNext = -1;
		}
		runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();
	} while (isNext != -1);
	*iter = step;
	// if(step == 0)
	// outPutSFCM("demo/sfcm_final.txt", V, U, N, 0);

	// printf("SFCM:   Finish with iter %d\n",step);
	outPutV(fileout, step, runtime, 1);
	printf("AS3FCPC:   OutputV ok with %ld steps numrm = %d/%d\n", step, numrm, Nrm);

	// free(listW);
	free(min);
	free(max);
	free(markX);
	free(U1);
	free(vold);
	free(Ula);
	// printf("free 1/2\n");
	for (i = 0; i < NLA; i++)
	{
		free(Xla[i]);
	}

	free(CLASSLA);
	return runtime;
}

int main(int argc, char **argv)
{

	char s1[50], s[50];
	FILE *f;
	int i, n = 10, iter, averageIter = 0, count = 0;
	;
	double runtime = 0, runti, ifv = 0, db = 0, ma = 0, ca = 0, nmi1 = 0, ti, ta, tg, runtime1;
	double db1, ri1, asw1, asw, pbm1, pbm;
	double ma1 = 0, ma2 = 0, ri = 0, ca1 = 0, tgg;
	if (label == NULL)
		label = (int *)malloc(N * sizeof(int));
	//	ALPHA = 0.5;
	//	EPS = 0.01;

	srand(time(NULL));

	if (argc == 4)
	{
		sprintf(s1, "data/%s.txt", argv[1]);
		n = atoi(argv[2]);
		if (n < 1)
			n = 1;
		PERCENT_WRONG = atoi(argv[3]);
		// ALPHA = (double)atof(argv[3]);
	}
	else
	{
		printf("invalid input parameter ifs\n");
		exit(0);
	}
	// Allocate must_link and cannot_link
	must_link = (int **)malloc(N * sizeof(int *));
	cannot_link = (int **)malloc(N * sizeof(int *));
	for (i = 0; i < N; i++)
	{
		must_link[i] = (int *)calloc(N, sizeof(int));
		cannot_link[i] = (int *)calloc(N, sizeof(int));
	}

	runti = input(s1);
	//	normalizeX();
	//	printf("input ok!\n");
	sprintf(s1, "result/As3fcpc/As3fcpc_result_%s_%d.csv", argv[1], PERCENT_WRONG);
	f = fopen(s1, "w");

	if (!f)
	{
		printf("Can't open file %s!", s1);
		exit(0);
	}
	fprintf(f, "Result_file_%s:\n", argv[1]);
	fprintf(f, "Time,ca,ca_label,db,asw,pbm,runtime\n");
	for (i = 1; i <= n; i++)
	{
		printf("\nTime %d :\n", i);
		sprintf(s, "result/1/ts3fcm_%s_%d_time%d.txt", argv[1], PERCENT_WRONG, i);
		runtime1 = FCM(runti, s, &iter);
		runtime += runtime1;

		db1 = DB(s);
		db += db1;
		tgg = validity_kuhn();
		ca1 += tgg;
		free(clust);
		clust = NULL;
		tg = validity_kuhn();
		ca += tg;
		asw1 = ASWC(s);
		asw += asw1;
		pbm1 = PBM(s);
		pbm += pbm1;

		fprintf(f, "%d,%10.5lf,%10.5lf,%10.5lf,%10.5lf,%10.5lf,%10.5lf\n", i, tg, tgg, db1, asw1, pbm1, runtime1);

		printf("ca = %10.5lf asw=%10.5lf pbm=%10.5lf runtime=%10.5lf\n", tg, asw1, pbm1, runtime1);
		printf("ca_label = %10.5lf DB = %10.5lf", tgg, db1);
		// printf("RI=%10.5lf",ri1);
		averageIter += iter;
	}

	fprintf(f,"\nRI: %15.10f\n",ri/n);
	fprintf(f,"NMI: %15.10f\n",tg);
	fprintf(f, "DB: %15.10f\n", db / n);
	fprintf(f, "F1-score: %15.10f\n", F1_SCORE(s1) / n);
//	fprintf(f, "ASW: %15.10f\n", asw / n);
//	fprintf(f, "PBM: %15.10f\n", pbm / n);
	//	if(count){
	//		tg = nmi1/count;
	//	}else tg = -1;
	//	fprintf(f,"IFV: %15.10f\n",ifv/n);
//	fprintf(f, "Ca: %15.10f\n", ca / n);
//	fprintf(f, "Ca_label: %15.10f\n", ca1 / n);
//	fprintf(f, "Time: %15.10f\n", (runtime / n) + runti);
	fprintf(f, "average_iteration: %15.10f\n", (double)averageIter / n);

	fclose(f);
	free(X);
	free(U);
	free(V);
	free(cluster);
	// Free must_link and cannot_link
	for (i = 0; i < N; i++)
	{
		free(must_link[i]);
		free(cannot_link[i]);
	}
	free(must_link);
	free(cannot_link);
	free(label);
	return 0;
}
