#include"mylib.h"
int PERCENT_SELECT = 20;
int PERCENT_WRONG = 0;

int Bool[100];
int permu[100];
double caval = 0;
int neighbor = 5;
double lamda1 = 1;
double lamda2 = 10;
int *clust = NULL;

void xuat() {
	double max, tmp;
	int i,j,k,count = 0;
	
	if(clust == NULL){	
		for(i = 0;i<N;i++){
			max = 0;
			k = 0;
			for(j = 0;j<C;j++){
				if(max < U[i][j]){
					max = U[i][j];
					k = j;	
				}			
			}
			if(cluster[i] == permu[k]){
				count++;
			}
		}
		tmp = (double)count/N;
		if(tmp > caval) caval = tmp;
	}else{
		for(i = 0;i<N;i++){
			max = 0;
			k = 0;
			for(j = 0;j<C;j++){
				if(max < U[i][j]){
					max = U[i][j];
					k = j;	
				}			
			}
			if(clust[i] == permu[k]){
				count++;
			}
		}
		tmp = (double)count/NLA;
		if(tmp > caval) caval = tmp;
	}
}
 
void trypermu(int k) {	
	int i;
	for (i = 0; i < C; i++) {
		if (!Bool[i]) {
			permu[k] = i;
			Bool[i] = 1;
			if (k == C - 1) xuat();
			else trypermu(k + 1);
			Bool[i] = 0;
		}
	}
}

double validity_kuhn(){
	int i,j,k;
	int count = 0, num = 0;
	for(i = 0;i<C;i++){
		Bool[i] = 0;
		permu[i] = i;	
	}
	caval = 0;
	trypermu(0);
	return caval;
}

double checkCentroid(double **vold, double **vnew){
	int i,j;
	double tg;
	for(j = 0;j<C;j++){
		tg = 0;
		for(i = 0;i<D;i++){
		   tg += (double)pow(vnew[j][i] - vold[j][i],2);       
		}
		if (tg>= EPS) return tg;
	}
	return -1;
}

void ramdomLabel(){
	int ita = 0,i,j,k,l,ij;
	int thres,tg,*markcla;
	if(clust == NULL) clust = (int*)malloc(N*sizeof(int)); 
	 
	NLA = (int)(PERCENT_SELECT * N / 100);	
	thres = (int)ceil(NLA/C);
	
	Xla = (double**)malloc(NLA*sizeof(double*));
	Ula = (double**)malloc(NLA*sizeof(double*));
	CLASSLA = (int*)malloc(NLA*sizeof(int));
	for(i = 0;i<NLA;i++){
		Xla[i] = (double*)malloc(D*sizeof(double));
		Ula[i] = (double*)malloc(C*sizeof(double));
	}
	
	for(i = 0;i<N;i++){
		clust[i] = -1;	
	}
		
	// Random select label
	//printf("Begin random select label\n");
	j = 0;
	while(ita < NLA){
		k = rand() % CLASS[j];
		l = 0;
		for(i = 0;i<N;i++){
			if(cluster[i] == j && l < k && clust[i] == -1){
				for(ij = 0;ij < D;ij++){
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
		if(j == C) j = 0;
	}
	
	// Make wrong label
	//printf("Begin make wrong label\n");
	if(PERCENT_WRONG > 0){
		tg = (int)(PERCENT_WRONG * NLA / 100);
		markcla = (int*)malloc(NLA*sizeof(int)); 
		for(i = 0;i<NLA;i++){
			markcla[i] = 0;
		}
		ij = 0;		
		while(ij < tg){
			k = rand() % tg;		
			j = 0;
			for(i = 0;i<NLA;i++){
				if(markcla[i] == 0){										
					if(j == k){
						if(C == 2) CLASSLA[i] = (CLASSLA[i] + 1) % 2;
						else CLASSLA[i] = ((rand() % (C - 1)	+ 1) + CLASSLA[i]) % C;
						markcla[i] = 1;
					}else{
						j++;
					} 
				}							
			}	
			ij++;
		}
		free(markcla);
	}
}

double calculate_P(int *clust, double *sk, double **weight, int i, int j){
	int k,l;
	double tgg,tg1,tg2,tg3,tg4;	
	tgg = calcX_subtract_V2(X[i],V[j]);
	if(clust[i] == j){
		tg1 = lamda1 * sk[i] * tgg;							
	}else{
		tg1 = 0;
	}            		
	tg2 = 0;
	for(l = 0;l<N;l++){
		if(clust[l] == -1){                   
			tg2 += weight[i][l] * U[i][j];			
		}
	}
	if(tg1 + sk[i] != 0 && sk[i] != 0) return tg2 * lamda2 / sk[i];
	else return 0;
	
}

double calculate_Q(int *clust, double *sk, double **weight, int i, int j){
	int k,l;
	double tgg,tg1,tg2,tg3,tg4;	
	tgg = calcX_subtract_V2(X[i],V[j]);
	tg1 = tgg + lamda1 * sk[i] * tgg;	
	tg2 = 0;
	for(l = 0;l<N;l++){
		if(clust[l] == -1){                   
			tg2 += weight[i][l];
		}
	}
	if(tg1 + sk[i] != 0 && sk[i] != 0) return tg2 * lamda2 / sk[i];
	else return 0;
}

double calculate_Z(int *clust, double *sk, double **weight, int i, int j){
	int k,l;
	double tgg,tg1,tg2,tg3,tg4;	
	tg1 = 0;
	for(l = 0;l<N;l++){
		if(clust[l] > -1 && sk[i] != 0){                   
			tg1 += weight[i][l] * U[l][j] / sk[l];			
		}
	}
	return lamda2 * tg1;
}

double calculate_T(int *clust, double *sk, double **weight, int i, int j){
	int k,l;
	double tgg,tg1,tg2,tg3,tg4;	
	tg1 = calcX_subtract_V2(X[i],V[j]);
	for(l = 0;l<N;l++){
		if(clust[l] > -1 && sk[i] != 0){                   
			tg1 += lamda2 * weight[i][l] / sk[i];			
		}
	}
	return tg1;
}

double cs3fcm(double prepareTime, char *fileout, int *iter){
	int i,j,k,l,step = 0;
	double isNext = 1, lambda = 1;
	double tg1 = 0,tgg,tg2 = 0, tg3,tg4,runtime = 0;
	double **vold;
	char s[50];
	clock_t t;
	t = clock();        
	// find minmax
	double *min, *max;
	double **U1;
	ramdomLabel();
	long naver;
	
	double radius, sigma1 = -1;
	
	//printf("finish random label!\n");
	
	min = (double*)malloc(D*sizeof(double));
	max = (double*)malloc(D*sizeof(double));
	vold = (double**)malloc(C*sizeof(double*));
	double **per = (double**)malloc(C*sizeof(double*));
	double *sk = (double*)malloc(N*sizeof(double));
	int *ll1 = (int*)malloc(C*sizeof(int));
	int *ll2 = (int*)malloc(C*sizeof(int));
	int *ll3 = (int*)malloc(C*sizeof(int));
	int *ll4 = (int*)malloc(C*sizeof(int));
	for(j = 0;j<C;j++){
		per[j] = (double*)malloc(C*sizeof(double));
		ll1[j] = 0;
		ll2[j] = 0;
		ll3[j] = j;
		ll4[j] = j;
		for(i = 0;i<C;i++){
			per[j][i] = 0.0;
		}
	}
	
	int *markX = (int*)malloc(N*sizeof(int));
	
	int **nei = (int**)malloc(NLA * sizeof(int*));
	double *distan = (double*)malloc(N * sizeof(double));
	double **weight = (double**)malloc(N * sizeof(double*));
	for(i = 0;i<NLA;i++){
		nei[i] = (int*)malloc((neighbor + 1) * sizeof(int));
		for(j = 0;j<neighbor + 1;j++){
			nei[i][j] = -1;
		}
	}
	for(i = 0;i<N;i++){
		weight[i] = (double*)malloc(N * sizeof(double));
		for(j = 0;j<N;j++){
			weight[i][j] = 0;
		}
	}
			
	//printf("Calculate min and max X\n");    	
	for(j = 0;j<D;j++){
		min[j] = Xla[0][j];
		max[j] = Xla[0][j];
		for(i = 1;i<N;i++){
			if(min[j] > X[i][j]){
				min[j] = X[i][j];
			}
			if(max[j] < X[i][j]){
				max[j] = X[i][j];
			}
		}
	}
			
	// randomV
	for(i = 0;i<C;i++){
		vold[i] = (double*)malloc(D*sizeof(double));  
		for(j = 0;j<D;j++){
			V[i][j] = min[j] + (max[j] - min[j]) * (double)(rand() % 100) / 100;
			vold[i][j] = V[i][j];			
		}
	}
	free(min);
	free(max);
	//calculate average distance
	naver = 0;
	tg1 = 0;
	for(i = 0;i<N - 1;i++){
		for(j = i+1;j<N;j++){
			tg2 = 0;
			for(l = 0;l < D;l++){
				tg2 += pow(X[i][l] - X[j][l],2);
			}
			//printf("%10.5lf\n",sqrt(tg2));
			tg1+=sqrt(tg2);			
		}
	}
	sigma1 = 2*tg1/(N*(N-1));	
	//Calculate the neighbor of label
	k = 0;
	for(i = 0;i<N;i++){
		if(clust[i] > -1){
			nei[k][0] = i;			
			for(j = 0;j<N;j++){
				if(j != i){
					tg2 = 0;
					for(l = 0;l < D;l++){
						tg2 += pow(X[i][l] - X[j][l],2);
					}
					distan[j] = tg2;
				}
			}			
			for(l = 1;l<=neighbor;l++){
				for(j = 0;j<N;j++){
					if(j != i){
						if(l == 1){
							if(nei[k][l] == -1){
								nei[k][l] = j;	
							}else{
								if(distan[j] <= distan[nei[k][l]] && j != nei[k][l]){
									nei[k][l] = j;	
								}
							}							
						}else{
							if(nei[k][l] == -1 && distan[j] >= distan[nei[k][l-1]] && j != nei[k][l-1]){
								nei[k][l] = j;
							}else{
								if(nei[k][l] > -1 && distan[j] <= distan[nei[k][l]] && distan[j] >= distan[nei[k][l-1]] && j != nei[k][l] && j != nei[k][l-1]){
									nei[k][l] = j;	
								}
							}
						}
					}
				}	
			}							
			k++;
		}
	}
		
	// FCM for all data
	do{       
		// calculate U[i][j]  
		for(j = 0;j<C;j++){              
			for(i = 0;i<N;i++){
				tg2 = 0;
				for(l = 0;l<C;l++){                        
					tg2 += pow(calcX_subtract_V2(X[i],V[j])/calcX_subtract_V2(X[i],V[l]),2/(M-1));
				}                        
				U[i][j] = 1/tg2;   
			}   
		}
				 
		//calculate V[i][j]  
		for(j = 0;j<C;j++){
			for(i = 0;i<D;i++){
				tg1 = 0;
				tg2 = 0;  
				for(k = 0;k<N;k++){
					tg1 += pow(U[k][j],M) * X[k][i];
					tg2 += pow(U[k][j],M);        
				}
				V[j][i] = tg1/tg2;
			}               
		}     
		isNext = checkCentroid(vold,V);
		if(isNext != -1 && step < MAXSTEPS){
			// Saving old V
			for(j = 0;j<C;j++){
				for(i = 0;i<D;i++){
					vold[j][i] = V[j][i];       
				}
			}
			if(step % 10 == 0) printf("FCM:   iter %d ok with eps = %10.8lf\n",step,isNext);
			step++;  
			//printf("step label=%d\n",step);          
		} 
		runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();                        
	}while(isNext != -1);
	//outPutSFCM("demo/cs3fcm_finish_fcm.txt", V, U, N, step);
	// defuzzied FCM
	for(i = 0;i<N;i++){
		tg1 = U[i][0];
		tg2 = 0;
		for(j = 1;j<C;j++){
			if(tg1 < U[i][j]){
				tg1 = U[i][j];
				tg2 = j;
			}
		}
		markX[i] = tg2;
		if(markX[i] >= C) markX[i] = 0;
	}
	
	for(i = 0;i<N;i++){
		if(clust[i] > -1){
			for(j = 0;j<C;j++){		
				if(markX[i] == j){
					ll2[j]++;					
				}
				if(clust[i] == j){
					ll1[j]++;
				}	
			}
		}
	}
	
	for(i = 0;i<C - 1;i++){
		for(j = i+1;j<C;j++){
			if(ll1[i] < ll1[j]){
				k = ll1[i];
				ll1[i] = ll1[j];
				ll1[j] = k;
				k = ll3[i];
				ll3[i] = ll3[j];
				ll3[j] = k;
			}	
			if(ll2[i] < ll2[j]){
				k = ll2[i];
				ll2[i] = ll2[j];
				ll2[j] = k;
				k = ll4[i];
				ll4[i] = ll4[j];
				ll4[j] = k;
			}			
		}
	}
	
	for(j=0;j<C;j++){
		for(i = 0;i<N;i++){			
			if(cluster[i] == ll3[j]){			
				if(clust[i] == ll3[j]){
					clust[i] = ll4[j] + C;					 
				}
				cluster[i] = ll4[j] + C;				
			}
		}		 
	}
	
	for(i = 0;i<N;i++){
		if(clust[i] >= C){
			clust[i] = clust[i] - C;					 
		}
		if(cluster[i] >= C){
			cluster[i] = cluster[i] - C;			
		}		
	}
		
	//Compute weight Sk
	for(i = 0;i<N;i++){
		if(clust[i] > -1){						
			per[clust[i]][markX[i]] = per[clust[i]][markX[i]] + 1;			
		}
	}
	for(i = 0;i<C;i++){
		tg3 = 0;
		for(j = 0;j<C;j++){
			tg3 += per[i][j];
		}
		if(tg3 > 0){		
			//printf("tg3 = %10.5lf\n",tg3);
			for(j = 0;j<C;j++){
				per[i][j] = (double)per[i][j] / tg3;
				//printf("per[%d][%d] = %10.5lf\n",i,j,per[i][j]);
			}
		}		
	}	
//	for(i = 0;i<C;i++){
//		for(j = 0;j<C;j++){
//			printf("%5.3lf ",per[i][j]);	
//		}
//		printf("\n");
//	}
	
	runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
	t = clock();
	//printf("Compute Sk begin runtime = %10.5lf!\n",runtime);
	
	for(i = 0;i<N;i++){
		//printf("i = %d markX=%d clust=%d\n",i,markX[i],clust[i]);
		k = markX[i];
		l = clust[i];
		if(k == l){
			//printf("i=%d bang nhau clust=%d markX=%d u = %10.5lf",i,l,k,U[i][l]); 
			//printf(" per = %10.5lf\n",per[k][l]);
			sk[i] = per[k][l] * U[i][l];	
		}else{
			if(l > -1){
				//printf("i=%d khac nhau clust=%d markX=%d u = %10.5lf",i,l,k,U[i][l]); 
				//printf(" per = %10.5lf\n",per[k][l]);
				sk[i] = per[k][l] * (1 - U[i][l]);	
				//printf("cal sk ok\n");
			}else{
				sk[i] = 0;
			}
		}
		//printf("i = %d sk=%10.5lf markX=%d clust=%d\n",i,sk[i],markX[i],clust[i]);
	}	
	//printf("Compute Sk complete!\n");	
	//getch();
	
	// construct graph weight
	for(i = 0;i<NLA;i++){
		for(j = 1;j<=neighbor;j++){
			tg2 = 0;
			for(l = 0;l < D;l++){
				tg2 += pow(X[nei[i][0]][l] - X[nei[i][j]][l],2);
			}
			//printf("\ni = %d  tg2 = %10.5lf sigma1=%10.5lf weight=",i,tg2,sigma1);
			weight[nei[i][0]][nei[i][j]] = exp(-tg2/pow(sigma1,2));
			//printf("%10.5lf ",weight[nei[i][0]][nei[i][j]]);
		}
		//printf("\n");
	}	
	//printf("Compute weight complete!\n");
	
	// Initial center V of label
	for(j = 0;j<C;j++){
		for(l = 0;l<D;l++){
			V[j][l] = 0;
		}
		k = 0;
		for(i = 0;i<N;i++){	
			if(clust[i] == j){
				for(l = 0;l<D;l++){
					V[j][l] = ((double)k / (k + 1)) * V[j][l] + X[i][l] / (k + 1);
				}
				k++;	
			}				
		}
		for(l = 0;l<D;l++){
			vold[j][l] = V[j][l];
		}
	}
	
	runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
	t = clock();
	step = 0;
	isNext = 1;	
	//Cs3FCM
	do{       
		// calculate U[i][j]  
		for(j = 0;j<C;j++){              
			for(i = 0;i<N;i++){
				if(clust[i] > -1){
					tg1 = calculate_P(clust,sk,weight,i,j);
					tg2 = calculate_Q(clust,sk,weight,i,j);	
					tg3 = 0;
					tg4 = 0;
					for(l = 0;l<C;l++){
						tgg = calculate_Q(clust,sk,weight,i,l);
						if(tgg != 0){
							tg3 += calculate_P(clust,sk,weight,i,l) / tgg;
							tg4 += 1 / tgg;	
						} 
					}					
					//printf("label tg1=%10.5lf tg2=%10.5lf tg3=%10.5lf tg4=%10.5lf\n\n",tg1,tg2,tg3,tg4);
					if(tg4 != 0 && tg2 != 0){
						U[i][j] = (tg1 + (1 - tg3) / tg4) / tg2;
					}else{
						U[i][j] = 1.0 / C;
					}
				}else{
					tg1 = calculate_Z(clust,sk,weight,i,j);
					tg2 = calculate_T(clust,sk,weight,i,j);
					tg3 = 0;
					tg4 = 0;
					for(l = 0;l<C;l++){
						tgg = calculate_T(clust,sk,weight,i,l);
						if(tgg != 0){
							tg3 += calculate_Z(clust,sk,weight,i,l) / tgg;
							tg4 += 1 / tgg;
						}
					}
					//printf("unlabel tg1=%10.5lf tg2=%10.5lf tg3=%10.5lf tg4=%10.5lf\n\n",tg1,tg2,tg3,tg4);
					if(tg4 != 0 && tg2 != 0){
						U[i][j] = (tg1 + (1 - tg3) / tg4) / tg2;
					}else{
						U[i][j] = 1.0 / C;
					}
				}
				//printf("%10.5lf ",U[i][j]);
			} 
			//printf("\n");  
		}
		
		//calculate V[i][j]  
		for(j = 0;j<C;j++){
			for(i = 0;i<D;i++){
				tg1 = 0;
				tg2 = 0;
				tg3 = 0;
				tg4 = 0;  
				for(k = 0;k<N;k++){
					tg1 += pow(U[k][j],M) * X[k][i];
					tg2 += pow(U[k][j],M);        
					if(clust[k] > -1){
						if(clust[k] == j){
							tg3 += sk[k] * pow(U[k][j] - 1,2) * X[k][i];
							tg4 += sk[k] * pow(U[k][j] - 1,2);
						}else{
							tg3 += sk[k] * pow(U[k][j],2) * X[k][i];
							tg4 += sk[k] * pow(U[k][j],2);
						}
					}
				}
				tg1 += lamda1 * tg3;
				tg2 += lamda1 * tg4;
				if(tg2 != 0) V[j][i] = tg1/tg2;
				else{
					printf("tg1=%10.5lf  tg2=%10.5lf\n",tg1,tg2);
					//getch();
				}
				//printf("%10.5lf(%10.5lf,%10.5lf)\n",V[j][i],tg1,tg2);
			}               
			//printf("\n");
		}     
		//getch();
		//if(step == 0) outPutSFCM("demo/cs3fcm_calc_u_first.txt", V, U, N, step);		 
		isNext = checkCentroid(vold,V);
		//printf("isNext=%5.3lf  step=%d\n",isNext,step);
		if(isNext >= 0 && step < MAXSTEPS){
			// Saving old V
			//printf("Saving old V\n");
			for(j = 0;j<C;j++){
				for(i = 0;i<D;i++){
					vold[j][i] = V[j][i];       
				}
			}
			if(step % 10 == 0) printf("CS3FCM:   iter %d ok with eps = %10.8lf\n",step,isNext);
			step++;  
			//printf("step label=%d\n",step);          
		} 
		runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();                        
	}while(isNext >= 0);	
	//outPutSFCM("demo/cs3fcm_finish.txt", V, U, N, step);		 
	outPutV(fileout, step, runtime,1); 
	printf("SFCM:   OutputV ok with %ld steps\n",step);
	
	runtime += (double)(clock() - t) / CLOCKS_PER_SEC;
//	for(i = 0;i<C;i++){
//		free(per[i]);	
//	}
	free(per);	
//	for(i = 0;i<NLA;i++){
//		free(nei[i]);	
//	}
	free(nei);
	free(weight);
	free(ll1);
	free(ll2);
	free(ll3);
	free(ll4);
	free(distan);
	free(markX);	
	free(Xla);
	free(Ula);
	free(CLASSLA);
	
	return runtime;
}

int main(int  argc, char **argv){
	
	char s1[50],s[50];
	FILE *f;
	int i,n = 10,iter,averageIter = 0,count = 0;;
	double runtime = 0,runti,ifv = 0,db = 0,ma = 0,ca=0,nmi1 = 0,ti,ta,tg,runtime1;
	double db1,ri1,asw1,asw,pbm1,pbm;
	double ma1 = 0,ma2 = 0,ri = 0,ca1=0,tgg;
//	ALPHA = 0.5;
//	EPS = 0.01;
	
	srand(time(NULL));
	
	if (argc >= 3){ 
		sprintf(s1,"data/%s.txt",argv[1]);
		n = atoi(argv[2]);
		if(n < 1) n = 1; 
		PERCENT_WRONG = atoi(argv[3]);
		//ALPHA = (double)atof(argv[3]);		
				
	}else{
		printf("invalid input parameter ifs\n");
		exit(0);
	}
			
	runti = input(s1);	
//	normalizeX();
//	printf("input ok!\n");
	sprintf(s1,"result/cs3fcm/CS3FCM_result_%s_%d.csv",argv[1],PERCENT_WRONG);
	f = fopen(s1,"w");
	if(!f){
		printf("Can't open file %s!",s1);		
		exit(0);
	}
	fprintf(f,"Result_file_%s:\n",argv[1]);
	fprintf(f,"Time,ca,ca_label,db,asw,pbm,runtime\n");
	for(i = 1;i<=n;i++){		
		printf("\nTime %d :\n",i);
		sprintf(s,"result/1/CS3FCM_%s_%d_time%d.txt",argv[1],PERCENT_WRONG,i);
		runtime1 = cs3fcm(runti,s,&iter);		
		runtime += runtime1; 
		
//		printf("calculate DB\n");
		db1 = DB(s);
		db += db1;
//		ifv += IFV(s1);		
//		printf("calculate RI\n");
		//ri1 = RI(s);	
		//ri += ri1;
//		ta = nmi(s1,&ti); 
//		if(ta != -1000){
//			count++;
//			nmi1 += ta; 
//		}
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
		
		fprintf(f,"%d,%10.5lf,%10.5lf,%10.5lf,%10.5lf,%10.5lf,%10.5lf\n",i,tg,tgg,db1,asw1,pbm1,runtime1);
		//fprintf(f,"%d,%10.5lf,%10.5lf,%10.5lf,%10.5lf\n",i,tg,tgg,db1,runtime1);
				
		printf("ca = %10.5lf asw=%10.5lf pbm=%10.5lf runtime=%10.5lf\n",tg,asw1,pbm1,runtime1);
		printf("ca_label = %10.5lf  DB = %10.5lf",tgg,db1);
		//printf("RI=%10.5lf",ri1);
		averageIter += iter;	
	}
	
	
	//fprintf(f,"\nRI: %15.10f\n",ri/n);

	// Print clustering metrics
	double ri = RI(s1);
	double f1 = F1_SCORE();
	double nmi_time = 0.0;
	double nmi_val = nmi(s1, &nmi_time);
	double db_val = DB(s1);
	printf("RI = %10.5lf\n", ri);
	printf("F1 = %10.5lf\n", f1);
	printf("NMI = %10.5lf\n", nmi_val);
	printf("DB = %10.5lf\n", db_val);
	// fprintf(f,"DB: %15.10f\n",db/n);
	// fprintf(f,"ASW: %15.10f\n",asw/n);	
	// fprintf(f,"PBM: %15.10f\n",pbm/n);		
//	if(count){
//		tg = nmi1/count;		
//	}else tg = -1;
//	fprintf(f,"NMI: %15.10f\n",tg);
//	fprintf(f,"IFV: %15.10f\n",ifv/n);
	// fprintf(f,"Ca: %15.10f\n",ca/n);
	// fprintf(f,"Ca_label: %15.10f\n",ca1/n);
	// fprintf(f,"Time: %15.10f\n",(runtime/n) + runti);
	// fprintf(f,"average_iteration: %15.10f\n",(double)averageIter/n);
	
	fclose(f);
	free(X);
	free(U);	
	free(V);
	free(cluster);
	return 0;
}

