#include"mylib.h"
int PERCENT_SELECT = 20;
int PERCENT_WRONG = 0;

int Bool[100];
int permu[100];
double caval = 0;

void xuat() {
	double max, tmp;
	int i,j,k,count = 0;
	
//    for(i = 0;i<C;i++){
//    	printf("%d ",permu[i]);	
//    }
//    printf("\n");
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

int *ramdomLabel(){
	int ita = 0,i,j,k,l,ij;
	int thres,tg,*idc,*markcla;
	idc = (int*)malloc(N*sizeof(int)); 
	 
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
		idc[i] = -1;	
	}
		
	// Random select label
	//printf("Begin random select label\n");
	j = 0;
	while(ita < NLA){
		k = rand() % CLASS[j];
		l = 0;
		for(i = 0;i<N;i++){
			if(cluster[i] == j && l < k && idc[i] == -1){
				for(ij = 0;ij < D;ij++){
					Xla[ita][ij] = X[i][ij];
				}
				CLASSLA[ita] = j;
				idc[i] = j;
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
	//printf("random label ok!\n");
	return idc;
}

double FCM(double prepareTime, char *fileout, int *iter){
	int i,j,k,l,step = 0;
	double isNext = 1, lambda = 1;
	double tg1 = 0,tgg,tg2 = 0, tg3, runtime = 0;
	double **vold;
	char s[50];
	clock_t t;
	t = clock();        
	// find minmax
	double *min, *max;
	//clust = ramdomLabel();
			
	//printf("Calculate min and max X\n");
	min = (double*)malloc(D*sizeof(double));
	max = (double*)malloc(D*sizeof(double));	
	vold = (double**)malloc(C*sizeof(double*));	
	for(j = 0;j<D;j++){
		min[j] = X[0][j];
		max[j] = X[0][j];
		for(i = 1;i<N;i++){
			if(min[j] > X[i][j]){
				min[j] = X[i][j];
			}
			if(max[j] < X[i][j]){
				max[j] = X[i][j];
			}
		}
	}
	
	printf("cal distance\n");
	
	tg1 = -1;
	for(i = 0;i<N - 1;i++){
		for(j = i+1;j<N;j++){
			tg2 = 0;
			for(l = 0;l < D;l++){
				tg2 += pow(X[i][l] - X[j][l],2);
			}
			//printf("%10.5lf\n",sqrt(tg2));
			if(tg1 == -1 || tg1 > tg2){
				tg1 = tg2;
			}
		}
	}	
	EPS = tg1/2;
	if(EPS < 0.001) EPS = 0.001;
			
	// randomV
	for(i = 0;i<C;i++){
		vold[i] = (double*)malloc(D*sizeof(double));  
		for(j = 0;j<D;j++){
			V[i][j] = min[j] + (max[j] - min[j]) * (double)(rand() % 100) / 100;
			vold[i][j] = V[i][j];			
		}
	}
	printf("begin fcm\n",step);
	// FCM for label
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
		printf("step:%ld\n",step);    
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
	*iter = step;    
	//printf("SFCM:   Finish with iter %d\n",step);	    
	outPutV(fileout, step, runtime,1); 
	printf("FCM:   OutputV ok with %ld steps\n",step);
		
	free(min);
	free(max);	
	free(vold);			
	free(CLASSLA);	
	return runtime;
}


int main(int  argc, char **argv){
	
	char s1[50],s[50];
	FILE *f;
	int i,n = 10,iter,averageIter = 0,count = 0;;
	double runtime = 0,runti,ifv = 0,db = 0,ma = 0,ca=0,nmi1 = 0,ti,ta,tg, runtime1;
	double db1,ri1,asw1,asw,pbm1,pbm;
	double ma1 = 0,ma2 = 0,ri = 0;
//	ALPHA = 0.5;
//	EPS = 0.01;
	
	srand(time(NULL));
	
	if (argc >= 3){ 
		sprintf(s1,"data/%s.txt",argv[1]);
		n = atoi(argv[2]);
		if(n < 1) n = 1; 
		//ALPHA = (double)atof(argv[3]);		
				
	}else{
		printf("invalid input parameter ifs\n");
		exit(0);
	}
			
	runti = input(s1);	
//	normalizeX();
//	printf("input ok!\n");
	sprintf(s1,"result/fcm/FCM_origin_result_%s.csv",argv[1]);
	f = fopen(s1,"w");
	
	if(!f){
		printf("Can't open file %s!",s1);		
		exit(0);
	}
	fprintf(f,"Result_file_%s:\n",argv[1]);
	fprintf(f,"Time,ca,db,asw,pbm,runtime\n");
	for(i = 1;i<=n;i++){		
		// printf("\nTime %d :\n",i);
		// sprintf(s,"result/1/fcm_origin_%s_time%d.txt",argv[1],i);
		runtime1 = FCM(runti,s,&iter);		
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
		tg = validity_kuhn();
		ca += tg;
		asw1 = ASWC(s);
		asw += asw1;
		pbm1 = PBM(s);
		pbm += pbm1;
		
		// fprintf(f,"%d,%10.5lf,%10.5lf,%10.5lf,%10.5lf,%10.5lf\n",i,tg,db1,asw1,pbm1,runtime1);
				
		// printf("ca = %10.5lf asw=%10.5lf pbm=%10.5lf runtime=%10.5lf\n",tg,asw1,pbm1,runtime1);
		// printf("DB = %10.5lf",db1);
		//printf("RI=%10.5lf",ri1);
		averageIter += iter;	
	}
		
	//fprintf(f,"\nRI: %15.10f\n",ri/n);
	// fprintf(f,"DB: %15.10f\n",db/n);
	// fprintf(f,"ASW: %15.10f\n",asw/n);	
	// fprintf(f,"PBM: %15.10f\n",pbm/n);			
//	if(count){
//		tg = nmi1/count;		
//	}else tg = -1;
//	fprintf(f,"NMI: %15.10f\n",tg);
//	fprintf(f,"IFV: %15.10f\n",ifv/n);
	// fprintf(f,"Ca: %15.10f\n",ca/n);
	// fprintf(f,"Time: %15.10f\n",(runtime/n) + runti);
	// fprintf(f,"average_iteration: %15.10f\n",(double)averageIter/n);
	
	fclose(f);
	free(X);
	free(U);	
	free(V);
	free(cluster);

	// Print clustering metrics (RI, F1-score, NMI, DB)
	printf("\n==== Clustering Metrics ====");
	printf("\nRI: %10.5lf", RI(s1));
	printf("\nF1-score: %10.5lf", F1_SCORE(s1));
	printf("\nNMI: %10.5lf", nmi(s1, NULL));
	printf("\nDB: %10.5lf\n", DB(s1));
	return 0;
}

