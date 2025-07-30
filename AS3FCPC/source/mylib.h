
// Macro F1-score for clustering (refactored: pass all needed arguments)
double F1_SCORE(int *trueLabels, double **U, int N, int C) {
	int *predLabels = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++) {
		double maxU = U[i][0];
		int maxIdx = 0;
		for (int j = 1; j < C; j++) {
			if (U[i][j] > maxU) { maxU = U[i][j]; maxIdx = j; }
		}
		predLabels[i] = maxIdx;
	}
	double f1Sum = 0.0;
	for (int c = 0; c < C; c++) {
		int tp = 0, fp = 0, fn = 0;
		for (int i = 0; i < N; i++) {
			if (predLabels[i] == c && trueLabels[i] == c) tp++;
			if (predLabels[i] == c && trueLabels[i] != c) fp++;
			if (predLabels[i] != c && trueLabels[i] == c) fn++;
		}
		double precision = tp + fp > 0 ? (double)tp / (tp + fp) : 0.0;
		double recall = tp + fn > 0 ? (double)tp / (tp + fn) : 0.0;
		double f1 = precision + recall > 0 ? 2 * precision * recall / (precision + recall) : 0.0;
		f1Sum += f1;
	}
	free(predLabels);
	return f1Sum / C;
}

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<conio.h>

#include <ctype.h>


int read_mask_to_CLASS(const char *filename, int N) {
	FILE *f = fopen(filename, "r");
	if (!f) {
		printf("Can't open mask file %s!\n", filename);
		return -1;
	}
	char header[10];
	int is_pgm = 0, width = 0, height = 0, maxval = 0;
	// Đọc header PGM
	if (fgets(header, sizeof(header), f)) {
		if (strncmp(header, "P2", 2) == 0) is_pgm = 1;
	}
	if (CLASS) free(CLASS);
	CLASS = (int*)malloc(N * sizeof(int));
	if (is_pgm) {
		// Bỏ qua comment
		char line[256];
		do {
			if (!fgets(line, sizeof(line), f)) break;
		} while (line[0] == '#');
		// Đọc width height
		sscanf(line, "%d %d", &width, &height);
		fscanf(f, "%d", &maxval);
		int npix = width * height;
		if (npix != N) {
			printf("Mask size %d không khớp với N=%d!\n", npix, N);
			fclose(f);
			return -1;
		}
		for (int i = 0; i < N; i++) {
			int pix = 0;
			fscanf(f, "%d", &pix);
			// 0 hoặc 255 (hoặc maxval) -> nhãn 0 hoặc 1
			CLASS[i] = (pix >= maxval/2) ? 1 : 0;
		}
	} else {
		// Đọc txt: mỗi dòng là 1 pixel hoặc 1 dòng = 1 hàng ảnh
		rewind(f);
		int count = 0;
		char line[2048];
		while (fgets(line, sizeof(line), f) && count < N) {
			if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
			char *p = strtok(line, " ,\t\r\n");
			while (p && count < N) {
				int val = atoi(p);
				CLASS[count] = (val >= 128) ? 1 : 0;
				count++;
				p = strtok(NULL, " ,\t\r\n");
			}
		}
		if (count != N) {
			printf("Mask txt size %d không khớp với N=%d!\n", count, N);
			fclose(f);
			return -1;
		}
	}
	fclose(f);
	return 0;
}


int read_image_to_X(const char *filename, int *out_N, int *out_D) {
	FILE *f = fopen(filename, "r");
	if (!f) {
		printf("Can't open image file %s!\n", filename);
		return -1;
	}
	char header[10];
	int width = 0, height = 0, maxval = 0;
	int is_pgm = 0;
	// Đọc header PGM
	if (fgets(header, sizeof(header), f)) {
		if (strncmp(header, "P2", 2) == 0) is_pgm = 1;
	}
	if (is_pgm) {
		// Bỏ qua comment
		char line[256];
		do {
			if (!fgets(line, sizeof(line), f)) break;
		} while (line[0] == '#');
		// Đọc width height
		sscanf(line, "%d %d", &width, &height);
		fscanf(f, "%d", &maxval);
		*out_N = width * height;
		*out_D = 1;
		// Cấp phát X
		X = (double**)malloc((*out_N) * sizeof(double*));
		for (int i = 0; i < *out_N; i++) X[i] = (double*)malloc(sizeof(double));
		for (int i = 0; i < *out_N; i++) {
			int pix = 0;
			fscanf(f, "%d", &pix);
			X[i][0] = (double)pix / maxval;
		}
	} else {
		// Đọc txt: mỗi dòng là 1 pixel hoặc 1 dòng = 1 hàng ảnh
		rewind(f);
		int n = 0, d = 0;
		char line[2048];
		// Đếm số dòng và số cột
		while (fgets(line, sizeof(line), f)) {
			if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
			int cnt = 0;
			char *p = strtok(line, " ,\t\r\n");
			while (p) { cnt++; p = strtok(NULL, " ,\t\r\n"); }
			if (cnt > d) d = cnt;
			n++;
		}
		*out_N = n;
		*out_D = d;
		X = (double**)malloc(n * sizeof(double*));
		for (int i = 0; i < n; i++) X[i] = (double*)malloc(d * sizeof(double));
		rewind(f);
		int row = 0;
		while (fgets(line, sizeof(line), f)) {
			if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
			int col = 0;
			char *p = strtok(line, " ,\t\r\n");
			while (p && col < d) {
				X[row][col] = atof(p);
				col++;
				p = strtok(NULL, " ,\t\r\n");
			}
			// Nếu thiếu giá trị, điền 0
			while (col < d) { X[row][col] = 0.0; col++; }
			row++;
		}
	}
	fclose(f);
	return 0;
}

long N,NLA,MAXSTEPS = 1000;
int D,C;
double EPS = 0.001;
double ALPHA = 0.5;
double M,M1,M2,A,B,BETA;
double **U, **V, **T, **W, **X;
double **Ula, **Vla, **Xla;
double ASWC_EPS = 0.01;
double alpha = 1;
int *cluster,*CLASS,*CLASSLA;
int isInit = 0;
int *num_must_link;
int *num_cannot_link;
int **must_link;
int **cannot_link;
int *label;
int *queried;
double calcX_subtract_V(int k, int j){
	double tg = 0;
	int i;
	for(i = 0;i<D;i++){
		tg += pow(X[k][i] - V[j][i],2);                     
	}   
	return tg;
} 

double calcX_subtract_V2(double *x2, double *v2){
	double tg = 0;
	int i;
	for(i = 0;i<D;i++){
		tg += pow(x2[i] - v2[i],2);                     
	}   
	return tg;
} 

double input(char *filename){
	FILE *f;
	int i,j,k,type,nclass,ila,jla;
	double tg;
	char s[50],*p;
	clock_t t;
	t = clock();

	f = fopen(filename,"r");
	if(!f){
		  printf("Can't open file %s in output V!",filename);
	//              scanf("%d",&i);
		  exit(0);
	} 
	
	fscanf(f,"N=%ld\n",&N);
	fscanf(f,"R=%d\n",&D);
	fscanf(f,"C=%d\n",&C);
	fscanf(f,"m=%lf",&M);

	if (N + D + C < 3){
		printf("Incorrect data!");
		exit(0);
	}
	fscanf(f,"\nclass=%d\n",&nclass);	
	CLASS = (int*)malloc(nclass*sizeof(int));
	for(i =0;i<nclass;i++){
		fscanf(f,"%d ",&j);
		CLASS[i] = j;
	}
	
	fscanf(f,"\ntype=%d\n",&type);
	fscanf(f,"%s\n",s); 
	
	if(type == 0){		
		//if (cluster == NULL) 
		cluster = (int*)malloc(N*sizeof(int));		
		//if(X == NULL)  
		X = (double**)malloc(N*sizeof(double*));
								 
		for(i = 0;i<N;i++){
			//cluster[i] = 0;
			//if(X[i] == NULL) 
			X[i] = (double*)malloc(D*sizeof(double));			
			fscanf(f,"%d ",&k);			
			cluster[i] = k - 1;				
			for(j = 0;j<D;j++){
				fscanf(f,"%lf",&tg);
				X[i][j] = tg;
			}
		}
		
		U = (double**)malloc(N*sizeof(double*));
				
		for(i = 0;i<N;i++){
			U[i] = (double*)malloc(C*sizeof(double));			
		}
		
		Vla = (double**)malloc(C*sizeof(double*));
		V = (double**)malloc(C*sizeof(double*));
		for(j = 0;j<C;j++){
			Vla[j] = (double*)malloc(D*sizeof(double));				
			V[j] = (double*)malloc(D*sizeof(double));				
		}	
		
	}else{	
		if(U == NULL)  U = (double**)malloc(N*sizeof(double*));		
		fscanf(f,"%s\n",s); 
		for(i = 0;i<N;i++){
			if(U[i] == NULL)  U[i] = (double*)malloc(C*sizeof(double));			
			for(j = 0;j<C;j++){
				fscanf(f,"%lf",&tg);
				U[i][j] = tg;
			}
		}
		
		if(type > 1){
			if(T == NULL)  T = (double**)malloc(N*sizeof(double*));		
			fscanf(f,"%s\n",s); 
			for(i = 0;i<N;i++){
				if(T[i] == NULL)  T[i] = (double*)malloc(C*sizeof(double));			
				for(j = 0;j<C;j++){
					fscanf(f,"%lf",&tg);
					T[i][j] = tg;                   
				}
			}

			if(type == 3){
				if(W == NULL)  W = (double**)malloc(N*sizeof(double*));
				fscanf(f,"%s\n",s); 
				for(i = 0;i<N;i++){
					if(W[i] == NULL)  W[i] = (double*)malloc(C*sizeof(double));			
					for(j = 0;j<C;j++){
						fscanf(f,"%lf",&tg);
						W[i][j] = tg;                   
					}
				}
			}
		}		
		if(V == NULL)  V = (double**)malloc(C*sizeof(double*));
		fscanf(f,"%s\n",s);
		for(i = 0;i<C;i++){
			if(V[i] == NULL)  V[i] = (double*)malloc(C*sizeof(double));
			for(j = 0;j<D;j++){
				fscanf(f,"%lf",&tg);
				V[i][j] = tg;       				
			}
		}
	}	   
	fclose(f);
	isInit = 1;

//	printf("input is ok!\n");
	   
	return (double)(clock() - t)/CLOCKS_PER_SEC; 
}

double **maxminValueAtt(){
	double **tg, max, min;
	int i,j;
			   
	tg = (double**)malloc(D*sizeof(double*));
	for(i = 0;i<D;i++){
		tg[i] = (double*)malloc(2*sizeof(double));                
	}

	for(i = 0;i<D;i++){
		max = X[0][i];
		min = X[0][i];  
		for(j = 1;j<N;j++){
		  if(max < X[j][i]) max = X[j][i];
		  if(min > X[j][i]) min = X[j][i];
		}
		tg[i][0] = max;
		tg[i][1] = min;               
	}         

	return tg;
}

double randomV(int option){
	 int i,j,k,tg;
	 double **maxmin = maxminValueAtt();
	 clock_t t;
	 t = clock();
	 U = (double**)malloc(N*sizeof(double*));
	 if(option == 1) T = (double**)malloc(N*sizeof(double*));
	 for(i = 0;i<N;i++){
		U[i] = (double*)malloc(C*sizeof(double));
		if(option == 1) 
			T[i] = (double*)malloc(C*sizeof(double));      
	 }      
		  
	 srand(time(NULL));
	 if (isInit == 1){
		 V = (double**)malloc(C*sizeof(double*));
		 for(i = 0;i<C;i++){
			   V[i] = (double*)malloc(D*sizeof(double));                
			   for(j = 0;j<D;j++){
					  V[i][j] = ((double)(rand() % 10001) / 10000) * (maxmin[j][0] - maxmin[j][1]);
					  V[i][j] += maxmin[j][1];     
			   }
		 }
		 
	 }else{
		 printf("The data input is empty!");
		 exit(0);        
	 }
	 free(maxmin);
	 return (double)(clock() - t)/CLOCKS_PER_SEC;   
}

double randomUTW(int option){
	 int i,j,k,tg;
	 clock_t t;
	 t = clock();
//     srand(time(NULL));
	 if (isInit == 1){
		 if(option == 1){
//            if(U == NULL) 
				U = (double**)malloc(N*sizeof(double*));
//			if(T == NULL) 
				T = (double**)malloc(N*sizeof(double*));
			for(i = 0;i<N;i++){
//				if(U[i] == NULL) 
					U[i] = (double*)malloc(C*sizeof(double));
//				if(T[i] == NULL) 
					T[i] = (double*)malloc(C*sizeof(double));
				for(j = 0;j<C;j++){
					   U[i][j] = ((double)(rand() % 9999) + 1) / 10000;
					   T[i][j] = (1 - U[i][j]) * ((double)(rand() % 101) / 100);
				}            
			}
			isInit = 3;                        
		 }else if (option == 2){
//			if(U == NULL) 
				U = (double**)malloc(N*sizeof(double*));
//			if(T == NULL) 
				T = (double**)malloc(N*sizeof(double*));
//			if(W == NULL) 
				W = (double**)malloc(N*sizeof(double*));
			for(i = 0;i<N;i++){
//				if(U[i] == NULL) 
					U[i] = (double*)malloc(C*sizeof(double));
//				if(T[i] == NULL) 
					T[i] = (double*)malloc(C*sizeof(double));
//				if(W[i] == NULL) 
					W[i] = (double*)malloc(C*sizeof(double));
				for(j = 0;j<C;j++){
					   U[i][j] = ((double)(rand() % 9999) + 1) / 10000;
					   T[i][j] = (1 - U[i][j]) * (((double)(rand() % 9999) + 1) / 10000);
					   W[i][j] = (1 - U[i][j] - T[i][j]) * (((double)(rand() % 9999) + 1) / 10000);
				}            
			}
			isInit = 4;     
		 }else{
			U = (double**)malloc(N*sizeof(double*));
			for(i = 0;i<N;i++){
				U[i] = (double*)malloc(C*sizeof(double));
				for(j = 0;j<C;j++){
					U[i][j] = (double)(rand() % 10001) / 10000;
			//                           printf("%10.5f ",U[i][j]);scanf("%d",&k);
				}            
			}
			   isInit = 2;                   
		 } 
		 
//      if(V == NULL) 
			V = (double**)malloc(C*sizeof(double*));
		for(i = 0;i<C;i++){
//          if(V[i] == NULL) 
				V[i] = (double*)malloc(D*sizeof(double));
			for(j = 0;j<D;j++){
				V[i][j] = 0;
			}
		}
	}else if(isInit == 0){
		 printf("Incorect data!");        
		 exit(0);
	}
	 return (double)(clock() - t)/CLOCKS_PER_SEC; 
}

void normalizeX(){
	int i,j,k;
//	printf("id%d begin normalize\n",id);
	double **tg = maxminValueAtt();	
//	printf("id%d begin normalize(2)\n",id);
	for(i = 0;i<N;i++){
		for(j = 0;j<D;j++){
			if(tg[j][0] - tg[j][1] != 0){
				X[i][j] = (X[i][j] - tg[j][1]) / (tg[j][0] - tg[j][1]);
				if(X[i][j] < 0) X[i][j] = 0;	
			}
		}	
	}	
	free(tg);		
}

// Output of U and V
void outPutV(char *filename, int iter, double running_time, int type){
	 int i,j,k;
	 FILE *f;
	 f = fopen(filename,"w");
	 if(!f){
		  printf("Can't open file %s!",filename);
		  //scanf("%d",&i);
		  exit(0);
	 }
	fprintf(f,"N=%d\n",N);
	fprintf(f,"R=%d\n",D);
	fprintf(f,"C=%d\n",C);
	fprintf(f,"m=%lf\n",M);
	fprintf(f,"class=%d\n",C);
	for(i = 0;i<C;i++){
		fprintf(f,"%d ",CLASS[i]);
	}
	fprintf(f,"\ntype=%d\n",type);
	fprintf(f,"U:\n");
	 
	for(i = 0;i<N;i++){
		for(j = 0;j<C;j++){
			fprintf(f,"%7.6lf ",U[i][j]);          
		}      
		fprintf(f,"\n");
//		printf("i = %d\n",i);
	}	
	if(type > 1){
		fprintf(f,"T:\n");	
//		printf("T:\n");
		for(i = 0;i<N;i++){
			for(j = 0;j<C;j++){
				fprintf(f,"%7.6lf ",T[i][j]);          
			}      
			fprintf(f,"\n");
			//printf("i = %d\n",i);
		}
		if(type == 3){
			fprintf(f,"W:\n");		 
//			printf("W:\n");
			for(i = 0;i<N;i++){
				for(j = 0;j<C;j++){
					fprintf(f,"%7.6lf ",W[i][j]);          
				}      
				fprintf(f,"\n");
				//printf("i = %d\n",i);
			}
		}
	} 
	fprintf(f,"V:\n");

	for(i = 0;i<C;i++){
		for(j = 0;j<D;j++){
			fprintf(f,"%7.6lf ",V[i][j]);          
		}      
		fprintf(f,"\n");
	}

	fprintf(f,"Leap= %d\n",iter);
	fprintf(f,"Time= %10.5f\n",running_time);
	fprintf(f,"Max iterations= %d\n",MAXSTEPS);          
	fclose(f);          
}

int *repmat(int *a, int n, int m, int n1, int m1){
	int *b = (int*)malloc(n*n1*m*m1*sizeof(int));
	int i,j;
	for(i = 0;i<n*n1;i++){		
		for(j = 0;j<m*m1;j++){
			b[i * m * m1 + j] = a[(i%n) * m + (j%m)];
		}
	}	
	return b;
}

double *checkMatrix(int *a, int *b, int n, int m){
	int i,j;	
	for(i = 0;i<n;i++){
		for(j = 0;j<m;j++){
			if(a[i*m + j] == b[i*m + j])
				a[i*m + j] = 1;
			else
				a[i*m + j] = 0;
		}
	}
	
	double *sumMatrix = (double*)malloc(m*sizeof(double));
	
	for(j = 0;j<m;j++){
		sumMatrix[j] = 0;
		for(i = 0;i<n;i++){		
			sumMatrix[j] += a[i*m + j];
		}
	}
	return sumMatrix;
}

double nmi(char *filename, double *time){
	double result,max,tg,hl,hr,hlr,mi,v;
	*time = input(filename);
	//printf("time = %10.5lf\n",*time);
	int i,j,k;
		
	int *h,*g,*h_u,*g_u;	
	double *pl,*pr,*m;
	int *res_u = (int*)malloc(C*sizeof(int));
	
	int *res = (int*)malloc(N * sizeof(int));
	int *chek = (int*)malloc(C * sizeof(int));
	
	for(j = 0;j<C;j++){
		chek[j] = 0;
		res_u[j] = j;
	}
	
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			tg = U[i][j];
			if(max < tg){
				max = tg;
				k = j;
			}
		}		
		res[i] = k;
		chek[k] = 1;
	}
	
	k = 0;
	for(j = 0;j<C;j++){
		k += chek[j];
	}	
	free(chek);
	if(k < C){
		printf("incorrect data in %s, k = %d\n",filename,k);
		free(res_u);
		free(res);		
		return -1000;
	}else{		
		
		h = repmat(cluster, N, 1, 1, C);
		h_u = repmat(res_u,1,C,N,1);
		
		g = repmat(res, N, 1, 1, C);
		g_u = repmat(res_u,1,C,N,1);
		
		pl = checkMatrix(h,h_u,N,C);
		pr = checkMatrix(g,g_u,N,C);
		
		free(h_u);
		free(g_u);
		free(res);
		free(res_u);
			
		hl = 0;
		hr = 0;
		for(i = 0;i<C;i++){
			hl += -pl[i] * log(pl[i] + EPS/100)/log(2);
			hr += -pr[i] * log(pr[i] + EPS/100)/log(2);;
		}
		
		m = (double*)malloc(C*C*sizeof(double));
			
		for(j = 0;j<C;j++){
			for(k = 0;k<C;k++){
				m[j*C + k] = 0;
				for(i = 0;i<N;i++){		
					m[j*C + k] += h[i * C + j] * g[i * C + k];
				}
				m[j*C + k] = m[j*C + k]/N;
//				printf("%5.3lf ",m[j*C + k]);
			}
//			printf("\n");
		}
		hlr = 0;		
		for(j = 0;j<C;j++){
			for(k = 0;k<C;k++){
				hlr += -m[j*C + k] * (log(m[j*C + k] + EPS/100)/log(2));
//				printf("%5.3lf[%5.3lf-%5.3lf] ",hlr,m[j*C + k],log(m[j*C + k] + EPS)/log(2));
			}
		}	
		
//		printf("\nhlr = %10.5lf\n",hlr);
		
		mi = hl + hr - hlr;
		free(h);
		free(g);	
		free(m);
		free(pl);
		free(pr);	
		return sqrt((mi/hl) * (mi/hr));
	}
}

double MA(char *filename){
	double time,max;
	int i,j,k,count = 0,flag = 0,l,sub;	
	time = input(filename);
		
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			if(max < U[i][j]){ 
				max = U[i][j];
				k = j;
			}
		}	
		if(k == cluster[i]) count++;
	}
	return (count * 100) / N;
}

double MA1(char *filename){
	double time,max;
	int i,j,k,count = 0,flag = 0,l,sub;	
	time = input(filename);

	int *mark = (int*)malloc(N*sizeof(int));
	int *clas = (int*)malloc(C*sizeof(int));
		
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			if(max < U[i][j]){ 
				max = U[i][j];
				k = j;
			}
		}	
		mark[i] = k;		
	}
	
	for(j = 0;j<C;j++){
		clas[j] = 0;
		for(i = 0;i<N;i++){
			if(mark[i] == j) clas[j]++;
		}
	}
	
	max = 0;
	for(j = 0;j<C;j++){
		sub = abs(clas[j] - CLASS[j]);
		if(max < sub) max = sub;
	}
	free(mark);
	free(clas);
	return (N - max) * 100 / N;
}

double RI(char *filename){
	double time,max;
	int i,j,k,a = 0,b = 0,c = 0,d = 0;	
	time = input(filename);

	int *mark = (int*)malloc(N*sizeof(int));
	
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			if(max < U[i][j]){ 
				max = U[i][j];
				k = j;
			}
		}	
		mark[i] = k;				
	}
	
	for(i = 0;i<N-1;i++){
		for(j = i + 1;j<N;j++){
			if(cluster[i] == cluster[j] && mark[i] == mark[j]){
				a++;
			}else if(cluster[i] != cluster[j] && mark[i] != mark[j]){
				b++;			
			}else if(cluster[i] != cluster[j] && mark[i] == mark[j]){
				c++;
			}else if(cluster[i] == cluster[j] && mark[i] != mark[j]){
				d++;
			}
		}
	}
	
	free(mark);
	if(a + b + c + d == 0){
		printf("Loi rand index!\n");
		return -1;
	}else{
		return (double)(a + b) / (a + b + c + d);	
	}	
}

double calcDBS(int i,int *mark){
	int k = 0,j;
	double tg = 0;
	for(j = 0;j<N;j++){
		if(mark[j] == i){
			tg += calcX_subtract_V(j, i);		
			k++;
		}
	}
	if(k) return sqrt(tg/k);
	else {		
		printf("error division by zero at calcDBS!\n");
		//getch();
		return 0;
//		exit(0);		
	}
}

double calcDBM(int i,int j){
	int k;
	double tg = 0;
	for(k = 0;k<D;k++){
		tg += pow(V[i][k] - V[j][k],2);		
	}	
	return sqrt(tg);
}

double calcSDMax(){
	int i,j,k;
	double max = 0,tg;
	
	for(i = 0;i<C - 1;i++){
		for(j = i + 1;j<C;j++){
			tg = 0;
			for(k = 0;k<D;k++){
				tg += pow(V[i][k] - V[j][k],2);		
			}
			if(max < tg) max = tg;
		}
	}
	return max;
}

double sigmaD(){
	int i,j,k;
	double tg1 = 0,tg;
	for(j = 0;j<C;j++){
		tg = 0;
		for(i = 0;i<N;i++){
			tg += calcX_subtract_V(i, j);
		}
		tg1 += tg/N;
	}
	return tg1/C;
}

double DB(char *filename){
	double time,max,tg,tg1,tg2,a,tg3;
	int i,j,k,l,*mark;
	
	time = input(filename);
	mark = (int*)malloc(N*sizeof(int));
	
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			a = U[i][j];
			if(max < a){
				max = a;
				k = j;
			}		
		}
		mark[i] = k;
	}
	tg = 0;
	for(i = 0;i<C;i++){		
		max = 0;
		tg1 = calcDBS(i,mark);
		for(j = 0;j<C;j++){			
			if(i != j){
				tg3 = calcDBM(i,j);
				if(tg3 < 0.01) tg3 += 0.01;
				tg2 = (tg1 + calcDBS(j,mark)) / tg3;
				if(max < tg2) max = tg2;
			}
		}
		tg += max;
	}
	free(mark);
	if(tg/C > 100) return 0;
	else return tg / C;
}

double IFV(char *filename){
	double time,max,tg,tg1,tg2,sum;
	int i,j,k,l,*size;
		
	time = input(filename);
	sum = 0;
	for(j = 0;j<C;j++){	
		tg1 = 0;
		for(i = 0;i<N;i++){
			if(U[i][j] == 0) U[i][j] = EPS;
			if(U[i][j] == 1) U[i][j] = 1 - EPS;
			tg1 += log(U[i][j]) / log(2);
		}
		tg = pow(log(C)/log(2) - tg1/N,2);
		
		tg2 = 0;
		for(i = 0;i<N;i++){		
			tg2 += (double)pow(U[i][j],2);
		}
		tg2 = tg2/N;
		sum += tg * tg2;
//		printf("j = %d tg = %15.10lf sum = %15.10lf\n",j,tg,sum);	
	}
	
	return (sum * calcSDMax()) / (sigmaD() * C);	
}

void outPutSFCM(char *filename, double **v, double **u, int n, int iter){
	 int i,j,k;
	 FILE *f;
	 f = fopen(filename,"w");
	 if(!f){
		  printf("Can't open file %s!",filename);
		  scanf("%d",&i);
		  exit(0);
	 }
	fprintf(f,"N=%d\n",n);
	fprintf(f,"R=%d\n",D);
	fprintf(f,"C=%d\n",C);
	fprintf(f,"U:\n");     
	for(i = 0;i<n;i++){
		for(j = 0;j<C;j++){
			fprintf(f,"%7.6lf ",u[i][j]);          
		}      
		fprintf(f,"\n");
//		printf("i = %d\n",i);
	}	
	fprintf(f,"V:\n");
	for(i = 0;i<C;i++){
		for(j = 0;j<D;j++){
			fprintf(f,"%7.6lf ",v[i][j]);          
		}      
		fprintf(f,"\n");
	}
	fprintf(f,"Leap= %d\n",iter);
	fclose(f);          
}

double ASWC(char *filename){
	double time,max,tg,tg1,tg2,a,tg3;
	int i,j,k,l,*mark;	
	double sum, count = N;
	int markcount;
	
	int *countDif = (int*)malloc(C*sizeof(int));
	double *averageDif = (double*)malloc(C*sizeof(double));
	
	time = input(filename);
	
	mark = (int*)malloc(N*sizeof(int));
	
	for(i = 0;i<N;i++){
		max = U[i][0];
		k = 0;
		for(j = 1;j<C;j++){
			a = U[i][j];
			if(max < a){
				max = a;
				k = j;
			}		
		}
		mark[i] = k;
	}
		
	sum = 0;
//	printf("Begin ASWC C=%d\n",p->Cx);	
	for(i = 0;i<N;i++){
		markcount = 0;		
		for(j = 0;j<C;j++){
			averageDif[j] = 0;
			countDif[j] = 0;
		}
		for(j = 0;j<N;j++){
			if(i != j){
				tg = calcX_subtract_V2(X[i],X[j]);
				if(tg > 0) averageDif[mark[j]] += sqrt(tg);
				countDif[mark[j]]++;
			}
		}		
		max = -1;		
		for(j = 0;j<C;j++){
			if(countDif[j] > 0) averageDif[j] /= countDif[j];
			else averageDif[j] = 1;		
			if(j != mark[i]){
				if(max == -1 || max > averageDif[j]){
					max = averageDif[j];
				}			
			}
		}

		tg = max / (averageDif[mark[i]] + ASWC_EPS);
		sum += tg;
	}	
	
	free(averageDif);
	free(countDif);
	free(mark);
	return sum/N;	
}

double PBM(char *filename){
	int i,j,k;	
	double e1 = 0,ek = 0,dk = 0,tg;
	double time = input(filename);	
	
	double *value = (double*)malloc(D*sizeof(double));
	for(j = 0;j<D;j++){
		tg = 0;
		for(i = 0;i<N;i++){
			tg += X[i][j];	
		}
		value[j] = tg / N;	
	}
			
	for(i = 0;i<N;i++){
		tg = 0;
		for(j=0;j<D;j++){
			tg += pow(X[i][j] - value[j],2);
		}
		e1 += sqrt(tg);
	}

	for(i = 0;i<C;i++){
		for(j = 0;j<N;j++){
			ek += pow(U[j][i],M) * sqrt(calcX_subtract_V(j,i));			
		}		
	}
	
	for(i = 0;i<C - 1;i++){
		for(j = i + 1;j<C;j++){
			tg = 0;
			for(k = 0;k<D;k++){
				tg += pow(V[j][k] - V[i][k],2);
			}
			tg = sqrt(tg);
			if(dk < tg) dk = tg;
		}
	}
	free(value);
	if(ek == 0){
		printf("Loi PBM!\n");
		return -1;
		//exit(0);
	}
	return (e1 * dk)/(C * ek); 
}


