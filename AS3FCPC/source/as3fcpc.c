#include "mylib.h"

// Declare missing global variables
int must_link[N][N] = {0};
int cannot_link[N][N] = {0};
int queried[N] = {0};
double ALPHA = 0.5;
double convergenceThreshold = 1e-4;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define N 150
#define D 4
#define C 3
#define MAX_ITER 100
#define PERCENT_SELECT 20
#define PERCENT_WRONG 30

double data[N][D];
int label[N];
double U[N][C];
double V[C][D];

int NLA;
double **Xla;
double **Ula;
int *CLASSLA;
int *clust;

// ...existing code...

void randomLabel() {
    clust = (int*)malloc(N * sizeof(int));
    NLA = (int)(PERCENT_SELECT * N / 100);
    Xla = (double**)malloc(NLA * sizeof(double*));
    Ula = (double**)malloc(NLA * sizeof(double*));
    CLASSLA = (int*)malloc(NLA * sizeof(int));
    for (int i = 0; i < NLA; i++) {
        Xla[i] = (double*)malloc(D * sizeof(double));
        Ula[i] = (double*)malloc(C * sizeof(double));
    }
    for (int i = 0; i < N; i++) clust[i] = -1;
    int ita = 0, j = 0;
    while (ita < NLA) {
        int k = rand() % (N / C);
        int l = 0;
        for (int i = 0; i < N; i++) {
            if (label[i] == j && l < k && clust[i] == -1) {
                for (int ij = 0; ij < D; ij++) Xla[ita][ij] = data[i][ij];
                CLASSLA[ita] = j;
                clust[i] = j;
                ita++; l++; break;
            }
        }
        j++; if (j == C) j = 0;
    }
    // Make wrong label
    if (PERCENT_WRONG > 0) {
        int tg = (int)(PERCENT_WRONG * NLA / 100);
        int *markcla = (int*)malloc(NLA * sizeof(int));
        for (int i = 0; i < NLA; i++) markcla[i] = 0;
        int ij = 0;
        while (ij < tg) {
            int k = rand() % tg;
            int j = 0;
            for (int i = 0; i < NLA; i++) {
                if (markcla[i] == 0) {
                    if (j == k) {
                        if (C == 2) CLASSLA[i] = (CLASSLA[i] + 1) % 2;
                        else CLASSLA[i] = ((rand() % (C - 1) + 1) + CLASSLA[i]) % C;
                        markcla[i] = 1;
                    } else j++;
                }
            }
            ij++;
        }
        free(markcla);
    }
}

double calcX_subtract_V2(double *x, double *v) {
    double sum = 0.0;
    for (int i = 0; i < D; i++) sum += pow(x[i] - v[i], 2);
    return sqrt(sum);
}

void FCM_label() {
    double **Vla = (double**)malloc(C * sizeof(double*));
    for (int i = 0; i < C; i++) Vla[i] = (double*)malloc(D * sizeof(double));
    // random Vla
    for (int i = 0; i < C; i++)
        for (int j = 0; j < D; j++)
            Vla[i][j] = Xla[0][j] + (Xla[NLA-1][j] - Xla[0][j]) * (double)(rand() % 100) / 100;
    int step = 0;
    double isNext = 1;
    double vold[C][D];
    while (isNext != -1 && step < MAX_ITER) {
        // update Ula
        for (int j = 0; j < C; j++) {
            for (int i = 0; i < NLA; i++) {
                double tg2 = 0;
                for (int l = 0; l < C; l++)
                    tg2 += pow(calcX_subtract_V2(Xla[i], Vla[j]) / calcX_subtract_V2(Xla[i], Vla[l]), 2);
                Ula[i][j] = 1.0 / tg2;
            }
        }
        // update Vla
        for (int j = 0; j < C; j++) {
            for (int d = 0; d < D; d++) {
                double num = 0, denom = 0;
                for (int i = 0; i < NLA; i++) {
                    num += pow(Ula[i][j], 2) * Xla[i][d];
                    denom += pow(Ula[i][j], 2);
                }
                Vla[j][d] = num / denom;
            }
        }
        // check convergence
        isNext = -1;
        for (int j = 0; j < C; j++) {
            double diff = 0;
            for (int d = 0; d < D; d++) diff += pow(Vla[j][d] - vold[j][d], 2);
            if (diff > 1e-4) isNext = diff;
            for (int d = 0; d < D; d++) vold[j][d] = Vla[j][d];
        }
        step++;
    }
    // transfer Ula to U for labeled points
    int l = 0;
    for (int i = 0; i < N; i++) {
        if (clust[i] >= 0) {
            for (int j = 0; j < C; j++) U[i][j] = Ula[l][j];
            l++;
        }
    }
    for (int i = 0; i < C; i++) free(Vla[i]);
    free(Vla);
}

void FCM_all() {
    double vold[C][D];
    for (int i = 0; i < C; i++)
        for (int d = 0; d < D; d++) vold[i][d] = V[i][d];
    int step = 0;
    double isNext = 1;
    while (isNext != -1 && step < MAX_ITER) {
        // update U
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < C; j++) {
                double tg2 = 0;
                for (int l = 0; l < C; l++)
                    tg2 += pow(calcX_subtract_V2(data[i], V[j]) / calcX_subtract_V2(data[i], V[l]), 2);
                U[i][j] = 1.0 / tg2;
            }
        }
        // update V
        for (int j = 0; j < C; j++) {
            for (int d = 0; d < D; d++) {
                double num = 0, denom = 0;
                for (int i = 0; i < N; i++) {
                    num += pow(U[i][j], 2) * data[i][d];
                    denom += pow(U[i][j], 2);
                }
                V[j][d] = num / denom;
            }
        }
        // check convergence
        isNext = -1;
        for (int j = 0; j < C; j++) {
            double diff = 0;
            for (int d = 0; d < D; d++) diff += pow(V[j][d] - vold[j][d], 2);
            if (diff > 1e-4) isNext = diff;
            for (int d = 0; d < D; d++) vold[j][d] = V[j][d];
        }
        step++;
    }
}


int main() {
    srand(time(NULL));
    loadData();
    randomLabel();
    FCM_label();
    // initialize V for all data
    for (int j = 0; j < C; j++)
        for (int d = 0; d < D; d++) V[j][d] = Xla[0][d] + (Xla[NLA-1][d] - Xla[0][d]) * (double)(rand() % 100) / 100;
    FCM_all();
    printResults();

    // Output metrics
    char resultFile[100];
    sprintf(resultFile, "result/as3fcpc_result.csv");
    FILE *f = fopen(resultFile, "w");
    if (!f) {
        printf("Can't open result file!\n");
        return 1;
    }
    double ri = RI(resultFile);
    double db = DB(resultFile);
    double nmi_time = 0.0;
    double nmi = nmi(resultFile, &nmi_time);
    // Use F1-score from mylib.h
    double f1 = F1_SCORE();
    fprintf(f, "RI,NMI,F1,DB\n");
    fprintf(f, "%10.5lf,%10.5lf,%10.5lf,%10.5lf\n", ri, nmi, f1, db);
    fclose(f);
    printf("RI = %10.5lf\n", ri);
    printf("NMI = %10.5lf\n", nmi);
    printf("F1 = %10.5lf\n", f1);
    printf("DB = %10.5lf\n", db);

    // free memory
    for (int i = 0; i < NLA; i++) { free(Xla[i]); free(Ula[i]); }
    free(Xla); free(Ula); free(CLASSLA); free(clust);
    return 0;
}

void loadData() {
    FILE *fp = fopen("data.csv", "r");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++)
            fscanf(fp, "%lf,", &data[i][j]);
        label[i] = -1;
    }
    fclose(fp);
}

void initializeU() {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < C; j++) {
            U[i][j] = rand() % 100 + 1;
            sum += U[i][j];
        }
        for (int j = 0; j < C; j++)
            U[i][j] /= sum;
    }
}

double computeUbar(int i) {
    double sum = 0.0;
    for (int j = 0; j < C; j++)
        sum += U[i][j];
    return sum / C;
}

double computeDik(int i, int k) {
    double dist = 0.0;
    for (int d = 0; d < D; d++)
        dist += pow(data[i][d] - V[k][d], 2);
    return dist;
}

void updateCenters() {
    for (int k = 0; k < C; k++) {
        double numerator[D] = {0};
        double denominator = 0;
        for (int i = 0; i < N; i++) {
            double uik = U[i][k];
            double ubar = computeUbar(i);
            double coeff = pow(uik, 2) + pow(uik - ubar, 2);
            for (int d = 0; d < D; d++)
                numerator[d] += coeff * data[i][d];
            denominator += coeff;
        }
        for (int d = 0; d < D; d++)
            V[k][d] = numerator[d] / denominator;
    }
}

void updateMembership() {
    for (int i = 0; i < N; i++) {
        double ubar = computeUbar(i);
        for (int j = 0; j < C; j++) {
            double dik = computeDik(i, j);
            double num = 0, denom = 0;

            for (int l = 0; l < C; l++)
                denom += pow(U[i][l], 2) + pow(U[i][l] - ubar, 2);

            num = (2 * pow(U[i][j], 2) + 2 * pow(U[i][j] - ubar, 2)) * dik;

            double constraintTerm = 0.0;
            for (int x = 0; x < N; x++) {
                for (int l = 0; l < C; l++) {
                    if (must_link[i][x] && l != j) constraintTerm += U[i][l] * U[x][l];
                    if (cannot_link[i][x] && l == j) constraintTerm += U[i][l] * U[x][l];
                }
            }

            U[i][j] = num / (4 * dik * denom + ALPHA * constraintTerm + 1e-10);
        }
    }
}

void detectBoundaries() {
    for (int i = 0; i < N; i++) {
        double max1 = 0, max2 = 0;
        for (int j = 0; j < C; j++) {
            if (U[i][j] > max1) {
                max2 = max1;
                max1 = U[i][j];
            } else if (U[i][j] > max2) {
                max2 = U[i][j];
            }
        }
        double diff = max1 - max2;
        queried[i] = (diff < 0.2);
    }
}

void applyConstraints() {
    for (int i = 0; i < N; i++) {
        if (queried[i]) {
            label[i] = rand() % C;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (label[i] != -1 && label[j] != -1) {
                if (label[i] == label[j]) must_link[i][j] = 1;
                else cannot_link[i][j] = 1;
            }
        }
    }
}

void adjustMembership() {
    for (int i = 0; i < N; i++) {
        if (queried[i]) {
            double max1 = 0, max2 = 0;
            int idx1 = -1, idx2 = -1;
            for (int j = 0; j < C; j++) {
                if (U[i][j] > max1) {
                    max2 = max1; idx2 = idx1;
                    max1 = U[i][j]; idx1 = j;
                } else if (U[i][j] > max2) {
                    max2 = U[i][j]; idx2 = j;
                }
            }
            double epsilon = 0.01;
            if ((max1 - max2) > epsilon) {
                double mid = (U[i][idx1] + U[i][idx2]) / 2.0;
                U[i][idx1] = mid - epsilon / 2.0;
                U[i][idx2] = mid + epsilon / 2.0;
            }
        }
    }
}

int hasConverged(double oldV[C][D]) {
    for (int j = 0; j < C; j++) {
        for (int d = 0; d < D; d++) {
            if (fabs(V[j][d] - oldV[j][d]) > convergenceThreshold)
                return 0;
        }
    }
    return 1;
}

void printResults() {
    for (int i = 0; i < N; i++) {
        printf("Data %d: ", i);
        for (int j = 0; j < C; j++) {
            printf("%.2f ", U[i][j]);
        }
        printf("\n");
    }
}
