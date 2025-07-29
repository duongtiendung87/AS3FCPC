// Backup of original as3fcpc.c
// Date: 2025-07-28

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define N 150     // number of data points
#define D 4       // number of features
#define C 3       // number of clusters
#define MAX_ITER 100
#define ALPHA 0.8 // boundary weight
#define NQ 10     // number of queried boundary points

// Limitation-related constants (e.g., neighborhood radius not used)
double convergenceThreshold = 1e-4; // configurable

// Data structures
double data[N][D];         // dataset
int label[N];              // expert-labeled data (partially filled)
double U[N][C];            // membership matrix
double V[C][D];            // cluster centers
int must_link[N][N];       // pairwise constraints
int cannot_link[N][N];
int queried[N];            // boundary query flag

// Function declarations
void loadData();
void initializeU();
void updateCenters();
void updateMembership();
void detectBoundaries();
void applyConstraints();
void adjustMembership();
int hasConverged(double oldV[C][D]);
void printResults();

double computeDik(int i, int k);
double computeUbar(int i);

int main() {
    int iter = 0;
    double oldV[C][D];

    loadData();
    initializeU();
    updateCenters();

    do {
        for (int i = 0; i < C; i++)
            for (int j = 0; j < D; j++)
                oldV[i][j] = V[i][j];

        updateMembership();
        detectBoundaries();
        applyConstraints();
        adjustMembership();
        updateCenters();
        iter++;
    } while (!hasConverged(oldV) && iter < MAX_ITER);

    printResults();
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
