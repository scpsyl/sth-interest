#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

void fft(double complex *X, int N)
{
if (N == 1) return;

double complex *X1 = (double complex*) malloc(N/2 * sizeof(double complex));
double complex *X2 = (double complex*) malloc(N/2 * sizeof(double complex));

for (int i = 0; i < N/2; i++) {
X1[i] = X[2*i];
X2[i] = X[2*i+1];
}

fft(X1, N/2);
fft(X2, N/2);

for (int k = 0; k < N/2; k++) {
double complex t = X1[k];
double complex u = cexp(-2 * PI * I * k / N) * X2[k];
X[k] = t + u;
X[k + N/2] = t - u;
}

free(X1);
free(X2);
}

void dht(double *x, int N)
{
    double *y = (double*) malloc(N * sizeof(double));
    double complex *X = (double complex*) malloc(N * sizeof(double complex));

    for (int k = 0; k < N; k++) {
        double s = 0;
        for (int n = 0; n < N; n++) {
            s += x[n] * cos(2 * PI * k * n / N);
        }
        y[k] = s;
    }

    for (int n = 0; n < N; n++) {
        X[n] = y[(N-n) % N] + I * y[n];
    }

    fft(X, N);

    for (int n = 0; n < N; n++) {
        x[n] = creal(X[n]) / N;
    }

    free(y);
    free(X);
}

int main()
{
    int N = 8;
    double x[] = {1, 2, 3, 4, 4, 3, 2, 1};
    printf("%d\n", sizeof(x)/ sizeof(x[0]));
    printf("Input:\n");
    for (int n = 0; n < N; n++) {
        printf("%.2f ", x[n]);
    }
    printf("\n");

    dht(x, N);

    printf("Output:\n");
    for (int n = 0; n < N; n++) {
        printf("%.2f ", x[n]);
    }
    printf("\n");

    return 0;
}
