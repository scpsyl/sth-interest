#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

void fft(double* x, double* y, int n, int step) {
    if (n == 1) {
        return;
    }
    fft(x, y, n/2, 2*step);
    fft(x+step, y+step, n/2, 2*step);
    for (int i = 0; i < n/2; i++) {
        double t1 = x[i*2*step + step];
        double t2 = y[i*2*step + step];
        double t3 = cos(2*PI*i/n);
        double t4 = sin(2*PI*i/n);
        double t5 = x[i*2*step] - t1;
        double t6 = y[i*2*step] - t2;
        x[i*2*step] += t1;
        y[i*2*step] += t2;
        x[i*2*step + step] = t5*t3 - t6*t4;
        y[i*2*step + step] = t5*t4 + t6*t3;
    }
}

void DIF_FFT(double* x, double* y, int n) {
    fft(x, y, n, 1);
}

int main() {
    double x[8] = {1, 2, 3, 4, 4, 3, 2, 1};
    double y[8] = {0};
    int n = 8;

    DIF_FFT(x, y, n);

    for (int i = 0; i < n; i++) {
        printf("%f + %fi\n", x[i], y[i]);
    }

    return 0;
}
