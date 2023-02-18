#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PI 3.14159265358979323846

// FIR椭圆滤波器
void fir_elliptic_filter(double *input, double *output, int input_size, int order, double passband_frequency, double stopband_frequency, double passband_ripple, double stopband_attenuation) {
    int i, j;
    double *h;
    double *x;
    double *y;
    double omega_p;
    double omega_s;
    double delta_p;
    double delta_s;
    double beta;
    double gamma;
    double epsilon;
    int N;
    double M;
    double *cosine_array;
    double *sine_array;
    double omega_c;
    double A;
    double *b;
    double *a;
    double *w;
    double z;
    double y1, y2;
    double I=1;

    // 计算滤波器参数
    omega_p = 2 * PI * passband_frequency;
    omega_s = 2 * PI * stopband_frequency;
    delta_p = pow(10, passband_ripple / 20) - 1;
    delta_s = pow(10, -stopband_attenuation / 20);
    beta = asinh(sqrt(delta_s / delta_p) / 2) / order;
    gamma = asinh(sqrt(delta_s * delta_p)) / 2;
    epsilon = sinh(order * beta) * sinh(gamma);
    N = ceil(acosh(sqrt((1 / delta_p) - 1) / sqrt(delta_s)) / acosh((omega_s / omega_p)));
    M = ceil((N % 2 == 0 ? N + 1 : N) / 2.0);
    cosine_array = (double*) malloc(M * sizeof(double));
    sine_array = (double*) malloc(M * sizeof(double));
    h = (double*) malloc((order + 1) * sizeof(double));
    x = (double*) malloc((order + 1) * sizeof(double));
    y = (double*) malloc((order + 1) * sizeof(double));
    b = (double*) malloc((order + 1) * sizeof(double));
    a = (double*) malloc((order + 1) * sizeof(double));
    w = (double*) malloc((order + 1) * sizeof(double));
    omega_c = omega_p / cosh(acosh(sqrt((1 / delta_p) - 1) / sqrt(delta_s)) / N);
    A = pow(10, -passband_ripple / 20);
    for (i = 0; i < M; i++) {
        cosine_array[i] = cos((2 * i + 1) * PI / (2 * N));
        sine_array[i] = sin((2 * i + 1) * PI / (2 * N));
    }
    for (i = 0; i <= order; i++) {
        z = exp(2 * PI * i / (order + 1) * I);
        y1 = 1;
        y2 = 0;
        for (j = 0; j < M; j++){
            y2 = y1;
            y1 = z * y1 - cosine_array[j] * y2 - I * sine_array[j] * x[i];
            x[i] = z * x[i] - cosine_array[j] * y2 + I * sine_array[j] * y1;
    }
        h[i] = epsilon * (sinh(order * beta) * sin(omega_c * (i - order / 2)) / (PI * (i - order / 2)));
        h[i] = h[i] * (A + ((-1) ^ i) * epsilon);
        b[i] = h[i];
        a[i] = 0;
    }
    a[0] = 1;
    for (i = 0; i <= order; i++) {
    w[i] = 0;
}

    // 进行滤波
    for (i = 0; i < input_size; i++) {
    y[i] = b[0] * input[i] + w[0];
    for (j = 1; j <= order; j++) {
    y[i] += b[j] * input[i - j] - a[j] * y[i - j];
    }
    w[0] = b[order + 1] * input[i] - a[order + 1] * y[i];
    for (j = 1; j <= order; j++) {
    w[j] = b[order + 1 + j] * input[i] - a[order + 1 + j] * y[i] + w[j - 1];
    }
    output[i] = y[i];
    }

    // 释放内存
    free(cosine_array);
    free(sine_array);
    free(h);
    free(x);
    free(y);
    free(b);
    free(a);
    free(w);
    }

int main() {
    // 测试
    int input_size = 10;
    double input[] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
    double output[input_size];
    int order = 10;
    double passband_frequency = 0.2;
    double stopband_frequency = 0.3;
    double passband_ripple = 0.1;
    double stopband_attenuation = 60;
    fir_elliptic_filter(input, output, input_size, order, passband_frequency, stopband_frequency, passband_ripple, stopband_attenuation);

    for (int i = 0; i < input_size; i++) {
        assert(output[i]); 
        printf("%f ", output[i]);
    }
    printf("\n");
    return 0;
}
