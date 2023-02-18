#include <stdio.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

// DIT-FFT 算法实现
void fft(double *x,  double  *y, int n) {
    // 如果 n = 1，结束递归
    if (n == 1) {
        y[0] = x[0];
        return;
    }

    // 分治递归
    double even[n/2], odd[n/2], even_fft[n/2], odd_fft[n/2];
    for (int i = 0; i < n/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i+1];
    }
    fft(even, even_fft, n/2);
    fft(odd, odd_fft, n/2);

    // 合并结果
    for (int k = 0; k < n/2; k++) {
        double complex t = cexp(-I * 2 * PI * k / n) * odd_fft[k];
        y[k] = even_fft[k] + t;
        y[k + n/2] = even_fft[k] - t;
    }
}

int main() {
    // 输入信号
    int n = 4;
    double x[] = {2,1,3,4};
    double  y[n];

    // 计算 DIT-FFT
    fft(x, y, n);

    // 输出结果
    for (int i = 0; i < n; i++) {
        printf("%f + %f i\n", creal(y[i]), cimag(y[i]));
    }

    return 0;
}
