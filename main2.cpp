#include <iostream>
#include <cmath>
#include <numeric>

const double eps = 1e-9;
const double pi = 3.141592653;

// f(x) = sinx/x
double f1(double x) {
    if (abs(x) < eps)
        return 1.0;
    return sin(x) / x;
}

double f2(double x) {
    return sin(x) * sin(x);
}

// 复化梯形公式
double trapezoid_int(double (*func)(double x), double a, double b, int n) {
    double h = (b - a) / n;
    double res = (func(a) + func(b));
    for (int i = 1; i <= n - 1; i++)
        res += 2 * func(a + i * h);
    return (h / 2) * res;
} 

// 复化simpson公式
double simpson_int(double (*func)(double x), double a, double b, int n) {
    double h = (b - a) / n;
    double res = (func(a) + func(b));
    for (int i = 1; i <= n - 1; i++)
        res += 2 * func(a + i * h);
    for (int i = 1; i <= n; i++)
        res += 4 * func(a + (i - 0.5) * h);
    return (h / 6) * res;
} 

int main()
{
    printf("trapezoid: %.10lf\n", trapezoid_int(f1, 0, 1, 100));
    printf("simpson: %.10lf\n", simpson_int(f1, 0, 1, 100));

    printf("trapezoid: %.10lf\n", trapezoid_int(f2, 0, pi / 2, 100));
    printf("simpson: %.10lf\n", simpson_int(f2, 0, pi / 2, 100));
    return 0;
}