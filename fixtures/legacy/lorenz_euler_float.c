/* Lorenz attractor -- forward Euler, single precision (legacy C). */
#include <stdio.h>

float x, y, z;
float h = 0.01f;

int main(void)
{
    float sigma = 10.0f, rho = 28.0f, beta = 8.0f / 3.0f;
    float dx, dy, dz;
    int i;

    x = 1.0f; y = 1.0f; z = 1.0f;
    for (i = 0; i < 50000; i++) {
        dx = sigma * (y - x);
        dy = x * (rho - z) - y;
        dz = x * y - beta * z;
        x = x + h * dx;
        y = y + h * dy;
        z = z + h * dz;
    }
    printf("%f %f %f\n", x, y, z);
    return 0;
}
