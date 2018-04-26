#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RK_order 4
#define KAPPA1 1
#define KAPPA2 2
#define ALPHA 2
#define BETA 1
#define GAMMA1 0.25
#define GAMMA2 0.25
#define SIGMA1 0.25
#define SIGMA2 0.25


/* define prototypes */

void virus_dynamic(float *, float *);
void rk4sys(int, float, float *, int);

void virus_dynamic(float *X, float *f) {
  f[0] = KAPPA1 * X[0] * X[2] * (1 - X[0] - X[1] - GAMMA1);
  f[1] = KAPPA2 * X[1] * X[2] * (1 - X[0] - X[1] - GAMMA2);
  f[2] = ALPHA * X[0] * (1 - X[2] - X[3]) - SIGMA1 * X[2];
  f[3] = BETA * X[1] * (1 - X[2] - X[3]) - SIGMA2 * X[3];
}

void rk4sys(int n, float h, float *X, int nsteps) {
  int i, step;
  float *Y, **F;

  F = calloc((RK_order+1), sizeof(float *));
  Y = calloc((n+1), sizeof(float));

  for (i = 0; i <= RK_order; i++){
    F[i] = calloc((n+1), sizeof(float));
  }

  for (step = 0; step < nsteps; step++) {
    virus_dynamic(X, F[1]);
    for (i = 0; i <= n; i++)
      Y[i] = X[i] + 0.5 * h * F[1][i];
    virus_dynamic(Y, F[2]);
    for (i = 0; i <= n; i++)
      Y[i] = X[i] + 0.5 * h * F[2][i];
    virus_dynamic(Y, F[3]);
    for (i = 0; i <= n; i++)
      Y[i] = X[i] + h * F[3][i];
    virus_dynamic(Y, F[4]);
    for (i = 0; i <= n; i++)
      X[i] = X[i] + (h/6) * (F[1][i] + 2 * F[2][i] + 2 * F[3][i] + F[4][i]);

    // for (i = 0; i <= n; i++){
    //   printf("step = %2d, X[%d] = %14f\n", step, i, X[i]);
    // }

    printf("%2d, %14f, %14f, %14f, %14f\n", step, X[0], X[1], X[2], X[3]);
  }

  for (i = 0; i <= RK_order; i++)
    free(F[i]);

  free(F);
  free(Y);

}

int main(int argc, char *argv[]) {

  const int n = 3; // 4 equations
  const int nsteps = 10000;
  // int i;
  float a = 0, b= 1.0, h;
  // float X[]= {10.0, 1.0, 0, 0};
  float X[]= {atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4])};

  h = (b - a) / nsteps;
  rk4sys(n, h, X, nsteps);

  return 0;
}



