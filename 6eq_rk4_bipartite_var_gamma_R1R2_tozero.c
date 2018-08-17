#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define RK_order 6

/* define prototypes */

void virus_dynamic(double *, double *, float *);
void rk4sys(int, float, float *, double *, int);


void virus_dynamic(double *X, double *f, float *param) {
	f[0] = param[0] * X[0] * X[2] * (1 - X[0] - X[1]) - param[8] * X[0] * X[3] - param[4] * X[0];
	f[1] = param[1] * X[1] * X[2] * (1 - X[0] - X[1]) - param[9] * X[1] * X[3] - param[5] * X[1];
	f[2] = param[2] * X[0] * (1 - X[2] - X[3]) - param[6] * X[2];
	f[3] = param[3] * X[1] * (1 - X[2] - X[3]) - param[7] * X[3];
	f[4] = param[8] * X[0] * X[3] - param[10] * X[4];
	f[5] = param[9] * X[1] * X[3] - param[11] * X[5];
}

void rk4sys(int n, float h, float *param, double *X, int nsteps) {
	int i, step;
	double *Y, **F;
	// int tozero

	F = calloc((RK_order+1), sizeof(double *));
	Y = calloc((n+1), sizeof(double));

	for (i = 0; i <= RK_order; i++){
		F[i] = calloc((n+1), sizeof(double));
	}

	for (step = 0; step < nsteps; step++) {
		virus_dynamic(X, F[1], param);
		for (i = 0; i <= n; i++)
			Y[i] = X[i] + 0.5 * h * F[1][i];
		virus_dynamic(Y, F[2], param);
		for (i = 0; i <= n; i++)
			Y[i] = X[i] + 0.5 * h * F[2][i];
		virus_dynamic(Y, F[3], param);
		for (i = 0; i <= n; i++)
			Y[i] = X[i] + h * F[3][i];
		virus_dynamic(Y, F[4], param);
		for (i = 0; i <= n; i++)
			X[i] = X[i] + (h/6) * (F[1][i] + 2 * F[2][i] + 2 * F[3][i] + F[4][i]);

		// for (i = 0; i <= n; i++){
		//   printf("step = %2d, X[%d] = %14f\n", step, i, X[i]);
		// }
	}

	for (i = 0; i <= RK_order; i++)
		free(F[i]);

	free(F);
	free(Y);

	// return(tozero);
}

int main(int argc, char *argv[]) {

	const int n = 5; // 4 equations
	const int nsteps = atoi(argv[1]);
	const float h = 0.1;

	float gamma = 0.0;
	float R1_init = 0.0;
	float R2_init = 0.0;

	int partitions_gamma = atoi(argv[20]); // default: 55
	int partitions_R1 = atoi(argv[21]);    // default: 1000
	int partitions_R2 = atoi(argv[22]);    // default: 10

	int j, k1, k2;
	float gamma_step = (0.55 - gamma)/ partitions_gamma;
	float R1_step = (1.0 - R1_init) / partitions_R1;
	float R2_step = (1.0 - R2_init) / partitions_R2;

	// Parameters
	float *params = calloc(12, sizeof(float));

	for (j = 0; j <= partitions_gamma; j++) {
		// Initial conditions
		fprintf(stderr, "Gamma = %0.3f\n", gamma);
		for (k2 = 0; k2 <= partitions_R2; k2++) {
		// R2 loop
			for (k1 = 0; k1 <= partitions_R1; k1++) {
			// R1 loop
				// Initialize values
				double X[]= {R1_init, R2_init, atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7])}; // R1.init, R2.init, p.init, s.init, v1.init, v2.init

				params[0] = atof(argv[8]);   // kappa1
				params[1] = atof(argv[9]);   // kappa2
				params[2] = atof(argv[10]);  // alpha
				params[3] = atof(argv[11]);  // beta
				params[4] = gamma;           // gamma1
				params[5] = gamma;           // gamma2
				params[6] = atof(argv[14]);  // sigma1
				params[7] = atof(argv[15]);  // sigma2
				params[8] = atof(argv[16]);  // epsilon1
				params[9] = atof(argv[17]);  // epsilon2
				params[10] = atof(argv[18]); // delta1
				params[11] = atof(argv[19]); // delta2

				rk4sys(n, h, params, X, nsteps);

				printf("%.3f, %.3f, %.2f, %.15f, %.15f\n", gamma, R1_init, R2_init, X[0], X[1]);

				// Update R1
				R1_init += R1_step;
			} // end R1 loop
			R1_init = 0;
			// Update R2
			R2_init += R2_step;
		} // end R2 loop
		R2_init = 0;
		// Update gamma
		gamma += gamma_step;
	} // end gamma loop

	free(params);

	return 0;
}



