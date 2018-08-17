#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define RK_order 4

/* define prototypes */

void virus_dynamic(float *, float *, float *);
void rk4sys(int, float, float *, float *, int);


void virus_dynamic(float *X, float *f, float *param) {
	f[0] = param[0] * X[0] * X[1] * (1 - X[0]) - param[2] * X[0];
	f[1] = param[1] * X[0] * (1 - X[1]) - param[3] * X[1];
}

void rk4sys(int n, float h, float *param, float *X, int nsteps) {
	int i, step;
	float *Y, **F;
	// int tozero

	F = calloc((RK_order+1), sizeof(float *));
	Y = calloc((n+1), sizeof(float));

	for (i = 0; i <= RK_order; i++){
		F[i] = calloc((n+1), sizeof(float));
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

}

int main(int argc, char *argv[]) {

	const int n = 1; // 2 equations
	const int nsteps = atoi(argv[1]);
	const float h = 0.1;

	float gamma = 0;
	float R1_init = 0;
	float p_init = 0;

	int partitions_gamma = atoi(argv[8]); // default: 55
	int partitions_R1 = atoi(argv[9]);    // default: 1000
	int partitions_p = atoi(argv[10]);    // default: 10

	int j, k1, k2;
	float gamma_step = 0.55/partitions_gamma;
	float R1_step = 1.0/partitions_R1;
	float p_step = 1.0/partitions_p;

	// Parameters
	float *params = calloc(4, sizeof(float));

	for (j = 0; j <= partitions_gamma; j++) {
		// Initial conditions
		for (k2 = 0; k2 <= partitions_p; k2++) {
		// P loop
			for (k1 = 0; k1 <= partitions_R1; k1++) {
			// R1 loop
				// Initialize values
				float X[]= {R1_init, p_init}; // R1.init, R2.init, p.init, s.init, v1.init, v2.init

				params[0] = atof(argv[4]);   // kappa
				params[1] = atof(argv[5]);   // omega
				params[2] = gamma;           // gamma
				params[3] = atof(argv[7]);   // sigma

				rk4sys(n, h, params, X, nsteps);

				printf("%.3f, %.3f, %.2f, %6f, %6f\n", gamma, R1_init, p_init, X[0], X[1]);

				// Update R1
				R1_init += R1_step;
			} // end R1 loop
			R1_init = 0;
			// Update p
			p_init += p_step;
		} // end p loop
		p_init = 0;
		// Update gamma
		gamma += gamma_step;
	} // end gamma loop

	free(params);

	return 0;
}



