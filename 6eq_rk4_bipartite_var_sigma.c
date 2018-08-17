#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define RK_order 6

/* define prototypes */

void virus_dynamic(float *, float *, float *);
void rk4sys(int, float, float *, float *, int);


void virus_dynamic(float *X, float *f, float *param) {
	f[0] = param[0] * X[0] * X[2] * (1 - X[0] - X[1]) - param[8] * X[0] * X[3] - param[4] * X[0];
	f[1] = param[1] * X[1] * X[2] * (1 - X[0] - X[1]) - param[9] * X[1] * X[3] - param[5] * X[1];
	f[2] = param[2] * X[0] * (1 - X[2] - X[3]) - param[6] * X[2];
	f[3] = param[3] * X[1] * (1 - X[2] - X[3]) - param[7] * X[3];
	f[4] = param[8] * X[0] * X[3] - param[10] * X[4];
	f[5] = param[9] * X[1] * X[3] - param[11] * X[5];
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
	}

	for (i = 0; i <= RK_order; i++)
		free(F[i]);

	free(F);
	free(Y);

}

int main(int argc, char *argv[]) {

	const int n = 5; // 4 equations
	const int nsteps = atoi(argv[1]);
	const float h = 0.1;

	float sigma = 0;
	float R1_init = 0;
	float R2_init = 0;
	float p_init = 0;

	// bool tozero;

	float a, b, c;
	float p_minus, R1_minus, R2_minus;
	float R1_plus, R2_plus, p_plus;

	int partitions_sigma = atoi(argv[20]);
	int partitions_vars = atoi(argv[21]);
	int j, k1, k2, k3;
	float sigma_step = 1.0/partitions_sigma;
	float vars_step = 1.0/partitions_vars;

	// Parameters
	float *params = calloc(12, sizeof(float));

	for (j = 0; j <= partitions_sigma; j++) {
		// Initial conditions
		for (k3 = 0; k3 <= partitions_vars; k3++) {
		// P loop
			for (k2 = 0; k2 <= partitions_vars; k2++) {
			// R2 loop
				for (k1 = 0; k1 <= partitions_vars; k1++) {
				// R1 loop
					// Initialize values
					float X[]= {R1_init, R2_init, p_init, atof(argv[5]), atof(argv[6]), atof(argv[7])}; // R1.init, R2.init, p.init, s.init, v1.init, v2.init

					params[0] = atof(argv[8]);   // kappa1
					params[1] = atof(argv[9]);   // kappa2
					params[2] = atof(argv[10]);  // alpha
					params[3] = atof(argv[11]);  // beta
					params[4] = atof(argv[12]);  // gamma1
					params[5] = atof(argv[13]);  // gamma2
					params[6] = sigma;  // sigma1
					params[7] = sigma;  // sigma2
					params[8] = atof(argv[16]);  // epsilon1
					params[9] = atof(argv[17]);  // epsilon2
					params[10] = atof(argv[18]); // delta1
					params[11] = atof(argv[19]); // delta2

					rk4sys(n, h, params, X, nsteps);

					// Fixed points separatrix (Formulation 1)
					a = (params[6] / params[2] + 1 ) * params[0];
					b = (params[6] / params[2]) * params[0] * X[3] - params[0] * (1 - X[3]) - params[8] * X[3] - params[4];
					c = (params[8] * X[3] + params[4]) * (1 - X[3]);

					p_minus = (- b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
					R1_minus = (params[0] * p_minus - params[8] * X[3] - params[4]) / (params[0] * (p_minus + X[3]));
					R2_minus = R1_minus * X[3] / p_minus;

					p_plus = (- b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
					R1_plus = (params[0] * p_plus - params[8] * X[3] - params[4]) / (params[0] * (p_plus + X[3]));
					R2_plus = R1_plus * X[3] / p_plus;

					printf("%6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f, %6f\n", sigma, R1_init, R2_init, p_init, X[0], X[1], X[2], X[3], R1_plus, R2_plus, p_plus, R1_minus, R2_minus, p_minus);

					// Update R1
					R1_init += vars_step;
				} // end R1 loop
				R1_init = 0;
				// Update R2
				R2_init += vars_step;
			} // end R2 loop
			R2_init = 0;
			// Update p
			p_init += vars_step;
		} // end p loop
		p_init = 0;
		// Update sigma
		sigma += sigma_step;
	} // end sigma loop

	free(params);

	return 0;
}



