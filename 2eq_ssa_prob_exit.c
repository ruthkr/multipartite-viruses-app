#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float uniform_dist(float, float);
void ssa_gillespie(float *, float *);

// Random Uniform Distribution
float uniform_dist(float min, float max) {
	rand(); // Fix for random
	float scale = rand() / (float) RAND_MAX; /* [0, 1.0] */
	return min + scale * ( max - min );      /* [min, max] */
}

void ssa_gillespie(float *X, float *params) {
	float r1 = uniform_dist(0, 1);
	float r2 = uniform_dist(0, 1);

	// Transition rates
	float w1 = params[2] * X[1];
	float w2 = params[0] * X[1] * X[1] * X[2] / pow(1000, 2);
	float w3 = params[3] * X[2];
	float w4 = params[1] * X[1] * X[2] / 1000;
	float w5 = params[1] * X[1];
	float w6 = params[0] * X[1] * X[2] / 1000;
	float w0 = w1 + w2 + w3 + w4 + w5 + w6;

	// Timestep tau
	float tau = 1/w0 * log(1/r1);

	// printf("%f, %f\n", r2, w0);

	// Compute at time + tau
	if ((0 <= r2) & (r2 < w1/w0)) {
		X[1] -= 1;
	} else if ( (w1/w0 <= r2) & (r2 < (w1+w2)/w0) ) {
		X[1] -= 1;
	} else if ( ((w1+w2)/w0 <= r2) & (r2 < (w1+w2+w3)/w0)) {
		X[2] -= 1;
	} else if ( ((w1+w2+w3)/w0 <= r2) & (r2 < (w1+w2+w3+w4)/w0)) {
		X[2] -= 1;
	} else if ( ((w1+w2+w3+w4)/w0 <= r2) & (r2 < (w1+w2+w3+w4+w5)/w0)) {
		X[2] += 1;
	} else if ( ((w1+w2+w3+w4+w5)/w0 <= r2) & (r2 < (w1+w2+w3+w4+w5+w6)/w0)) {
		X[1] += 1;
	}
	// Update time step
	X[0] += tau;
}


int main(int argc, char const *argv[]) {

	// Maximum Time
	float max_time = atof(argv[1]);

	// Initial Condition
	float X[] = {0, 0, 0};

	// Parameters
	float *params = calloc(4, sizeof(float));
	params[0] = atof(argv[4]);  // kappa
	params[1] = atof(argv[5]);  // omega
	params[2] = atof(argv[6]);  // gamma
	params[3] = atof(argv[7]);  // sigma

	// Seed
	srand((unsigned)time(NULL));

	int j;
	int count_prob = 0;
	int nsim = atoi(argv[8]);

	// Loop
	for(j = 0; j < nsim; j++) {
		unsigned long int step = 0;

		// Reset Initial Condition
		X[0] = 0;
		X[1] = atof(argv[2]);
		X[2] = atof(argv[3]);

		while (X[0] < max_time) {
			ssa_gillespie(X, params);
			step += 1;
			// printf("%f, %f, %f\n", X[0], X[1], X[2]);
		}

		if ((X[1] <= 0) & (X[2] <= 0)) {
			count_prob += 1;
		}
		// printf("%f, %d, %d, %lu\n", params[2], count_prob, j, step);
	}

	float prob_exit = (float) count_prob / nsim;
	printf("%f, %f\n", params[2], prob_exit);

	return 0;
}
