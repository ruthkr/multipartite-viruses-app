#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float uniform_dist(float, float);
void ssa_gillespie(float[], float *);

// Random Uniform Distribution
float uniform_dist(float min, float max) {
	rand(); // Fix for random
	float scale = rand() / (float) RAND_MAX; /* [0, 1.0] */
	return min + scale * ( max - min );      /* [min, max] */
}

void ssa_gillespie(float X[], float *params) {
	float r1 = uniform_dist(0, 1);
	float r2 = uniform_dist(0, 1);

	// Transition rates
	float *w = calloc(7, sizeof(float));
	float *wsum = calloc(7, sizeof(float));

	w[1] = params[2] * X[1];
	w[2] = params[0] * X[1] * X[1] * X[2] / pow(1000, 2);
	w[3] = params[3] * X[2];
	w[4] = params[1] * X[1] * X[2] / 1000;
	w[5] = params[1] * X[1];
	w[6] = params[0] * X[1] * X[2] / 1000;

	// Total Number of all of the transition rates above
	for (int i = 1; i <= 6; i++) {
		w[0] += w[i];
	}

	// Partial sums of the rates (to be evaluated in the conditionals)
	wsum[0] = 0;

	for (int i = 1; i <= 6; i++) {
		for (int j = 1; j <= i; j++) {
			wsum[i] += w[j];
		}
		wsum[i] /= w[0];
	}

	// Timestep tau
	float tau = 1/w[0] * log(1/r1);

	// Compute at time + tau
	if ( (wsum[0] <= r2) & (r2 < wsum[1]) ) {
		X[1] -= 1;
	} else if ( (wsum[1] <= r2) & (r2 < wsum[2]) ) {
		X[1] -= 1;
	} else if ( (wsum[2] <= r2) & (r2 < wsum[3]) ) {
		X[2] -= 1;
	} else if ( (wsum[3] <= r2) & (r2 < wsum[4]) ) {
		X[2] -= 1;
	} else if ( (wsum[4] <= r2) & (r2 < wsum[5]) ) {
		X[2] += 1;
	} else if ( (wsum[5] <= r2) & (r2 < wsum[6]) ) {
		X[1] += 1;
	}
	// Update time step
	X[0] += tau;

	free(w);
	free(wsum);
}


int main(int argc, char const *argv[]) {

	// Maximum Time
	float max_time = atof(argv[1]);

	// Initial Condition
	float X[] = {0, 0, 0};

	// Parameters
	float *params = calloc(4, sizeof(float));
	params[0] = atof(argv[4]);  // kappa default: 1
	params[1] = atof(argv[5]);  // omega default: 1
	params[2] = atof(argv[6]);  // gamma
	params[3] = atof(argv[7]);  // sigma default: 0.1

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
			// printf("Step: %ld\n", step);
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
	printf("%f\n", prob_exit);

	free(params);
	return 0;
}
