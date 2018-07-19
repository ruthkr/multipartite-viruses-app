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

	// X[1] = R1, X[2] = R2, X[3] = p, X[4] = s, X[5] = v1, X[6] = v2

	// Transition rates
	// Reaction Rates for RNA1
	float w1 = params[0] * X[1] * X[3] / params[6];
	float w2 = params[0] * X[1] * X[1] * X[3] / pow(params[6],2);
	float w3 = params[0] * X[1] * X[2] * X[3] / pow(params[6],2);
	float w4 = params[4] * X[1] * X[4] / params[6];
	float w5 = params[2] * X[1];

	// Reaction Rates for RNA2
	float w6 = params[0] * X[2] * X[3] / params[6];
	float w7 = params[0] * X[2] * X[2] * X[3] / pow(params[6],2);
	float w8 = params[0] * X[1] * X[2] * X[3] / pow(params[6],2);
	float w9 = params[4] * X[2] * X[4] / params[6];
	float w10 = params[2] * X[2];

	// Reaction Rates for viral replicase
	float w11 = params[1] * X[1];
	float w12 = params[1] * X[1] * X[3] / params[6];
	float w13 = params[1] * X[1] * X[4] / params[6];
	float w14 = params[3] * X[3];

	// Reaction Rates for coat protein
	float w15 = params[1] * X[2];
	float w16 = params[1] * X[2] * X[3] / params[6];
	float w17 = params[1] * X[2] * X[4] / params[6];
	float w18 = params[3] * X[4];

	// Reaction Rates for virion1
	float w19 = params[4] * X[1] * X[4]/ params[6];
	float w20 = params[5] * X[5];

	// Reaction Rates for virion2
	float w21 = params[4] * X[2] * X[4]/ params[6];
	float w22 = params[5] * X[6];

	// Total Number of all of the transition rates above
	float w0 = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20 + w21 + w22;

	// Timestep tau
	float tau = 1/w0 * log(1/r1);

	// printf("%f, %f\n", r2, w0);

	// Compute at time + tau
	if ((0 <= r2) & (r2 < w1 / w0)) {
		X[1] += 1;
	} else if ( (w1 / w0 <= r2) & (r2 < (w1 + w2) / w0) ) {
		X[1] -= 1;
	} else if ( ((w1 + w2)/w0 <= r2) & (r2 < (w1 + w2 + w3) / w0)) {
		X[1] -= 1;
	} else if ( ((w1 + w2 + w3)/w0 <= r2) & (r2 < (w1 + w2 + w3 + w4) / w0)) {
		X[1] -= 1;
	} else if ( ((w1 + w2 + w3 + w4) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5) / w0)) {
		X[1] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6) / w0)) {
		X[2] += 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7) / w0)) {
		X[2] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8) / w0)) {
		X[2] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 +w5 + w6 + w7 + w8 + w9) / w0)) {
		X[2] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10) / w0)) {
		X[2] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10) / w0 <= r2) & (r2 < (w1 + w2 +w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11) / w0)) {
		X[3] += 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12) / w0)) {
		X[3] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13) / w0)) {
		X[3] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14) / w0)) {
		X[3] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15) / w0)) {
		X[4] += 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15)/w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16) / w0)) {
		X[4] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16)/w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17) / w0)) {
		X[4] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 +w12 + w13 + w14 + w15 + w16 + w17 + w18) / w0)) {
		X[4] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17 + w18) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 +w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19) / w0)) {
		X[5] += 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 +w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20) / w0)) {
		X[5] -= 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 +w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20 + w21) / w0)) {
		X[6] += 1;
	} else if ( ((w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 + w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20 + w21) / w0 <= r2) & (r2 < (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10 + w11 +w12 + w13 + w14 + w15 + w16 + w17 + w18 + w19 + w20 + w21 + w22) / w0)) {
		X[6] -= 1;
	}
	// Update time step
	X[0] += tau;
}


int main(int argc, char const *argv[]) {
	// Parameters
	float *params = calloc(7, sizeof(float));
	params[0] = atof(argv[8]);  // kappa
	params[1] = atof(argv[9]);  // omega
	params[2] = atof(argv[10]);  // gamma
	params[3] = atof(argv[11]);  // sigma
	params[4] = atof(argv[12]);  // epsilon
	params[5] = atof(argv[13]);  // delta
	params[6] = atof(argv[14]);  // current capacity

	srand((unsigned)time(NULL));

	float X[] = { 0, atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]) };
	float max_time = atof(argv[1]);

	printf("%f, %f, %f, %f, %f, %f, %f\n", X[0], X[1], X[2], X[3], X[4], X[5], X[6]);

	unsigned long int step = 0;

	while (X[0] < max_time) {
		ssa_gillespie(X, params);
		step += 1;
		printf("%f, %f, %f, %f, %f, %f, %f\n", X[0], X[1], X[2], X[3], X[4], X[5], X[6]);
	}

	return 0;
}


