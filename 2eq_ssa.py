#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Ruth Kristianingsih <ruth.kristianingsih30@gmail.com>

from math import log, sqrt
#from random import random
import numpy as np

#MAX_TIME = 50
SELECT_FLAG = True

# Given variables
kappa = 1
alpha = 1
gamma = 0.1
sigma = 0.1


#def uniform_dist(min, max):
#	return min + (max - min) * random()

def ssa_gillespie(max_time, R_init, p_init):
	# Initialize populations
	R = R_init
	p = p_init
	# Initialize time
	time = 0.0
	step = 0

	with open("results-ssa.csv", "w") as file:
		# Print initial values to file
		file.write("%f %d %d\n" % (time, R, p))
		# Main loop
		while (time < max_time):
			# Random numbers
			r1 = np.random.uniform(0, 1)
			r2 = np.random.uniform(0, 1)
#			print(r1, r2)

			# Transition rates
			w1 = gamma * R
			w2 = sigma * p
			w3 = alpha * R
			w4 = alpha * R * p
			w5 = kappa * R * p
			w6 = kappa * R * R * p
			w0 = w1 + w2 + w3 + w4 + w5 + w6
			print(R, p, w0)

			# Timestep tau
			tau = 1/w0 * log(1/r1)

			# Compute at time + tau
			if 0 <= r2 < w1/w0:
				R -= 1
			elif w1/w0 <= r2 < (w1+w2)/w0:
				p -= 1
			elif (w1+w2)/w0 <= r2 < (w1+w2+w3)/w0:
				p += 1
			elif (w1+w2+w3)/w0 <= r2 < (w1+w2+w3+w4)/w0:
				p -= 1
			elif (w1+w2+w3+w4)/w0 <= r2 < (w1+w2+w3+w4+w5)/w0:
				R += 1
			elif (w1+w2+w3+w4+w5)/w0 <= r2 < (w1+w2+w3+w4+w5+w6)/w0:
				R -= 1

			# Update time step
			time += tau
			step += 1
			# Print values to file
			file.write("%f %d %d\n" % (time, R, p))
	# Close file
	while not file.closed:
		file.close()

ssa_gillespie(100, 0.5, 0)



#def qssa_gillespie(X1_init):
#	# Initialize populations
#	X1 = X1_init
#	n1 = X1/S
#	P = (n1**2)/( (n1**2) + nu11)
#	X2 = np.random.binomial(E, P)
#	# Initialize time
#	time = 0.0
#	step = 0
#
#	with open("results-qssa.csv", "w") as file:
#		# Print initial values to file
#		file.write("%f %d %d\n" % (time, X1, X2))
#		# Main loop
#		while (time < MAX_TIME):
#			# Compute P and X2
#			P = (n1**2)/( (n1**2) + nu11)
#			X2 = np.random.binomial(E, P)
#
#			# Random numbers
#			z1 = uniform_dist(0,1)
#			z2 = uniform_dist(0,1)
#
#			# Transition rates
#			w1 = E * (E * R + kappa1 * X2 )
#			w2 = E**2 * (kappa2 * n1)
#			w3 = (n1**2) * (E - X2) * E
#			w4 = nu11 * X2 * E
#			w0 = w1 + w2 + w3 + w4
#
#			# Timestep tau
#			tau = 1/w0 * log(1/z1)
#			if tau < 0:
#				print(tau, z1, "It's going back in time.")
#				break
#
#			# Compute at time + tau
#			if 0 <= z2 < w1/w0:
#				n1 += 1
#			elif w1/w0 <= z2 < (w1+w2)/w0:
#				n1 -= 1
#			elif (w1+w2)/w0 <= z2 < (w1+w2+w3)/w0:
#				n1 -= 2
#			elif (w1+w2+w3)/w0 <= z2 <= 1:
#				n1 += 2
#
#			# Update X1
#			X1 = n1 * S
#			# Update time step
#			time += tau
#			step += 1
#			# Print values to file
#			file.write("%f %d %d\n" % (time, X1, X2))
#	# Close file
#	while not file.closed:
#		file.close()


#assumption = input("Which simulation do you want to run? (SSA/QSSA): ")
#
#if assumption.upper() == "SSA" or assumption.upper() == "S":
#	# Stochastic Simulation Algorithm Gillespie
#	print("Simulating Stochastic Simulation Algorithm Gillespie...")
#	# Input initial values
#	if SELECT_FLAG:
#		X1_init = float(input("Select the initial population for X1: "))
#		X2_init = float(input("Select the initial population for X2: "))
#	else:
#		X1_init = 0
#		X2_init = 0
#	# Run the simulation
#	ssa_gillespie(X1_init, X2_init)
#	print("\nThe results can be found on results-ssa.csv")
#
#elif assumption.upper() == "QSSA" or assumption.upper() == "Q":
#	# Quasi-Steady State Assumption Gillespie
#	print("Simulating Quasi-Steady State Assumption Gillespie...")
#	# Input initial values
#	if SELECT_FLAG:
#		X1_init = float(input("Select the initial population for X1: "))
#	else:
#		X1_init = 0
#	# Run the simulation
#	qssa_gillespie(X1_init)
#	print("\nThe results can be found on results-qssa.csv")
#
#else:
#	print("That's not a valid answer")



