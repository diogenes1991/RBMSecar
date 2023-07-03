#%%
import math
import array
import random
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.special import rel_entr
import subprocess
import matplotlib.pyplot as plt

import cosyHelp as cosy

import time

#==============================================================================================================
#============================================ MAGNET INPUT ====================================================
#==============================================================================================================

n = 100
magnets = np.linspace(1.0-0.5, (1.0+0.5-1/n), n)
magnets_names = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'B1', 'B2', 'S1']
# print(magnets)

#============================================ BEAM INPUT ======================================================

# Distribution at target
xh = 0 # center of x distribution in mm
yh = 0 # center of x distribution in mm
widthX = 1.5/2 # half x- beam spot in mm
widthY = 1.5/2 # half y- beam spot in mm
dE = 0.0025 # half of total energy spread in fraction 1% = 0.01
#dE = 0.00 # half of total energy spread in fraction 1% = 0.01
aX = 5 # mrad
aY = 5 # mrad

numberMC = 1000 # number of beam particles. Was 800

#=============================================================================================================
# Function that generates positions and angles at the target within the input parameters
#=============================================================================================================

def generateInitialDistribution(widthX, widthY, aX, aY):
	beam = np.zeros([numberMC, 8], dtype="float64", order="F")

	for j in range(numberMC):
		# Sampling within ellipse of possible positions
		r = widthX * np.sqrt(random.uniform(0, 1))
		theta = random.uniform(0, 1) * 2 * np.pi
		x = xh + r * np.cos(theta)
		y = yh + widthY/widthX *  r * np.sin(theta)
		beam[j][0] = x/1000 # in m
		beam[j][2] = y/1000 # in m

		# Sampling within ellipse of possible angles
		r = aX * np.sqrt(random.uniform(0, 1))
		theta = random.uniform(0, 1) * 2 * np.pi
		angleX = r * np.cos(theta)
		if aX != 0:
			angleY = aY/aX *  r * np.sin(theta)
		else:
			angleY = aY * np.sqrt(random.uniform(0, 1))
		beam[j][1] = angleX/1000 # in rad
		beam[j][3] = angleY/1000 # in rad
		beam[j][5] = dE
	return (beam)


#=============================================================================================================
# Main program
#=============================================================================================================

if __name__ == '__main__':

	print('Start running')

	## Print Matrix elements for the magnets in magnet_names list and for scaling factors in magnets list.
	# for i in range(len(magnets_names)):
	# 	if magnets_names[i]== 'B1' or magnets_names[i]== 'B2':
	# 		start_time = time.time()
	# 		cosy.COSYSingleRun([1.0], magnets_names, 0, i)
	# 		mv_cmd = 'mv "%s" "%s"' % ("fort.50",  'Results/{}_{:.2f}'.format(magnets_names[i], 1.0))
	# 		subprocess.call(mv_cmd, shell=True)
	# 		print(i, "--- %s seconds ---" % (time.time() - start_time))
	# 	else:
	# 		for j in range(len(magnets)):
	# 			start_time = time.time()
	# 			cosy.COSYSingleRun(magnets, magnets_names, j, i)
	# 			mv_cmd = 'mv "%s" "%s"' % ("fort.50",  'Results/{}_{:.2f}'.format(magnets_names[i], magnets[j]))
	# 			subprocess.call(mv_cmd, shell=True)
	# 			print(i, "--- %s seconds ---" % (time.time() - start_time))

	## x, ax, y, ay, dL, dE, dM, dZ
	beam_in = generateInitialDistribution(widthX, widthY, aX, aY)
	ax_in = beam_in.transpose()[3]
	plt.hist(ax_in, color='blue', bins=np.linspace(-5e-2, 5e-2, 500))

	# test_scale = 0.5626
	# cosy.COSYSingleRun([test_scale], ['Q1'], 0, 0) # args = magnet scale values array, magnet name array, index in magent scale array, index in magnet name array
	# mv_cmd = 'mv "%s" "%s"' % ("fort.50",  'test_map0.txt')
	# subprocess.call(mv_cmd, shell=True)
	# beam_out = cosy.transportTotal(beam_in, 'test_map0.txt') # args = input beam distribution array, transport map file name
	# ax_out = beam_out.transpose()[3]
	# plt.hist(ax_out, color='k', bins=np.linspace(-5e-2, 5e-2, 500), alpha=0.5)

	beam_out = cosy.transportTotal(beam_in, '0.5626')  # args = input beam distribution array, transport map file name
	ax_out = beam_out.transpose()[3]
	plt.hist(ax_out, color='r', bins=np.linspace(-5e-2, 5e-2, 500))

	plt.show()

