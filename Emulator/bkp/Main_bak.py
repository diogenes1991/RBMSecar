import math
import array
import random
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
from scipy.special import rel_entr

import cosyHelp as cosy

#==============================================================================================================
#============================================ INPUT -- MODIFY HERE ============================================
#==============================================================================================================


FP = 'FP1' # Focal plane to consider


magnets = np.linspace(0.5,1.5,100)
# print(magnets)


#=============================================================================================================
# Main program
#=============================================================================================================

if __name__ == '__main__':

	print('Start running')

	for i in range(len(magnets)):
		cosy.COSYSingleRun(magnets, i)     # Testing of a single cosy run

	#cosy.generateMatrices(magnets, FP)    # Command to generate matrices for all magnet settings

	#beam0 = generateInitialDistribution(widthX, widthY, aX, aY)
	#KDE, POS =  generateKDE(saveReference = True, beam0 = beam0)
	#KDE, POS =  generateKDE(coeff = [[0, 1, 2.0]], beam0 = beam0)
	#plot.plotDistribution( KDE, POS, magnets = magnets, dimVD = dimVD )    # Displays KDE
	#plot.plotDistribution( magnets = magnets, dimVD = dimVD )    # Displays KDE
	#print ( compareKDE( coeff = [[0, 1, 1.294853909271475]], beam0 = beam0 ) )
	#Grid1DKL(beam0)
	#plot.plot1DKL('M33')
	#Grid2DKL(beam0)
	#plot.plot2DKL()

	# KL_minimize()
