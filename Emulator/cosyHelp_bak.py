# Functions to create/modify/run cosy files

import subprocess
from os.path import exists
import multiprocessing
import numpy as np

#Set the path to COSY executable
pathToCOSY='/Users/fernandomontes/Research/COSY/cosy'

#Name of original COSY file
cosyfile = 'SECAR_GammaOptics.fox' #the file you will execute

#Name of COSY output file
cosyoutput = 'output'


# Function to replace of text in a file
def replace_line(file_nameOrig, file_nameTemp, line_num, text):
	lines = open(file_nameOrig, 'r').readlines()
	lines[line_num] = text
	out = open(file_nameTemp, 'w')
	out.writelines(lines)
	out.write('\n')
	out.close()

# Function for adding a multiple lines of text in a file
def add_lines(file_nameOrig, file_nameTemp, line_num, text):
	lines = open(file_nameOrig, 'r').readlines()
	for i in range(numberMC):
		lines.append('\n')
	lines[(line_num+len(text)-1):] = lines[(line_num-1):]
	for i in range(len(text)):
		lines[line_num-1+i] = text[i]
	out = open(file_nameTemp, 'w')
	out.writelines(lines)
	out.write('\n')
	out.close()

# Function that uses COSY to generate matrices in a folder
def COSYSingleRun(magnets, i = 0, FP = 'FP1No'):

	#Name of temp COSY file
	cosyfileTemp = 'SECAR_GammaOpticsTEMP' + str(i) + ".fox"
	print("Running COSY i = ", i)

	replace_line(cosyfile, cosyfileTemp, 0, "INCLUDE '../../COSY';\n")
	replace_line(cosyfileTemp, cosyfileTemp, 248, \
		'M5 0.250 Q1*SC*' + str(magnets[i]) + ' Q1H*SC*' + str(magnets[i]) + ' 0 -0.00318*Q1*SC*' + \
			str(magnets[i]) + ' 0 0.055;        {Q1+Hex}\n')

	#mkDir_cmd = 'mkdir "%s"' % (str(i))
	#subprocess.call(mkDir_cmd, shell=True)
	temp = open('temp.txt', "w") # temporary file for writing COSY terminal output
	subprocess.call([pathToCOSY, cosyfileTemp], stdout=temp)
	mv_cmd = 'mv "%s" "%s"' % ("fort.50",  'Results/' + str(magnets[i]))
	subprocess.call(mv_cmd, shell=True)
	temp.close()

# Parallelization of COSY runs for the different magnet settings
def COSYParallelRun(magnets, FP):
	try:
		pool = multiprocessing.Pool()  # Take as many processes as possible
	except:
		for c in multiprocessing.active_children():
			os.kill(c.pid, signal.SIGKILL)
		pool = multiprocessing.Pool(4)  # Number of processes to be used
	for i in range(len(magnets)):
		pool.apply_async( COSYSingleRun, [magnets, i] )
	pool.close()
	pool.join()

def generateMatrices(magnets, FP):
	COSYParallelRun(magnets, FP)










#####################

# Read matrices and generate coefficients and power arrays
def read_map(filename):
	"""
	The full precision cosy maps are a pain in the ass to read.
	"""

	total_coeff = []
	total_power = []
	max_len = 0

	with open(filename, "r") as f:
		coeff = []
		power = []

		for _ in range(3):
			next(f)
		for line in f:
			if "---" in line:
				if len(coeff) > max_len:
					max_len = len(coeff)
				total_coeff.append(coeff)
				total_power.append(power)
				coeff = []
				power = []
			elif "I" in line:
				pass
			else:
				line = line.split()
				temp = line[3:11]
				coeff.append(float(line[1]))
				power.append(np.asarray(temp, dtype="float64"))
	# Make this a square array, pad with zeros
	coeff_array = np.zeros((8, max_len))
	power_array = np.zeros((8, max_len, 8))

	for i in range(8):
		length = len(total_coeff[i])
		coeff_array[i, :length] = total_coeff[i]
		power_array[i, :length, :] = total_power[i]

	coeff_array = np.asarray(coeff_array)
	power_array = np.asarray(power_array)

	return coeff_array, power_array

# Function that performs multiplication of transport matrix by beam array
# x0 has dimensions [number of samples, beam array dimension]
def transport(x0, map_coeff, map_power):
	x = np.zeros(x0.shape, dtype="float64", order="F")
	samples = x0.shape[0]
	n = map_coeff.shape[0]  # Beam array dimensions
	m = map_coeff.shape[1]  # Number of coefficients

	for l in range(samples):  # Loop over the samples
		for i in range(n):  # Loop over beam array dimensions
			j = 0
			while j < m and map_coeff[i][j] != 0:  # Loop over coefficients
				power_temp = 1
				for k in range(n):    # Loop over beam array dimensions
					power_temp = power_temp*x0[l][k]**map_power[i][j][k]
				x[l][i] = x[l][i] + map_coeff[i][j]*power_temp
				j = j + 1
	return x

# Performs beam transport to end given a beam array
def transportTotal(x0, i, coeff):
	coeff_array, power_array = read_map(str(i) + '/fort.50')  # DL1
	x1 = transport(x0, coeff_array, power_array)
	coeff_array, power_array = read_map(str(i) + '/fort.51')  # Q1
	for c in coeff:  # Modifying transport matrix element by element
		coeff_array[c[0]][c[1]] = c[2]
	x2 = transport(x1, coeff_array, power_array)
	coeff_array, power_array = read_map(str(i) + '/fort.52')  # toFP1
	x3 = transport(x2, coeff_array, power_array)
	return (x3)
