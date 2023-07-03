# Functions to create/modify/run cosy files

import subprocess
from os.path import exists
import numpy as np

#Set the path to COSY executable
pathToCOSY='/Users/ruchigarg/Desktop/Codes/COSY/cosy'

#Name of original COSY file
cosyfilebase = 'SECAR_GammaOptics_' #the file you will execute

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
def COSYSingleRun(magnets, magnet_names, i=0, j=0):

	#Name of temp COSY file
	cosyfileTemp = 'SECAR_GammaOpticsTEMP' + ".fox"
	print("Running COSY i = ", i)

	cosyfile = cosyfilebase + magnet_names[j] + '.fox'

	replace_line(cosyfile, cosyfileTemp, 0, "INCLUDE '/Users/ruchigarg/Desktop/Codes/COSY/cosy';\n")

	if magnet_names[j]=='Q1':
		replace_line(cosyfileTemp, cosyfileTemp, 248-1, \
		'M5 0.250 Q1*SC*' + str(magnets[i]) + ' Q1H*SC*' + str(magnets[i]) + ' 0 -0.00318*Q1*SC*' + \
			str(magnets[i]) + ' 0 0.055;        {Q1+Hex}\n')
		 
	elif magnet_names[j]=='Q2':
		replace_line(cosyfileTemp, cosyfileTemp, 268-1, \
		'MQ 0.2979 Q2*SC*1.0*' + str(magnets[i]) + ' 0.068;                   			{Q2}\n')

	elif magnet_names[j]=='Q3':
		replace_line(cosyfileTemp, cosyfileTemp, 353-1, \
		'MQ 0.3499 Q3*SC*' + str(magnets[i]) + ' 0.11;                                    {Q3}')

	elif magnet_names[j]=='Q4':
		replace_line(cosyfileTemp, cosyfileTemp, 373-1, \
		'M5 0.3467 Q4*SC*' + str(magnets[i]) + ' 0 0 0 0 0.08;                           {Q4}')

	elif magnet_names[j]=='Q5':
		replace_line(cosyfileTemp, cosyfileTemp, 393-1, \
		'MQ 0.3466 Q5*SC*' + str(magnets[i]) + ' 0.06;                                   {Q5}')

	elif magnet_names[j]=='S1':
		replace_line(cosyfileTemp, cosyfileTemp, 333-1, \
		'MH 0.26 H1*SC*' + str(magnets[i]) + ' 0.11;                                   {HEX1}')


	#mkDir_cmd = 'mkdir "%s"' % (str(i))
	#subprocess.call(mkDir_cmd, shell=True)
	temp = open('temp.txt', "w") # temporary file for writing COSY terminal output
	subprocess.call([pathToCOSY, cosyfileTemp], stdout=temp)
	temp.close()

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
def transportTotal(x0, mapfile):
	coeff_array, power_array = read_map(mapfile)  # DL1
	x1 = transport(x0, coeff_array, power_array)
	return (x1)
