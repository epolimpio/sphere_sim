from os import listdir
from os.path import isfile, join
from myutils import extractVariablesFromString
import pickle

# parameters
path = './data'
file_start = 'sphere_data_analysis_'
file_pattern = 'sphere_data_analysis_{y}-{m}-{d}_{h}-{mm}-{s}.p'
outfile = './data/metadata.dat'
param_names = ['N','nu_0','J','eta_n','phi_pack','ftype','n_steps','dt','n_save']

# get all the files in the folder
allfiles = [ f for f in listdir(path)
			if (isfile(join(path,f)) and f[:len(file_start)] == file_start)]

# open the metadata file
ftowrite = open(outfile, 'w')

# write the heading
heading = ''
for param in param_names:
	if heading == '':
		heading = param
	else:
		heading += '\t'+param

heading += '\tdate_string\tfilename\n'
ftowrite.write(heading)

for fname in allfiles:
	data=pickle.load(open(join(path,fname), "rb" ) )
	parameters = data[0]
	line = ''
	for param in param_names:
		if param in parameters:
			string = str(parameters[param])
		else:
			string = ' '

		if line == '':
			line = string
		else:
			line += '\t'+string

	datetime_dict = extractVariablesFromString(fname, file_pattern)
	date_string = '{y},{m},{d},{h},{mm},{s}'.format(**datetime_dict)
	line += '\t' + date_string + '\t' + fname + '\n'
	ftowrite.write(line)

ftowrite.close()