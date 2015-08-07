from os import listdir
from os.path import isfile, join
from myutils import extractVariablesFromString
import pickle
from datetime import datetime

param_names = ['N','nu_0','J','eta_n','phi_pack','ftype','n_steps','dt','n_save']

def findFilesWithParameters(metadata, parameters):

    f = open(metadata, 'r')
    lines = f.readlines()[1:]
    files = []
    dates = []
    for line in lines:
        params = {}
        data = line.split('\t')
        match = True
        i = 0
        while match and i<len(param_names):
            p_name = param_names[i]
            if p_name in ['N', 'n_steps', 'n_save', 'ftype']:
                match = int(data[i].strip()) == int(parameters[p_name])
            elif p_name in ['nu_0','J','eta_n','phi_pack', 'dt']:
                match = float(data[i].strip()) == float(parameters[p_name])
            i += 1
        if match:
            files.append(data[-1].strip())
            date_list = [int(x) for x in data[-2].split(',')]
            dates.append(datetime(*date_list))

    f.close()
    return files, dates

def generateMetadata(path):
	"""
	Generate the metadata for reference of the parameters used (database)
	"""
	# parameters
	file_start = 'sphere_data_analysis_'
	file_pattern = 'sphere_data_analysis_{y}-{m}-{d}_{h}-{mm}-{s}.p'
	outfile = join(path,'metadata.dat')

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

if __name__ == "__main__":
    """
    Run the main when externally called
    """
    generateMetadata('./data')