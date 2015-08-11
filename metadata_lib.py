from os import listdir
from os.path import isfile, join
from myutils import extractVariablesFromString
import pickle
from datetime import datetime

def transformParameters(parameters):

    for key in parameters:
        if key in ['N', 'n_steps', 'n_save', 'ftype', 'N_fix']:
            parameters[key] = int(parameters[key])
        elif key in ['nu_0','J','eta_n','phi_pack', 'dx', 'fanisotropy', 'max_dist']:
            parameters[key] = float(parameters[key])
        elif key in ['update_nn', 'chemoatt', 'rotate_coord']:
            parameters[key] = int(parameters[key]) == 1

    return parameters

def findFilesWithParameters(metadata, parameters):

    f = open(metadata, 'r')
    all_lines = f.readlines()
    heading = all_lines[0]
    lines = all_lines[1:]
    files = []
    dates = []
    for line in lines:
        params = {}
        data = line.split('\t')
        match = True
        i = 0
        while match and i<len(param_names):
            p_name = param_names[i]
            if not data[i] == ' ':
                if p_name in ['N', 'n_steps', 'n_save', 'ftype']:
                    match = int(data[i].strip()) == int(parameters[p_name])
                elif p_name in ['nu_0','J','eta_n','phi_pack', 'dt']:
                    match = float(data[i].strip()) == float(parameters[p_name])
                elif p_name in ['update_nn']:
                    if parameters[p_name] == ' ':
                        match = int(data[i].strip()) == 1
                    else:
                        match = int(data[i].strip()) == int(parameters[p_name])
                        
            i += 1

        for p_name in ftype2_param:
            if (int(parameters['ftype']) == 2) and (not data[i] == ' '):
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
    if isfile(outfile):
        ftoread= open(outfile, 'r')
        lines = ftoread.readlines()[1:]
        ftoread.close()
        files_in_metadata = ['', ]*len(lines)
        for i, line in enumerate(lines):
            files_in_metadata[i] = line.split('\t')[-1].strip()
        ftowrite=open(outfile, 'a')
    else:
        files_in_metadata = []
        ftowrite=open(outfile, 'w')
        # write the heading
        heading = ''
        for param in param_names:
            if heading == '':
                heading = param
            else:
                heading += '\t'+param

        for param in ftype2_param:
            heading += '\t'+param

        heading += '\tdate_string\tfilename\n'
        ftowrite.write(heading)

    for fname in allfiles:
        # include only those not in the list
        if not fname in files_in_metadata:
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

            for param in ftype2_param:
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