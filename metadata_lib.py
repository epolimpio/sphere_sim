#! /usr/bin/env python3

from os import listdir
from os.path import isfile, join
from myutils import extractVariablesFromString, readConfigFile
import pickle
from datetime import datetime

def transformParameters(parameters):

    for key in parameters:
        if key in ['N', 'n_steps', 'n_save', 'ftype', 'N_fix']:
            parameters[key] = int(parameters[key])
        elif key in ['nu_0','J','eta_n','phi_pack', 'dx', 'fanisotropy', 'max_dist', 'J_chemo']:
            parameters[key] = float(parameters[key])
        elif key in ['update_nn', 'chemoatt', 'rotate_coord']:
            parameters[key] = int(parameters[key]) == 1

    return parameters

def findFilesWithParameters(metadata, parameters):

    parameters = transformParameters(parameters)
    f = open(metadata, 'r')
    all_lines = f.readlines()
    heading = []
    for key in all_lines[0].split('\t')[:-2]:
        heading.append(key.strip())
    lines = all_lines[1:]
    files = []
    dates = []
    for line in lines:
        param_comp = {}
        all_data = line.split('\t')
        for i, data in enumerate(all_data[:-2]):
            param_comp[heading[i]] = data.strip()

        match = True
        for key in ['N', 'n_steps', 'n_save', 'ftype', 'N_fix']:
            match = match and int(param_comp[key]) == parameters[key]

        for key in ['nu_0','J','eta_n','phi_pack', 'dx']:
            match = match and float(param_comp[key]) == parameters[key]

        for key in ['outfile_video','outfile_analysis','outfile_postrotation']:
            match = match and param_comp[key] == parameters[key]

        for key in ['update_nn', 'chemoatt', 'rotate_coord']:
            if parameters[key]:
                match = match and param_comp[key].lower() == 'true'
            else:
                match = match and param_comp[key].lower() == 'false'

        if parameters['chemoatt']:
            match = match and param_comp['chemo_points'] == parameters['chemo_points']
            match = match and float(param_comp['J_chemo']) == parameters['J_chemo']

        if parameters['ftype'] == 2:
            match = match and float(param_comp['fanisotropy']) == parameters['fanisotropy']
            match = match and float(param_comp['max_dist']) == parameters['max_dist']

        if match:
            files.append(all_data[-1].strip())
            date_list = [int(x) for x in all_data[-2].strip().split(',')]
            dates.append(datetime(*date_list))

    f.close()
    return files, dates

def generateMetadata(path, parameters_file = './parameters.ini'):
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
        all_lines = ftoread.readlines()
        heading = []
        for key in all_lines[0].split('\t')[:-2]:
            heading.append(key.strip())
        lines = all_lines[1:]
        ftoread.close()
        files_in_metadata = ['', ]*len(lines)
        for i, line in enumerate(lines):
            files_in_metadata[i] = line.split('\t')[-1].strip()
        ftowrite=open(outfile, 'a')
    else:
        files_in_metadata = []
        ftowrite=open(outfile, 'w')
        # write the heading
        heading = []
        heading_str = ''
        param_names = readConfigFile(parameters_file)
        for param in param_names:
            heading.append(param)
            if heading_str == '':
                heading_str = param
            else:
                heading_str += '\t'+param

        heading_str += '\tdate_string\tfilename\n'
        ftowrite.write(heading_str)

    for fname in allfiles:
        # include only those not in the list
        if not fname in files_in_metadata:
            data=pickle.load(open(join(path,fname), "rb" ) )
            parameters = data[0]
            line = ''
            for param in heading:
                string = str(parameters[param])

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