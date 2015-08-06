import numpy as np
import re
import os

def readConfigFile(path):
    """
    Read the standardized configuration file of this
    application.
    The output is a dictionary with the name of the variables
    and the string value to each of them
    """
    f = open(path, 'r')
    out = {}
    for line in f.readlines():
        if not line in ['\n', '\r\n']:
            if not line[0] == '#':
                aux = line.rstrip('\n').split('=')
                if len(aux) == 2:
                    name = aux[0].strip()
                    value = aux[1].strip()
                    out[name] = value
    f.close()
    return out

def calcRotationMatrix(axis_vec, angle):
    """
    Calculate the rotation matrix around an arbitrary axis
    """

    L = np.sqrt(np.sum(axis_vec**2,0))
    axis_vec = axis_vec/L

    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    u = axis_vec[0]
    v = axis_vec[1]
    w = axis_vec[2]

    R = np.array([[(u**2)*(1-cos_theta)+cos_theta, u*v*(1-cos_theta)-w*sin_theta, u*w*(1-cos_theta)+v*sin_theta],
                [u*v*(1-cos_theta)+w*sin_theta, (v**2)*(1-cos_theta)+cos_theta, v*w*(1-cos_theta)-u*sin_theta],
                [u*w*(1-cos_theta)-v*sin_theta, v*w*(1-cos_theta)+u*sin_theta, (w**2)*(1-cos_theta)+cos_theta]])

    return R

def extractVariablesFromString(string, pattern):
    """
    Given a pattern we extract the variables from string.
    Example:
    string = bbbb_aa_06_07_xy 
    patern = bbbb_{0}_{1}_{2}_{3}
    return {'0':'aa','1':'06','2':'07','3':'xy'}
    """

    # First we get the places of the variables in the string
    regPattern = r'\{([^}]+)\}'
    indexes = re.findall(regPattern, pattern)

    # for each of the places we write a group in the pattern
    newPattern = pattern
    for index in indexes:
        index_string = r'\{%s\}'%index
        # prepend 'name' to the group name to avoid error in case it is not letter
        group_regex = r'(?P<name%s>\w+)'%index
        newPattern = re.sub(index_string, group_regex, newPattern)

    matches = re.match(newPattern,string)
    matches_dict = matches.groupdict()

    # correct the gourp name and output it
    out = {}
    for key in matches_dict:
        out_key = key[4:]
        out[out_key] = matches_dict[key]

    return out

    