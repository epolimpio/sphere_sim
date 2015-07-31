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