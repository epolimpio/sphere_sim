from myutils import readConfigFile
import metadata_lib

parameters = readConfigFile('parameters.ini')
print(metadata_lib.findFilesWithParameters('./data/metadata.dat',parameters))