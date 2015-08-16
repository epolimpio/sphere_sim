#!/bin/bash
echo Insert the file to compile followed by [ENTER]:
read filename
corr_filename=${filename%%.*}
f2py3 ${filename} -m ${corr_filename} -h ${corr_filename}.pyf --overwrite-signature
f2py3 -c ${corr_filename}.pyf ${filename}