#!/bin/bash

# Check if the QP_ROOT environment variable is set. 
if [[ -z ${QP_ROOT} ]]
then
  print "The QP_ROOT environment variable is not set."
  print "Please reload the quantum_package.rc file."
  exit -1
fi
 
# Get the absolute path of the current directory.
currdir=${PWD}

file_list="`ls ../func_ecmd/*.irp.f`"

# Make a symbolic link for all scripts to be used in the ${QP_ROOT}/scripts/
# directory. 
 for i in $file_list
 do 
  echo $i
  ln --symbolic  ${currdir}/$i  ${QP_ROOT}/src/functionals/
 done
