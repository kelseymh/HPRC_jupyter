#!/bin/bash
#
# Usage: source ~kelsey/jupyter/setup.sh
#
# Adds my Jupyter directory to PYTHONPATH so that `import` actions work.
#
# 20241204  Michael Kelsey
# 20250314  Modify 'mydir' for ComputeCanada (Cedar)
# 20250317  Use $BASH_ARGV to discover path to setup file

mydir=.
[ -n "$BASH_VERSION" ] && mydir=`dirname ${BASH_ARGV[0]}`
mydir=`realpath $mydir`

export PYTHONPATH=${PYTHONPATH}:${mydir}
export PATH=${PATH}:${mydir}
