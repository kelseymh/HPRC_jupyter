#!/bin/bash
#
# Usage: source ~kelsey/jupyter/setup.sh
#
# Adds my Jupyter directory to PYTHONPATH so that `import` actions work.
#
# 20241204  Michael Kelsey

mydir=/scratch/user/kelsey/jupyter
export PYTHONPATH=${PYTHONPATH}:${mydir}
export PATH=${PATH}:${mydir}

# TODO: Should we load useful modules here as well?
