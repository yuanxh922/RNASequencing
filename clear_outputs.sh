#!/bin/bash
#
# Reduce the file size of Jupyter Notebooks in batches.
# Usage: ./clear_outputs.sh <folder>
#
# Xiaohong Yuan
# yuanxh922@gmail.com
#

echo "Processing all '*.ipynb' files under the folder $1"

# save and change IFS
OLDIFS=$IFS
IFS=$'\n'

NOTEBOOKS=($(find $1 -type f -name '*.ipynb' | grep -v ipynb_checkpoints))

# restore it
IFS=$OLDIFS

tLen=${#NOTEBOOKS[@]}
for (( i=0; i<${tLen}; i++ )); do
	NOTEBOOK="${NOTEBOOKS[$i]}"
	echo
	echo "Processing $NOTEBOOK ..."
	jupyter nbconvert  --ClearOutputPreprocessor.enabled=True --inplace "$NOTEBOOK"
done
