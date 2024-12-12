#!/usr/bin/bash

# Check if the directory gz_files exists
if [ ! -d "splitting_files" ]; then
    echo "Directory splitting_files does not exist."
    exit 1
fi

# Check if the file namestems.txt exists
if [ ! -f "namestems.txt" ]; then
    echo "File namestems.txt does not exist."
    exit 1
fi

# Check if the directory splitting_files exists, if not create it
if [ ! -d "matched_jets" ]; then
    mkdir matched_jets
fi

exe=/home/davidstewart/jet_and_thermal/bin/jet_and_thermal
dir=`pwd -P`
# Read the file namestems.txt line by line
while read -r line
do
    inname=${dir}/splitting_files/${line}.root
    outname=${dir}/matched_jets/${line}.root
    ls $inname
    # Run the command for each line in namestems.txt
    ${exe} ${inname} ${outname}&
done < "namestems.txt"