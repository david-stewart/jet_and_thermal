#!/home/davidstewart/penv-ML/bin/python3

import gunziptree as ggt
from glob import glob
from sys import argv
from os import path
import concurrent.futures
import subprocess

def process_file(file):
    ggt.gunziptree(file)
    print(f'Finished {file}')


if len(argv) > 1:
    in_path=argv[1]
else:
    in_path = "."

if not path.isdir(in_path):
    print(f'Error: {in_path} is not a directory')
    sys.exit(1)

files = glob(f'{in_path}/*.dat.gz')
print(f'{in_path}/*.dat.gz')
print(f'Running files: {files}')

with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_file, files)
