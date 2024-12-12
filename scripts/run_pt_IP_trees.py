#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
from os import path
import os
import concurrent.futures
import subprocess

# run this script from insize the gz_files directory

if not os.path.isdir('./IP_trees'):
    os.mkdir('./IP_trees')

exe = '/home/davidstewart/jet_and_thermal/bin/pt_IP_trees'
def process_file(file):
    print(file)
    f_out = './IP_trees/'+path.basename(file)
    print('a2')
    subprocess.run([exe, file, f_out])

if len(argv) > 1:
    in_path=argv[1]
else:
    in_path = "."

files = glob('*.root')

if False:
  for file in files:
    print(file)
    f_out = './IP_trees/'+path.basename(file)
    print(f_out)
    subprocess.run([exe, file, f_out])

# print(f'Running files: {files}')
# print  (f'./JM_{path.basename(files[0])}')

if True: 
  with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_file, files)
