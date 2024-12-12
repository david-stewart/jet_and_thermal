#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
from os import path
import os
import concurrent.futures
import subprocess

# run this script from insize the gz_files directory

if not os.path.isdir('./IP_match_trees'):
    os.mkdir('./IP_match_trees')

exe = '/home/davidstewart/jet_and_thermal/bin/IP_match_trees'
def process_file(in_file):
    print(in_file)
    # in_file
    nevents='-1'
    out_file = './IP_match_trees/'+path.basename(in_file)
    in_hydro='/home/davidstewart/jetscape-docker/JETSCAPE/config/hydro/gz_files/hydro_31K.root'
    use_cpt='1'
    bkg_samples='20'
    seed='-100'
    subprocess.run([exe, in_file, nevents, out_file, in_hydro, use_cpt, bkg_samples, seed])

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
    subprocess.run([exe, file, f_out, ])

# print(f'Running files: {files}')
# print  (f'./JM_{path.basename(files[0])}')

if True: 
  with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_file, files)
