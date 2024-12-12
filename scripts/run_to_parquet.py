#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
from os import path
import concurrent.futures
import subprocess
import uproot as up
import pandas as pd
import to_parquet

if len(argv) > 1:
    in_path=argv[1]
else:
    in_path = "."

if not path.isdir(in_path):
    print(f'Error: {in_path} is not a directory')
    sys.exit(1)

files = glob(f'{in_path}/*.root')
print(f'Running files: {files}')

with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(to_parquet.process_file, files)

for f in files:
    print(f)
    with up.open(f) as file:
        tree = file['T']
        keys = tree.keys()
        df = tree.arrays(keys, library='pd')
        df.drop(columns=['truth_vec_cpt'],inplace=True,axis=1)
        df = df[df.matched_TtoR == True]
        o_name = path.basename(f)[:-5]+'.parquet'
        df.to_parquet(o_name,engine='fastparquet')

# print('a0')
# subprocess.run(['/home/davidstewart/jet_and_thermal/bin/jet_and_thermal', './brick_0_maxT_20_nEv_5K_pT_10_11.root', 'jetmatch_brick_0_maxT_20_nEv_5K_pT_10_11.root'])
# print('a1')
