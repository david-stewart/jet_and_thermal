#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
import sys
from os import path
import concurrent.futures
import subprocess
import uproot as up
import pandas as pd

def process_file(f, tree_name):
    with up.open(f) as file:
        tree = file[tree_name]
        keys = tree.keys()
        if 'truth_vec_cpt' in keys:
            keys.remove('truth_vec_cpt')
        # print(keys)
        # return
        df = tree.arrays(keys, library='pd')
        o_name = path.basename(f)[:-5]+'.parquet'
        df.to_parquet(o_name,engine='fastparquet')

# if not path.isfile(argv[1]):
    # print(f'Error: {argv[1]} is not a file')
    # sys.exit(1)

if __name__=="__main__":
    if len(argv) < 3:
        tree_name='T'
    else:
        tree_name=argv[2]
    process_file(argv[1], tree_name)
