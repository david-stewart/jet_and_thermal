#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
import os
import concurrent.futures
import subprocess

exe = '/home/davidstewart/jet_and_thermal/bin/jet_and_thermal'

if len(argv) > 1:
    n_embs = int(argv[1])
else:
    n_embs = 1

if not os.path.isdir("jet_matches"):
    os.mkdir("jet_matches")

def process_file(file):
    f_out = f'jet_matches/{os.path.basename(file)}'
    print((file, f_out))
    subprocess.run([exe, file, f_out, str(n_embs)])

files = [f for f in glob(f'./*.root')]
print(f'Running files: {files}')

if True:
  with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_file, files)

# for file in files:
    # f_out = f'jet_matches/{os.path.basename(file)}'
    # print(         [exe, file, f_out, n_embs])
    # subprocess.run([exe, file, f_out, str(n_embs)])
    # break
