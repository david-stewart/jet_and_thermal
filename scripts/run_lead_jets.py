#!/home/davidstewart/penv-ML/bin/python3

from glob import glob
from sys import argv
import os
import concurrent.futures
import subprocess

exe = '/home/davidstewart/jet_and_thermal/bin/IP_lead_jets'

if not os.path.isdir("lead_jets"):
    os.mkdir("lead_jets")

def process_file(file):
    f_out = f'lead_jets/{os.path.basename(file)}'
    todo =[exe, file, '-1', f_out]
    print(todo)
    subprocess.run(todo)

files = [f for f in glob(f'./*.root')]
print(f'Running files: {files}')

if True:
  with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(process_file, files)

# for file in files:
    # f_out = f'lead_jets/{os.path.basename(file)}'
    # print([exe, file, '-1', f_out])
    # subprocess.run([exe, file, '-1', f_out])
    # break
