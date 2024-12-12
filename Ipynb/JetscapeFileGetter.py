#!/home/davidstewart/penv-ML/bin/python3
from glob import glob
import copy
import pandas as pd
import numpy as np

''' names are like:
 config/lbt_brick/gz_files/IP_match_trees/brick_0_maxT_20_nEv_5K_pT_10_11.root
 config/matter_brick/gz_files/IP_match_trees/brick_9_maxT_20_nEv_5K_pT_10_11.root
'''
class JetscapeFile:
    # update for add hydro files
    def __init__(self, file, do_print=False):
        try:
            self.is_lbt =    True if ('lbt_brick' in file) else False
            self.is_matter = True if ('matter_brick' in file) else False
            self.is_hydro = True if ('hydro' in file) else False
            items = file.split('/')[-1].split('.')[0].split('_')
            if self.is_hydro:
                self.length = -1.
                self.T = -1.
                self.nEvents = int(1000*int(items[-1][:-1]))
                self.pt0 = float(items[2])
                self.pt1 = float(items[3])
                self.file = file
            else:
                self.length = float(items[1])
                self.T = float(items[3])
                self.nEvents = int(1000*int(items[5][:-1]))
                self.pt0 = float(items[7])
                self.pt1 = float(items[8])
                self.file  = file
        except:
            print(f'Error parsing file: {file}')
            raise
        if do_print:
            print(self.dict)

    def prefix(self):
        if self.is_lbt:
            return 'lbt'
        elif self.is_matter:
            return 'matter'
        else:
            return 'hydro'

    def tag(self):
        if self.is_hydro:
            return f'hydro_pthat_{self.pt0}_{self.pt1}'
        else:
            tag = 'lbt' if self.is_lbt else 'matter'
            return f'{self.prefix()}_{self.length}_maxT_{self.T}_pthat_{self.pt0}_{self.pt1}'

    def odir(self, stem):
        return f'{self.prefix()}/fit_{stem}/{self.tag()}'
    
    def get_normed_stats(self, col):
        df = pd.read_parquet(self.file)
        vals = np.array(df[col])
        staterr = np.sqrt(vals)
        bins = np.array(df['bin'])
        normval = vals.sum()*(bins[1]-bins[0])
        vals = vals/normval
        errs = staterr/normval
        return bins, vals, errs

class JetscapeFileGetter:
    def __init__(self, config_path='../config', which='both',file_type='parquet',
      sort='length', reverse=False, cut=None): # both, lbt, or matter
        files = []
        if which == 'both' or which == 'lbt':
            files.extend(glob(f'{config_path}/lbt_brick/gz_files/IP_match_trees/*.{file_type}'))
        if which == 'both' or which == 'matter':
            files.extend(glob(f'{config_path}/matter_brick/gz_files/IP_match_trees/*.{file_type}'))
        if which == 'hydro':
            files.extend(glob(f'{config_path}/hydro/gz_files/IP_match_trees/*.{file_type}'))
        print(files)
        self.files = [JetscapeFile(f) for f in files]

        if cut:
            self.files = [f for f in self.files if getattr(f,cut[0]) == cut[1]]

        self.files.sort(key=lambda x: getattr(x, sort), reverse=reverse)

    def __iter__(self):
        return iter(self.files)

    def __call__(self, select='is_lbt',val=True, sort='length', reverse=False):
        rdat = JetscapeFileGetter(config_path='none')

        for f in self:
            if getattr(f, select) == val:
                rdat.files.append(f)

        rdat.files.sort(key=lambda x: getattr(x, sort), reverse=reverse)

        return rdat

    def __getitem__(self, index):
        return self.files[index]

    

def main():
    getter = JetscapeFileGetter()
    for file in getter:
        print(file)

if __name__ == "__main__":
    files = JetscapeFileGetter(config_path='./config',cut=('pt0',10),which='lbt')

    for x in files:
        print(x.file)