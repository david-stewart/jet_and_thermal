from glob import glob
import os
import json
import pandas as pd
from operator import attrgetter

class JetMatchFile:
    def __init__(self, jet_match_file):
        self.filename = jet_match_file
        self.path = os.path.dirname(jet_match_file)
        self.base = os.path.basename(jet_match_file)
        self.stem, self.extension = os.path.splitext(self.base)

        self.items = self.stem.split('_')

        print(f'file {jet_match_file} ---  items: {self.items}')
        self.brick = float(self.items[1])
        self.maxT  = float(self.items[3])
        self.nEv   = int(1000*self.items[5][:-1])
        self.pT0   = float(self.items[7])
        self.pT1   = float(self.items[8])

        self.is_lbt = 'lbt_brick' in self.path
        self.brick_len = self.brick

    def tag(self):
        tag = 'lbt' if self.is_lbt else 'matter'
        return f'{tag}_{self.brick}_maxT_{self.maxT}_pthat_{self.pT0}_{self.pT1}'
    def odir(self, stem):
        prefix = 'lbt' if self.is_lbt else 'matter'
        return f'{prefix}/fit_{stem}/{self.tag()}'

class JetMatchesSet:
    def __init__(self, matched_dirs):
        if not isinstance(matched_dirs, list):
            matched_dirs = [matched_dirs]

        self.files = []
        for idir in matched_dirs:
            for f in glob(idir):
                self.files.append(JetMatchFile(f))

    def __iter__(self):
        return iter(self.files)

    def sort(self, key):
        self.files.sort(key=lambda x: getattr(x, key))

class ML_Result:
    def __init__(self, in_file):
        self.filename = in_file
        with open(in_file, 'r') as f:
            self.data = json.load(f)

        self.mean_rhoA = self.data['mean_delta_rhoA_all']
        self.mean_ML   = self.data['mean_delta_ML_all']
        self.std_rhoA  = self.data['std_delta_rhoA_all']
        self.std_ML    = self.data['std_delta_ML_all']

        self.items = self.filename.split('/')[-2].split('_')
        if self.items[0] == 'hydro':
            self.brick = 0.
            self.maxT = 0.
            self.pT0 = float(self.items[2])
            self.pT1 = float(self.items[3])
        else:
            self.brick = float(self.items[1])
            self.maxT  = float(self.items[3])
            self.pT0   = float(self.items[5])
            self.pT1   = float(self.items[6])

        self.is_hydro = 'hydro' in self.filename
        self.is_lbt = 'lbt' in self.filename
        self.pthat = self.pT0
        self.brick_len = self.brick


        if 'rhoAonly' in self.filename:
            self.features = 'rhoAonly'
        elif 'angularity' in self.filename:
            self.features = 'angularity'
        elif 'reco_all' in self.filename:
            self.features = 'reco_all'
        elif 'nconsts' in self.filename:
            self.features = 'nconsts'
        elif 'reco_cpt' in self.filename:
            self.features = 'reco_cpt'

        x_keys = None
        
    def dict(self):
        d = {key:value for key, value in self.__dict__.items() 
                if not key[:2] == ('__') 
                and not callable(key)
                and not key == 'data'
            }
        for key, val in self.data.items():
            d[key] = val
        return d

class MLResultsSets:
    def __init__(self, stems=['matter','lbt','hydro']):
        self.results = []
        for stem in stems:
            for path in glob(f'./{stem}/**/rss_score.json', recursive=True):
                self.results.append(ML_Result(path))

    def __iter__(self):
        return iter(self.results)

    def sort(self, *attr):
        self.results.sort(key=attrgetter(*attr))

    def dataframe(self):
        print(self.results)
        return pd.DataFrame([result.dict() for result in self.results])

    def trim_brick_gte(self, val):
        self.results = [result for result in self.results if result.brick < val]
