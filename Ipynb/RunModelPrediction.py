'''
This cell:
    read-in the Xsec_group_#.parquet files, correct, and write the short output locall
'''

# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ROOT as rt
import pickle
import os
from glob import glob
import json

# parameters
model_path = '../../NeuralNetwork/train/reco_cpt_W/'

# short class to get the model
class PickledModel:
    def __init__(self, model_path):
        self.model_path = model_path
        self.model = None
        self.parameters = None
        self.load_model()
        self.load_features()

    def load_model(self):
        with open(self.model_path + 'model.pkl', 'rb') as f:
            self.model = pickle.load(f)

    def load_features(self):
        with open(self.model_path + 'features.json', 'r') as f:
            self.features = json.load(f)

    def predict(self, x):
        return self.model.predict(x)

def runModelPrediction(model_path):
    # get the model for corrections
    ''' add a method to do the AB_method without ML '''
    if model_path == 'AB_method':
        model = 'AB_method'
    else:
        model = PickledModel(model_path)

    if not os.path.isdir('truth_reco_corr'):
        os.mkdir('truth_reco_corr')

    n_infiles = len(glob('../Xsec_groups/df_Xsec_*_clean.parquet'))
    cnt = 0

    for file in glob('../Xsec_groups/df_Xsec_*_clean.parquet'):
        print(f"Reading {cnt} of {n_infiles} : file {file}")
        cnt += 1


        i = int(file.split('_')[-2])
        outfile = f'truth_reco_corr/data_{i}_clean.parquet'
        
        if os.path.isfile(outfile):
            print(f"File {outfile} already exists. Skipping.")
            continue

        try:
            df = pd.read_parquet(file)
        except Exception as e:
            print(f"Error: {e}")
            break

        print(f'pre-match-cut size: {df.shape[0]}')
        df = df[(df.matched_TtoR==True)&(df.matched_IPtoT==True)]
        print(f'post-match-cut size: {df.shape[0]}')

        # print(corr)
        if model == 'AB_method':
            df['corr_pt'] = df['reco_pt']-df['rho_bkg_thermal']*df['reco_area']
            df = df[df['corr_pt']>0].copy()
        else:
            df['corr_pt'] = model.predict(df[model.features].values)

        print(df['corr_pt'].describe())
        df_short = df[['truth_pt', 'corr_pt', 'reco_pt', 'reco_area', 'rho_bkg_thermal']]
        df_short.to_parquet(outfile)
        
    print("Done.")
    with open('DONE_clean.txt', 'w') as f:
        f.write('done')