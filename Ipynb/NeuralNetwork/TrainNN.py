import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
import os
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras import mixed_precision
from tensorflow.keras.callbacks import CSVLogger
from tensorflow import distribute

from   tensorflow.python.client import device_lib
device_lib.list_local_devices()

import sys
sys.path.append('../')
from PlotAndJSON import PlotAndJSON

class TrainNN:
    '''
    - Trains a Random Forest on the data given
    - Saves to input directory:
    - model in model.pkl
     - json of features in features.json
     - json of model options in options.json
    '''
    PIECEWISE_FILE = '../IP_exp_decay/ExpPiecewiseData.json'

    mod_opt_def = {
        'optimizer' : 'adam',
        'loss' : 'mse',
        'epochs' : 4,
    }

    def __init__(self, data, odir, model_options={},
                x_keys=['reco_pt','rho_bkg_thermalandjet', 'reco_area'],
                test_size=0.2,
                epochs=4,
                # kFolds = 5, # not implemented yet
                # run_n_kFolds = 1, # not implemented yet
                random_state = 1,
                mini_NN=False,
                batch_size=1024
                ):
        ''' 
        Other inputs
        test_size=0.2, 
        n_estimators,  max_depth, max_features
        '''
        policy = mixed_precision.Policy('mixed_float16')
        mixed_precision.set_global_policy(policy)

        model_options['epochs'] = epochs
        model_options['test_size'] = test_size
        model_options['batch_size'] = batch_size

        odir_W = odir+'_W'

        # rhoA_resid = 'rhoA_resid_truthconst' if (y_key == 'reco_pt_truthconst') else 'rhoA_resid_truthfull'
        rhoA_resid = "rhoA_15to50"

        #update with default model options
        for k in TrainNN.mod_opt_def.keys():
            if (k not in model_options) and (TrainNN.mod_opt_def[k] is not None):
                model_options[k] = TrainNN.mod_opt_def[k]
        model_options['n_inputs'] = len(x_keys)

        for out in [odir, odir_W]:
            if not os.path.isdir(out):
                try:
                    os.makedirs(out)
                except:
                    exit(f'Could not create directory {out}')

            with open(os.path.join(out, 'model_options.json'), 'w') as f:
                json.dump(model_options, f)

            with open(os.path.join(out, 'features.json'), 'w') as f:
                json.dump(x_keys, f)
            
            with open(os.path.join(out, 'features.json'), 'w') as f:
                json.dump(x_keys, f)


            # with open(os.path.join(out, 'y_key.json'), 'w') as f:
                # json.dump([y_key,], f)


        y = data['truth_pt']
        X = data.drop(columns=['truth_pt'])

        if test_size == 0:
            X_train = X
            y_train = y
        else:
            X_train, X_test, y_train, y_test = train_test_split(X, y, 
                test_size=test_size, random_state=random_state)

        strategy = distribute.MirroredStrategy()

        with strategy.scope():
            model = Sequential()

            n_input = len(x_keys)
            print(f'input({n_input}): {x_keys}')

            if not mini_NN:
                model = Sequential()
                model.add(Dense(100, input_shape=(n_input,), activation='relu'))
                model.add(Dense(50, activation='relu'))
                model.add(Dense(50, activation='relu'))
            else:
                model = Sequential()
                model.add(Dense(10, input_shape=(n_input,), activation='relu'))
                model_options['epochs'] = epochs

            model.add(Dense(1))
            model.compile(optimizer=model_options['optimizer'], loss=model_options['loss'])

        logger_name = os.path.join(out, 'training_log.csv')
        csv_logger = CSVLogger(logger_name,append=True)
        model.fit(X_train[x_keys], y=y_train, epochs=model_options['epochs'],batch_size=batch_size,callbacks=[csv_logger])

        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        with open(os.path.join(out, 'TRAINING_STATUS.txt'), 'w') as f:
            f.write(f'Finished training at {current_time}\n')

        # with the update of keras, cannot now pickel the new keras models 
        print('FIXME HERE')
        model.save(os.path.join(odir_W, 'model.keras'))
        pickle.dump(model, open(os.path.join(odir_W, 'model.pkl'), 'wb'))

        # test and make output figures
        if test_size == 0:
            return
        else:
            y_pred = model.predict(X_test[x_keys])
            y_pred = y_pred.reshape(-1)
            X_test = X_test.join(y_test)
            PlotAndJSON(X_test, y_pred, X_test['truth_pt'], odir, 
                        odir, x_keys)
            PlotAndJSON(X_test, y_pred, X_test['truth_pt'], odir_W, 
                        odir_W, x_keys,
                        weightfn=self.PIECEWISE_FILE)

