""" 
  model_4_3.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script trains a version of our proposed RNN encoder-decoder 
  network that produced the best results for alpha-carbon distance
  matrix prediction.

  Note: This model currently does not have a model variant with MSA
  feature support.

"""

import numpy as np
from numpy import array
import argparse
import pickle
import gc
from tqdm import tqdm
import os, sys
from datetime import datetime
import matplotlib.pyplot as plt
import tensorflow as tf
import keras
from keras import backend as K
from keras.models import Model, Input, save_model, model_from_json
from keras.layers import dot, concatenate
from keras.layers import Dense, GRU, LSTM
from keras.layers import CuDNNGRU, CuDNNLSTM
from keras.layers import Lambda, Multiply
from keras.layers import Embedding, TimeDistributed, Bidirectional
from keras.layers import Average
from keras.backend.tensorflow_backend import set_session

from core import *
sys.path.append('..')
from utils.rgn_utils import *

MJ_FILE = '/home/arjunsrivatsa/results/rnn_torsion_dcalpha.json'
MW_FILE = '/home/arjunsrivatsa/results/weights_torsion_dcalpha.h5'

def load_model():
  with open(MJ_FILE, 'r') as json_file:
    loaded_model_json = json_file.read()

  custom_objects = {'tf': tf, 'reduce_l2_norm': reduce_l2_norm,
                    'dihedral_to_point': dihedral_to_point,
                    'point_to_coordinate': point_to_coordinate,
                    'pairwise_distance': pairwise_distance,
                    'theta': BOND_ANGLES,
                    'NUM_DIHEDRALS': NUM_DIHEDRALS,
                    'NUM_DIMENSIONS': NUM_DIMENSIONS,
                    'collections': collections}
  loaded_model = model_from_json(loaded_model_json, custom_objects)
  loaded_model.load_weights(MW_FILE)
  
  print('Model successfully loaded.')

  return loaded_model

# Returns the default model
# Note: This model is also used for Frobenius-norm-based loss variant
def get_model(log_flag=False, frobenius_flag=False):
  # AA Sequence
  aa_seq = Input(shape=(None,))
  embed_1 = Embedding(22,64)(aa_seq)

  # Q8 Sequence
  ss_seq = Input(shape=(None,))
  embed_2 = Embedding(9,64)(ss_seq)

  # PSSM
  pssm = Input(shape=(None,21))

  # Concatenation of Input Embedding
  merged_input = concatenate([embed_1, embed_2, pssm], axis=2)
  
  # Branched Bidirectional Recurrent Transform
  bigru_1 = Bidirectional(CuDNNGRU(units=1024, return_sequences=True))(merged_input)
  bigru_2 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_1)
  bigru_3 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_2)
  bigru_4 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_3)
  dense_3_1 = TimeDistributed(Dense(64, activation = 'relu'))(bigru_3)
  concat = concatenate([bigru_1, bigru_2, bigru_3, bigru_4, dense_3_1], axis=-1)

  # 3D Coordinate Prediction
  a2coords = TimeDistributed(Dense(500), name='alt2')(concat)
  new_coords = TimeDistributed(Dense(3), name='conversion_coords')(a2coords)

  # Distance Matrix Prediction
  dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(new_coords)

  # Log Distance Matrix
  if log_flag:
    log_dcalpha = Lambda(log_mat, name='log_dcalpha')(dcalpha)
    outputs = [log_dcalpha, log_dcalpha]

  # Frobenius-norm-based Loss Training
  elif frobenius_flag:
    outputs = [dcalpha, dcalpha]

  # Log-ratio Loss Training
  else:
    outputs = dcalpha

  return Model(inputs=[aa_seq, ss_seq, pssm], outputs=outputs)

# Returns the model with log-dcalpha prediction
def get_model_log_dcalpha():
  # AA Sequence
  aa_seq = Input(shape=(None,))
  embed_1 = Embedding(22,64)(aa_seq)

  # Q8 Sequence
  ss_seq = Input(shape=(None,))
  embed_2 = Embedding(9,64)(ss_seq)

  # PSSM
  pssm = Input(shape=(None,21))

  # Concatenation of Input Embedding
  merged_input = concatenate([embed_1, embed_2, pssm], axis=2)
  
  # Branched Bidirectional Recurrent Transform
  bigru_1 = Bidirectional(CuDNNGRU(units=1024, return_sequences=True))(merged_input)
  bigru_2 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_1)
  bigru_3 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_2)
  bigru_4 = Bidirectional(CuDNNGRU(units=512, return_sequences=True))(bigru_3)
  dense_3_1 = TimeDistributed(Dense(64, activation = 'relu'))(bigru_3)
  concat = concatenate([bigru_1, bigru_2, bigru_3, bigru_4, dense_3_1], axis=-1)

  # 3D Coordinate Prediction
  a2coords = TimeDistributed(Dense(500), name='alt2')(concat)
  new_coords = TimeDistributed(Dense(3), name='conversion_coords')(a2coords)

  # Distance Matrix Prediction
  dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(new_coords)
  
  # Log Distance Matrix
  log_dcalpha = Lambda(log_mat, name='log_dcalpha')(dcalpha)

  # TODO: Probably better to change this to single output for load_eval_model
  return Model(inputs=[aa_seq, ss_seq, pssm], outputs=[log_dcalpha, log_dcalpha])
  
# Trains a given model on folds (fold_start - fold_end)
def train_model(model, fold_start, fold_end, log_flag=False, frobenius_flag=False,
                ckpt_flag=False):
  # Free up memory
  gc.collect()
  
  # Load data
  training_data = load_data(fold_start, fold_end, radians_flag=True)
  
  # Preprocess data
  X_train, combined = preprocess_data(training_data, 'dcalpha')
  
  # Train
  for i in tqdm(range(len(X_train))):
    x = [np.array(X_train[i][0])[np.newaxis,:], np.array(X_train[i][1])[np.newaxis,:], 
         X_train[i][2][np.newaxis,:]]
    
    if log_flag:
      m = combined[i][np.newaxis,:] + EPSILON
      n = np.log(m)
      y = [n,n]
    elif frobenius_flag:
      m = combined[i][np.newaxis,:]
      y = [m,m]
    else:
      y = combined[i][np.newaxis,:]
    
    if (i % 50 == 0): 
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=1)
    else: 
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=0)
  
  # Save checkpoint
  if ckpt_flag:
    try:
      ckpt_path = os.path.join(RESULTS_PATH, 'model_4_3_{}.h5'.format(fold_end))
      model.save(ckpt_path)
      print('Model trained up to folds {} saved.'.format(fold_end))
    except:
      print('[Error] Checkpoint saving failed, ignoring.')
  
  # Free up memory
  del training_data
  del X_train
  del combined
    
  return model

# Trains the model according to the given specifications
def main(**kwargs):
  log_flag = kwargs['log']
  frobenius_flag = kwargs['frobenius']
  ckpt_flag = kwargs['ckpt']

  # Check if given options are valid
  if log_flag and frobenius_flag:
    raise RuntimeError('[Error] Conflicting options: "log" and "frobenius."')

  # Check path to save training results
  if not os.path.exists(RESULTS_PATH):
    try:
      mkdir(RESULTS_PATH)
    except:
      raise RuntimeError('[Error] Results path setting is incorrect.')

  # Define optimizer
  rmsprop = keras.optimizers.RMSprop(lr=0.00015, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)

  # Fetch model to train
  ## Log-dcalpha model
  if log_flag:
    model = get_model(log_flag=True)
    model.summary()
    rmsprop = keras.optimizers.RMSprop(lr=0.0002, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)
    model.compile(optimizer=rmsprop, loss=['mse', matrix_loss], loss_weights=[0.5, 1],
                  metrics=['mae', 'mse', matrix_loss])
  
  ## Frobenius-norm-based loss
  elif frobenius_flag:
    model = get_model(frobenius_flag=True)
    model.summary()
    rmsprop = keras.optimizers.RMSprop(lr=0.0002, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)
    model.compile(optimizer=rmsprop, loss=[matrix_loss, 'mae'], loss_weights=[1, 0.05],
                  metrics=['mae'])
  
  ## Log-ratio loss (Best Results)
  else:
    model = get_model()
    model.summary()
    rmsprop = keras.optimizers.RMSprop(lr=0.00015, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)
    model.compile(optimizer=rmsprop, loss=energetic_loss, metrics=['mae'])

  fold_pairs = [(1,1)]
  #fold_pairs = [(10*(i-1)+1, 10*i) for i in range(1,11)]
  #fold_pairs.append((101,104))

  start_time = datetime.now()

  for fold_pair in fold_pairs:
    fold_start, fold_end = fold_pair
    print('Training using folds: {}-{}...'.format(fold_start, fold_end))
    model = train_model(model, fold_start, fold_end, log_flag, frobenius_flag,
            ckpt_flag)

  print('Training completed.')
  end_time = datetime.now()

  # Save trained model and weights
  model.save(MODEL3_PATH)
  print('Model and weights saved.')

  # Report total training time
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Training Time: {}h {}m {}s'.format(hours, minutes, seconds))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--log', help='Option to train with log-dcalpha model', action='store_true',
                      default=False)
  parser.add_argument('--frobenius', help='Option to train with Frobenius-norm-based loss', 
                      action='store_true', default=False)
  parser.add_argument('--ckpt', help='Option to save model after every fold_pair training', 
                      action='store_true', default=False)
  args = parser.parse_args()

  main(**vars(args))
