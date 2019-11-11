""" 
  model_4_2.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script trains a version of our proposed RNN encoder-decoder 
  network that predicts torsion angles and generates alpha-carbon
  distance matrix predictions.

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
from keras.layers import Dense, GRU, LSTM, Conv2D
from keras.layers import CuDNNGRU, CuDNNLSTM
from keras.layers import Lambda, Multiply
from keras.layers import Embedding, TimeDistributed, Bidirectional
from keras.layers import Average

from core import *
sys.path.append('..')
from utils.rgn_utils import *
from sklearn import manifold

# Returns the model without MSA feature support
def get_model():
  # AA Sequence
  aa_seq = Input(shape=(None,))
  embed_1 = Embedding(22,32)(aa_seq)
  lstm_1 = CuDNNLSTM(units=75, return_sequences=True)(embed_1)
  bilstm_1 = Bidirectional(CuDNNLSTM(units=50, return_sequences=True))(lstm_1)
  bigru_1 = Bidirectional(CuDNNGRU(units=50, return_sequences=True))(bilstm_1)
  output_1 = CuDNNLSTM(units=64, return_sequences=True)(bigru_1)

  # Q8 Sequence
  ss_seq = Input(shape=(None,))
  embed_2 = Embedding(9,32)(ss_seq)
  bigru_2_1 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(embed_2)
  bigru_2_2 = Bidirectional(CuDNNGRU(units=50, return_sequences=True))(bigru_2_1)
  output_2 = Dense(units=64, activation='relu')(bigru_2_2)

  # PSSM
  pssm = Input(shape=(None,21))
  bigru_3_1 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(pssm)
  bigru_3_2 = Bidirectional(CuDNNGRU(units=64, return_sequences=True))(bigru_3_1)
  output_3 = Dense(units=64, activation='relu')(bigru_3_2)

  # Concatenation
  concat = concatenate([output_1, output_2, output_3], axis=-1)

  # Torsion Angle and Distance Matrix Prediction
  angles = TimeDistributed(Dense(3), name='3d_output')(concat)
  angles_dihedral = Lambda(slice, name = 'dihedrals')(angles)
  pts = Lambda(dihedral_to_point, name = 'pts_from_angle')(angles)
  cc = Lambda(point_to_coordinate, name = 'coords_from_pts')(pts)
  dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(cc)

  return Model(inputs=[aa_seq, ss_seq, pssm], outputs=[angles_dihedral, dcalpha])

# Returns the model with MSA feature support
def get_msa_model():
  # AA Sequence
  aa_seq = Input(shape=(None,))
  embed_1 = Embedding(22,32)(aa_seq)
  lstm_1 = CuDNNLSTM(units=75, return_sequences=True)(embed_1)
  bilstm_1 = Bidirectional(CuDNNLSTM(units=50, return_sequences=True))(lstm_1)
  bigru_1 = Bidirectional(CuDNNGRU(units=50, return_sequences=True))(bilstm_1)
  output_1 = CuDNNLSTM(units=64, return_sequences=True)(bigru_1)

  # Q8 Sequence
  ss_seq = Input(shape=(None,))
  embed_2 = Embedding(9,32)(ss_seq)
  bigru_2_1 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(embed_2)
  bigru_2_2 = Bidirectional(CuDNNGRU(units=50, return_sequences=True))(bigru_2_1)
  output_2 = Dense(units=64, activation='relu')(bigru_2_2)

  # PSSM
  pssm = Input(shape=(None,21))
  bigru_3_1 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(pssm)
  bigru_3_2 = Bidirectional(CuDNNGRU(units=64, return_sequences=True))(bigru_3_1)
  output_3 = Dense(units=64, activation='relu')(bigru_3_2)

  # MSA Covariance One Hot
  msa_onehot = Input(shape=(None,None,1))
  conv_4_1 = Conv2D(60, 5, data_format='channels_last', padding='same')(msa_onehot)
  conv_4_2 = Conv2D(30, 5, data_format='channels_last', padding='same')(conv_4_1)
  sliced_4 = Lambda(slice_msa, name='sliced_msa_onehot')(conv_4_2)
  bigru_4 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(sliced_4)
  output_4 = Bidirectional(CuDNNGRU(units=32, return_sequences=True))(bigru_4)

  # MSA Covariance Embedding
  msa_embed = Input(shape=(None,None,1))
  conv_5_1 = Conv2D(60, 5, data_format='channels_last', padding='same')(msa_embed)
  conv_5_2 = Conv2D(30, 5, data_format='channels_last', padding='same')(conv_5_1)
  sliced_5 = Lambda(slice_msa, name='sliced_msa_embed')(conv_5_2)
  bigru_5 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(sliced_5)
  output_5 = Bidirectional(CuDNNGRU(units=32, return_sequences=True))(bigru_5)

  # Concatenation
  concat = concatenate([output_1, output_2, output_3, output_4, output_5], axis=-1)

  # Torsion Angle and Distance Matrix Prediction
  angles = TimeDistributed(Dense(3), name='3d_output')(concat)
  angles_dihedral = Lambda(slice, name = 'dihedrals')(angles)
  pts = Lambda(dihedral_to_point, name = 'pts_from_angle')(angles)
  cc = Lambda(point_to_coordinate, name = 'coords_from_pts')(pts)
  dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(cc)

  return Model(inputs=[aa_seq, ss_seq, pssm, msa_onehot, msa_embed], outputs=[angles_dihedral, dcalpha])

# Returns the model without Q8 and MSA feature support
def get_no_q8_model():
  # AA Sequence
  aa_seq = Input(shape=(None,))
  embed_1 = Embedding(22,32)(aa_seq)
  lstm_1 = CuDNNLSTM(units=75, return_sequences=True)(embed_1)
  bilstm_1 = Bidirectional(CuDNNLSTM(units=50, return_sequences=True))(lstm_1)
  bigru_1 = Bidirectional(CuDNNGRU(units=50, return_sequences=True))(bilstm_1)
  output_1 = CuDNNLSTM(units=64, return_sequences=True)(bigru_1)

  # PSSM
  pssm = Input(shape=(None,21))
  bigru_2_1 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(pssm)
  bigru_2_2 = Bidirectional(CuDNNGRU(units=64, return_sequences=True))(bigru_2_1)
  output_2 = Dense(units=64, activation='relu')(bigru_2_2)

  # Concatenation
  concat = concatenate([output_1, output_2], axis=-1)

  # Torsion Angle and Distance Matrix Prediction
  angles = TimeDistributed(Dense(3), name='3d_output')(concat)
  angles_dihedral = Lambda(slice, name = 'dihedrals')(angles)
  pts = Lambda(dihedral_to_point, name = 'pts_from_angle')(angles)
  cc = Lambda(point_to_coordinate, name = 'coords_from_pts')(pts)
  dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(cc)

  return Model(inputs=[aa_seq, pssm], outputs=[angles_dihedral, dcalpha])

# Trains a given model on folds (fold_start - fold_end)
def train_model(model, fold_start, fold_end, msa_flag=False, no_q8_flag=False,
                ckpt_flag=False):
  # Free up memory
  gc.collect()

  # Load data
  training_data = load_data(fold_start, fold_end, msa_flag=msa_flag,
                            no_q8_flag=no_q8_flag, radians_flag=True)

  # Preprocess data
  X_train, combined = preprocess_data(training_data, 'combined', msa_flag=msa_flag)

  # Train
  for i in tqdm(range(len(X_train))):
    if msa_flag:
      x = [np.array(X_train[i][0])[np.newaxis,:], np.array(X_train[i][1])[np.newaxis,:], 
           X_train[i][2][np.newaxis,:], np.array(X_train[i][3])[np.newaxis,:], 
           np.array(X_train[i][4])[np.newaxis,:]]
    elif no_q8_flag:
      x = [np.array(X_train[i][0])[np.newaxis,:], X_train[i][2][np.newaxis,:]]
    else:
      x = [np.array(X_train[i][0])[np.newaxis,:], np.array(X_train[i][1])[np.newaxis,:], 
           X_train[i][2][np.newaxis,:]]
    
    y = [combined[i][0][np.newaxis,:], combined[i][1][np.newaxis,:]]

    if (i % 50 == 0):
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=1)
    else:
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=0)

  # Save checkpoint
  if ckpt_flag:
    try:
      ckpt_path = os.path.join(RESULTS_PATH, 'model_4_2_{}.h5'.format(fold_end))
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
  msa_flag = kwargs['msa']
  no_q8_flag = kwargs['no_q8']
  ckpt_flag = kwargs['ckpt']

  # Check if given options are valid
  if msa_flag and no_q8_flag:
    raise RuntimeError('[Error] Conflicting options: "msa" and "no_q8."')

  # Check path to save training results
  if not os.path.exists(RESULTS_PATH):
    try:
      mkdir(RESULTS_PATH)
    except:
      raise RuntimeError('[Error] Results path setting is incorrect.')

  # Fetch model to train
  if msa_flag:
    model = get_msa_model()
  elif no_q8_flag:
    model = get_no_q8_model()
  else:
    model = get_model()

  model.summary()
  rmsprop = keras.optimizers.RMSprop(lr=0.0002, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)
  model.compile(optimizer=rmsprop, loss='mse', metrics=['mae'])

  fold_pairs = [(1,1)]
  #fold_pairs = [(10*(i-1)+1, 10*i) for i in range(1,11)]
  #fold_pairs.append((101,104))

  start_time = datetime.now()

  for fold_pair in fold_pairs:
    fold_start, fold_end = fold_pair
    print('Training using folds: {}-{}...'.format(fold_start, fold_end))
    model = train_model(model, fold_start, fold_end, msa_flag, no_q8_flag, ckpt_flag)

  print('Training completed.')
  end_time = datetime.now()

  # Save trained model and weights
  model.save(MODEL2_PATH)
  print('Model and weights saved.')

  # Report total training time
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Training Time: {}h {}m {}s'.format(hours, minutes, seconds))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--msa', help='Option to train on MSA features', action='store_true',
                      default=False)
  parser.add_argument('--no_q8', help='Option to not use Q8 secondary structure sequences', 
                      action='store_true', default=False)
  parser.add_argument('--ckpt', help='Option to save model after every fold_pair training', 
                      action='store_true', default=False)
  args = parser.parse_args()

  main(**vars(args))
