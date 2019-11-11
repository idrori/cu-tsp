""" 
  train_rnn_coords.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script trains our proposed RNN encoder-decoder network for
  alpha carbon 3D coordinate prediction.

"""

import numpy as np
from numpy import array
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

from core import *
sys.path.append('..')
from utils.rgn_utils import *

# Returns our proposed RNN encoder-decoder model for 3D coordinate prediction
def get_coords_model():
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
  bigru_3 = Bidirectional(CuDNNGRU(units=100, return_sequences=True))(pssm)
  output_3 = Bidirectional(CuDNNGRU(units=64, return_sequences=True))(bigru_3)

  # Concatenation
  concat = concatenate([output_1, output_2, output_3], axis=-1)

  # 3D Coordinate Prediction
  coords = TimeDistributed(Dense(3), name='3d_output')(concat)
  #dcalpha = Lambda(pairwise_distance, name='3d_dcalpha_output')(coords)

  return Model(inputs=[aa_seq,ss_seq,pssm], outputs=coords)

# Trains a given model on folds (fold_start - fold_end)
def train_model(model, fold_start, fold_end):
  # Free up memory
  gc.collect()
  
  # Load data
  training_data = load_data(fold_start, fold_end, coords_flag=True)
  
  # Preprocess data
  X_train, train_coords = preprocess_data(training_data, 'coords')
  
  # Train
  for i in tqdm(range(len(X_train))):
    x = [np.array(X_train[i][0])[np.newaxis,:], np.array(X_train[i][1])[np.newaxis,:], X_train[i][2][np.newaxis,:]]
    y = train_coords[i][np.newaxis,:]
    if (i % 50 == 0): 
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=1)
    else: 
      history = model.fit(x, y, batch_size=BATCH_SIZE, epochs=N_EPOCHS, verbose=0)
      
  # Free up memory
  del training_data
  del X_train
  del train_coords
    
  return model

if __name__ == '__main__':
  start_time = datetime.now()

  recurrent_model = get_coords_model()
  recurrent_model.summary()
  rmsprop = keras.optimizers.RMSprop(lr=0.0002, rho=0.9, epsilon=None, decay=0.0, clipnorm=1.)
  recurrent_model.compile(optimizer=rmsprop, loss='mse', metrics=['mae'])

  #fold_pairs = [(1,1)]
  fold_pairs = [(10*(i-1)+1, 10*i) for i in range(1,11)]
  fold_pairs.append((101,105))

  for fold_pair in fold_pairs:
    fold_start, fold_end = fold_pair
    print('Training using folds: {}-{}...'.format(fold_start, fold_end))
    recurrent_model = train_model(recurrent_model, fold_start, fold_end)

  print('Training completed.')

  # Save trained model and weights
  try:
    rnn_json = recurrent_model.to_json()
    with open(os.path.join(COORDS_PATH, 'rnn_coords.json'), 'w') as fh:
      fh.write(rnn_json)
      recurrent_model.save_weights(os.path.join(COORDS_PATH, 'weights_coords.h5'))
      print('Model saved.')
  except:
    print('Error encountered while saving model.')

  end_time = datetime.now()
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Training Time: {}h {}m {}s'.format(hours, minutes, seconds))
