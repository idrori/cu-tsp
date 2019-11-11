"""
  core.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script contains the core functions and settings used for training
  and testing our proposed RNN encoder-decoder networks.

"""

import math
import numpy as np
from numpy import array
import pickle
import statistics
import gc
from tqdm import tqdm
import os
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

# CONSTANTS

AA_CODES = 'ACDEFGHIKLMNPQRSTVWXY'
Q8_CODES = 'GHITEBS-'

COORDS_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed
DCALPHA_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed
TORSION_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed
COORDS_DCALPHA_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed
TORSION_DCALPHA_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed
MDS_DCALPHA_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results' # Modify as needed

EPSILON = 1e-06
RAD = math.pi / 180
DEG = 180 / math.pi

# TRAINING

TRAINING_PATH = '/scratch/ass2186/train/training_100_{}.pkl' # Modify as needed

N_CHANNELS = 21
BATCH_SIZE = 1
N_EPOCHS = 6

# TESTING

#TEST_PATH = '/home/darshan/cu-tsp-dev/data/cuprotein_testing_casp13.pkl' # Modify as needed
#VALIDATION_PATH = '/home/darshan/cu-tsp-dev/data/cuprotein_testing_casp13.pkl' # Modify as needed
TEST_PATH = '/home/darshan/cu-tsp-dev/data/casp13_predictedq8.pkl' # Modify as needed
VALIDATION_PATH = '/home/darshan/cu-tsp-dev/data/casp13_predictedq8.pkl' # Modify as needed

RESULTS_PATH = '/home/darshan/cu-tsp-dev/model_4_branch/results'  # Modify as needed

CSV_OUT_PATH = os.path.join(RESULTS_PATH, 'predictions_casp13.csv')
CSV_DCALPHA_FORMAT = '{}_d_{}_{}'
CSV_PSI_FORMAT = '{}_psi_{}'
CSV_PHI_FORMAT = '{}_phi_{}'

PKL_OUT_PATH = os.path.join(RESULTS_PATH, 'predictions_casp13_log.pkl')

COORDS_MJ = os.path.join(COORDS_PATH, 'rnn_coords.json')
COORDS_MW = os.path.join(COORDS_PATH, 'weights_coords.h5')

DCALPHA_MJ = os.path.join(DCALPHA_PATH, 'rnn_dcalpha.json')
DCALPHA_MW = os.path.join(DCALPHA_PATH, 'weights_dcalpha.h5')

TORSION_MJ = os.path.join(TORSION_PATH, 'rnn_torsion.json')
TORSION_MW = os.path.join(TORSION_PATH, 'weights_torsion.h5')

COORDS_DCALPHA_MJ = os.path.join(COORDS_DCALPHA_PATH, 'rnn_coords_dcalpha.json')
COORDS_DCALPHA_MW = os.path.join(COORDS_DCALPHA_PATH, 'weights_coords_dcalpha.h5')

TORSION_DCALPHA_MJ = os.path.join(TORSION_DCALPHA_PATH, 'rnn_torsion_dcalpha.json')
TORSION_DCALPHA_MW = os.path.join(TORSION_DCALPHA_PATH, 'weights_torsion_dcalpha.h5')

# Returns a rotation matrix that maps v1 to v2
# NOTE: v1, v2 have to be unit vectors
def get_rot_matrix(v1, v2):
  v = np.cross(v1, v2)
  s = np.linalg.norm(v)
  c = np.dot(v1, v2)
  v_x = np.array([[0.,-v[2],v[1]],
                  [v[2],0.,-v[0]],
                  [-v[1],v[0],0.]])
  R = np.eye(3) + v_x + np.dot(v_x, v_x)*(1/(1+c))
    
  return R

# Reads a pickle file and standardizes the 3D coordinates
def read_pkl_std_coords(datafile):
  with open(datafile, 'rb') as fh:
    data = pickle.load(fh)

  casp_ids = list(data.keys())

  for casp_id in casp_ids:
    coords = data[casp_id]['coords']

    offset = coords[0]
    orig_coords = coords - offset # Linear translation
    unit = orig_coords[1] / (np.sqrt(np.sum(orig_coords[1]**2)) + EPSILON) # Unit vector to second amino acid
    R1 = get_rot_matrix(unit, np.array([0,0,1])) # Rotation matrix to [0,0,1]
    R2 = get_rot_matrix(unit, np.array([0,1,0])) # Rotation matrix to [0,1,0]

    R1_coords = np.array([np.matmul(R1,p) for p in orig_coords]) # R1 rotation applied
    R2_coords = np.array([np.matmul(R2,p) for p in R1_coords]) # R2 rotation applied

    data[protein]['coords'] = R2_coords

  return data 

# Loads and returns the training data in folds (fold_start - fold_end)
def load_data(fold_start, fold_end, coords_flag=None, radians_flag=None):
  train_aa_seq = [] # AA (Amino Acid) Sequence
  train_ss_seq = [] # Q8 Secondary Structure Sequence
  train_casp_id = [] # Protein CASP ID
  train_pssms = [] # Position-Specific Scoring Matrix (PSSM)
  train_dcalphas = [] # Ground Truth Distance Matrix
  train_psis = [] # Torsion Angle (psi)
  train_phis = [] # Torsion Angle (phi)
  train_coords = [] # Ground Truth 3D Coordinates

  for i in range(fold_start, fold_end+1):
    datafile = TRAINING_PATH.format(i)
    with open(datafile, 'rb') as fh:
      data = pickle.load(fh)

    casp_ids = list(data.keys())

    print('Number of proteins in training fold {}: {}'.format(i, len(casp_ids)))

    for casp_id in casp_ids:
      train_aa_seq.append(data[casp_id]['aa'])
      train_ss_seq.append(data[casp_id]['ss'])
      train_casp_id.append(casp_id)
      train_pssms.append(np.array(data[casp_id]['pssm']))
      train_dcalphas.append(data[casp_id]['dcalpha'])
      psi = data[casp_id]['psi']
      psi = np.array([0. if v is None else v for v in psi]) # Remove None objects
      if radians_flag: # Convert to radians
        psi = psi * RAD
      train_psis.append(np.array(psi))
      phi = data[casp_id]['phi']
      phi = np.array([0. if v is None else v for v in phi]) # Remove None objects
      if radians_flag: # Convert to radians
        phi = phi * RAD
      train_phis.append(np.array(phi))

      # 3D coordinate standardization
      if coords_flag:
        coords = data[casp_id]['coords']

        offset = coords[0]
        orig_coords = coords - offset # Linear translation
        unit = orig_coords[1] / (np.sqrt(np.sum(orig_coords[1]**2)) + EPSILON) # Unit vector to second amino acid
        R1 = get_rot_matrix(unit, np.array([0,0,1])) # Rotation matrix to [0,0,1]
        R2 = get_rot_matrix(unit, np.array([0,1,0])) # Rotation matrix to [0,1,0]

        R1_coords = np.array([np.matmul(R1,p) for p in orig_coords]) # R1 rotation applied
        R2_coords = np.array([np.matmul(R2,p) for p in R1_coords]) # R2 rotation applied

        train_coords.append(R2_coords)
      else:
        train_coords.append(np.array(data[casp_id]['coords']))

  training_data = (train_aa_seq, train_ss_seq, train_casp_id, train_pssms,
                   train_dcalphas, train_psis, train_phis, train_coords)

  return training_data

# Loads and returns the validation data
def load_val_data(coords_flag=None, radians_flag=None):
  val_aa_seq = [] # AA (Amino Acid) Sequence
  val_ss_seq = [] # Q8 Secondary Structure Sequence
  val_casp_id = [] # Protein CASP ID
  val_pssms = [] # Position-Specific Scoring Matrix (PSSM)
  val_dcalphas = [] # Ground Truth Distance Matrix
  val_psis = [] # Torsion Angle (psi)
  val_phis = [] # Torsion Angle (phi)
  val_coords = [] # Ground Truth 3D Coordinates

  with open(VALIDATION_PATH, 'rb') as fh:
    data = pickle.load(fh)

  casp_ids = list(data.keys())

  print('Number of proteins in validation set: {}'.format(len(casp_ids)))

  for casp_id in casp_ids:
    val_aa_seq.append(data[casp_id]['aa'])
    val_ss_seq.append(data[casp_id]['ss'])
    val_casp_id.append(casp_id)
    val_pssms.append(np.array(data[casp_id]['pssm']))
    val_dcalphas.append(data[casp_id]['dcalpha'])
    psi = data[casp_id]['psi']
    psi = [0. if v is None else v for v in psi] # Remove None objects
    if radians_flag: # Convert to radians
      psi = psi * RAD
    val_psis.append(np.array(psi))
    phi = data[casp_id]['phi']
    phi = [0. if v is None else v for v in phi] # Remove None objects
    if radians_flag: # Convert to radians
      phi = phi * RAD
    val_phis.append(np.array(phi))

    # 3D coordinate standardization
    if coords_flag:
      coords = data[casp_id]['coords']

      offset = coords[0]
      orig_coords = coords - offset # Linear translation
      unit = orig_coords[1] / (np.sqrt(np.sum(orig_coords[1]**2)) + EPSILON) # Unit vector to second amino acid
      R1 = get_rot_matrix(unit, np.array([0,0,1])) # Rotation matrix to [0,0,1]
      R2 = get_rot_matrix(unit, np.array([0,1,0])) # Rotation matrix to [0,1,0]

      R1_coords = np.array([np.matmul(R1,p) for p in orig_coords]) # R1 rotation applied
      R2_coords = np.array([np.matmul(R2,p) for p in R1_coords]) # R2 rotation applied

      val_coords.append(R2_coords)
    else:
      val_coords.append(np.array(data[casp_id]['coords']))

  validation_data = (val_aa_seq, val_ss_seq, val_casp_id, val_pssms,
                     val_dcalphas, val_psis, val_phis, val_coords)

  return validation_data

# Given a matplotlib axis object, returns it with its scales set to be same
def axis_equal_3D(ax):
  extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
  sz = extents[:,1] - extents[:,0]
  centers = np.mean(extents, axis=1)
  maxsize = max(abs(sz))
  r = maxsize / 2

  for ctr, dim in zip(centers, 'xyz'):
    getattr(ax, 'set_{}lim'.format(dim))(ctr-r, ctr+r)

def arjun_distance(A):
    A = tf.squeeze(A)
    dim = tf.shape(A)[0]
    r = tf.reduce_sum(A*A, 1)
    # turn r into column vector
    r = tf.reshape(r, [-1, 1])
    D = r - 2*tf.matmul(A, tf.transpose(A)) + tf.transpose(r)
    D = tf.math.sqrt(D)
    D = tf.matrix_set_diag(D, tf.zeros_like(D[0]), name=None)
    D = tf.expand_dims(D,0)
    D = tf.reshape(D, [1,dim, dim])
    return D
def energetic_loss(y_true, y_pred):
  y_true = tf.squeeze(y_true)
  dim = tf.shape(y_true)[0]
  y_pred = tf.squeeze(y_pred)
  ones = tf.ones_like(y_true)
  mask_a = tf.matrix_band_part(ones, 0, -1) # Upper triangular matrix of 0s and 1s
  mask_b = tf.matrix_band_part(ones, 0, 0)  # Diagonal matrix of 0s and 1s
  mask = tf.cast(mask_a - mask_b, dtype=tf.bool) # Make a bool mask
  mask = tf.reshape(mask, [dim, dim])
  upper_triangular_true = tf.boolean_mask(y_true, mask)
  upper_triangular_pred = tf.boolean_mask(y_pred, mask)
  Divided_t = tf.math.divide(upper_triangular_true, upper_triangular_pred)
  logged = tf.math.log(Divided_t)
  logged_square = tf.math.square(logged)
  loss = tf.reduce_sum(logged_square)
  return loss
# Given a sequence and corresponding CODE, generates the one-hot encodings
def one_hot(seq, code):
  # Input: sequence(str), CODE
  # Output: (n,L) matrix (L: Number of letters in the CODE)

  encoding = np.zeros((len(seq), len(code)), dtype=np.int64)
  for idx, seq_chr in enumerate(seq):
    tmp = np.zeros(len(code), dtype=np.int64)
    tmp[code.find(seq_chr)] = 1
    encoding[idx] = tmp

  return encoding

# Given a sequence and corresponding CODE, generates integer representations for each AA
def int_encoding(seq, code):
  # Input: sequence(str), CODE
  # Output: (n,1) vector

  encoding = np.zeros((len(seq),), dtype=np.int64)
  for idx, seq_chr in enumerate(seq):
    encoding[idx] = code.find(seq_chr)

  return encoding

# Takes an outer product of two vectors and returns the result
def outerprod(args):
  v1 = args[0]
  v2 = args[1]

  return K.batch_dot(v2, v1, axes=[2,2])

# Converts the 3D coordinates to a distance matrix
def coords_to_distmat(coords):
  # Input: (n,3) matrix
  # Output: (n,n) matrix

  n = coords.shape[0]
  dist_mat = np.zeros((n,n))

  for i in range(n):
    coords_i = coords[i]

    delta = coords - coords_i
    delta_squared = delta**2
    l2_squared = np.sum(delta_squared, axis=1)
    l2 = np.sqrt(l2_squared)

    dist_mat[i] = l2

  return dist_mat

# Converts a given matrix into an averaged symmetric matrix
def symmetrize_matrix(mat):
  for i in range(mat.shape[0]):
    for j in range(i):
      mat[i,j]  = mat[j,i] = (mat[j,i] + mat[i,j]) / 2

  for k in range(mat.shape[0]):
    mat[k,k] = 0

  return mat
def upper_tri_masking(A):
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]
# Computes the RMSD using validation/test set
def compute_mae(prediction, validation):
  n_proteins = len(prediction)
  mae_list = []
  for i in range(n_proteins):
    shape = np.shape(prediction[i])[0]
    total = shape**2
    diff = np.abs(prediction[i] - validation[i])
    summed = np.sum(diff)
    mae_i = summed/total
    mae_list.append(mae_i)
  print(mae_list)
  final_mae = statistics.mean(mae_list)
def compute_rmsd(prediction, validation, torsion_flag=False, arjun_flag =True):
  try:
    assert prediction.shape == validation.shape
  except:
    print('[Error] Mismatch in shape:')
    print('Prediction: {}'.format(prediction.shape))
    print('Ground Truth: {}'.format(validation.shape))
  if arjun_flag: 
    n_proteins = len(prediction)
    rmsds = []
    for i in range(n_proteins):
      n_aas = prediction[i].shape[0]
      total = n_aas **2 * 1.0
      total = total - n_aas * 1.0
      ind = np.triu_indices(n_aas)
      pred_triu = upper_tri_masking(prediction[i])
      val_triu = upper_tri_masking(validation[i])
      tmp = pred_triu - val_triu
      tmp = tmp**2
      squared = np.sum(tmp)
      mean_square = squared/total
      rmsd_i = np.sqrt(mean_square)
      rmsds.append(rmsd_i)
    print(rmsds)
    rmsd = statistics.mean(rmsds)
    return rmsd
  elif torsion_flag:
    n_aas = len(prediction) # Number of amino acids
    tmp = prediction - validation
    tmp = tmp**2
    se = sum(tmp)
    mse = se / n_aas
    rmsd = np.sqrt(mse)

  else:
    n_proteins = len(prediction)
    n = 0
    for i in range(n_proteins):
      n_aas = prediction[i].shape[0]
      n += n_aas**2
    tmp = prediction - validation
    tmp = tmp**2
    se = 0
    for i in range(n_proteins):
      se += np.sum(tmp[i])
    mse = se / n
    rmsd = np.sqrt(mse)

  return rmsd

# Preprocess the data from provided folds into desired format
def preprocess_data(training_data, y_flag):
  # Unpack training data
  train_aa_seq = training_data[0]
  train_ss_seq = training_data[1]
  train_casp_id = training_data[2]
  train_pssms = training_data[3]
  train_dcalphas = training_data[4]
  train_psis = training_data[5]
  train_phis = training_data[6]
  train_coords = training_data[7]

  n_proteins = len(train_aa_seq)

  # PSSM Processing
  train_pssms_processed = []
  for i in range(n_proteins):
    train_pssm_i = np.array(train_pssms[i])
    train_pssm_i = np.transpose(train_pssm_i)
    train_pssms_processed.append(train_pssm_i)

  # X_train Processing
  X_train = []
  for i in range(n_proteins):
    aa_encodings = int_encoding(train_aa_seq[i], AA_CODES)
    ss_encodings = int_encoding(train_ss_seq[i], Q8_CODES)
    pssm_processed = train_pssms_processed[i]

    X_train.append([aa_encodings, ss_encodings, pssm_processed])

  # Torsion Angle Processing
  train_angles = []
  for i in range(n_proteins):
    train_psis_i = np.array(train_psis[i])
    train_phis_i = np.array(train_phis[i])
    train_angles_i = np.transpose(np.vstack((train_psis_i, train_phis_i)))

    train_angles.append(train_angles_i)

  if y_flag == 'coords':
    Y_train = np.array(train_coords)
  elif y_flag == 'dcalpha':
    Y_train = np.array(train_dcalphas)
  elif y_flag == 'torsion':
    Y_train = np.array(train_angles)
  elif y_flag == 'combined':
    Y_train = []
    for i in range(n_proteins):
      Y_train.append([train_angles[i], np.array(train_dcalphas[i])])

  return X_train, Y_train
