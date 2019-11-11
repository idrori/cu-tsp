"""
  load_eval_model.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script loads a trained RNN encoder-decoder network, computes
  the dRMSD using the validation set, and makes predictions using the
  test set.

"""

import numpy as np
from numpy import array
import argparse
import pickle
import matplotlib.pyplot as plt
import os
import csv
from datetime import datetime
import sys
import keras
from tqdm import tqdm
import argparse
import tensorflow as tf
from keras.models import load_model

from core import *
sys.path.append('..')
from utils.rgn_utils import *

# CONSTANTS

ID_TO_MODEL = {1: MODEL1_PATH,
               2: MODEL2_PATH,
               3: MODEL3_PATH}

# Saves distance matrix predictions in .pkl format
def preds_to_pkl(test_dcalpha_preds):
  casp_to_preds = {}

  for i in range(len(test_dcalpha_preds)):
    casp_id = test_dcalpha_preds[i][0]
    dcalpha = test_dcalpha_preds[i][1]

    casp_to_preds[casp_id] = {'dcalpha': dcalpha}

  with open(PKL_OUT_PATH, 'wb') as fh:
    pickle.dump(casp_to_preds, fh)

# Loads and returns the validation data
def load_val_data(coords_flag=None, radians_flag=None, msa_flag=False):
  val_aa_seq = [] # AA (Amino Acid) Sequence
  val_ss_seq = [] # Q8 Secondary Structure Sequence
  val_casp_id = [] # Protein CASP ID
  val_pssms = [] # Position-Specific Scoring Matrix (PSSM)
  val_dcalphas = [] # Ground Truth Distance Matrix
  val_psis = [] # Torsion Angle (psi)
  val_phis = [] # Torsion Angle (phi)
  val_coords = [] # Ground Truth 3D Coordinates
  val_msa_onehot = [] # MSA Covariance One Hot
  val_msa_embeds = [] # MSA Covariance Embeddings

  print('Loading the validation protein data...')

  with open(VALIDATION_PATH, 'rb') as fh:
    data = pickle.load(fh)

  casp_ids = list(data.keys())

  for casp_id in tqdm(casp_ids):
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
    if msa_flag:
      # TODO: Fix catch
      try:
        val_msa_onehot.append(np.array(data[casp_id]['msacovonehot'])[:,:,np.newaxis]) # Expand dimension
        if data[casp_id]['msacovonehot'].shape[0] != len(data[casp_id]['aa']):
          print(casp_id)
          print(data[casp_id]['msacovonehot'].shape)
          raise RuntimeError
      except:
        n = len(data[casp_id]['aa'])
        val_msa_onehot.append(np.zeros((n,n,1), dtype=np.float64))
      try:
        val_msa_embeds.append(np.array(data[casp_id]['msacovembedding'])[:,:,np.newaxis]) # Expand dimension
        if data[casp_id]['msacovembedding'].shape[0] != len(data[casp_id]['aa']):
          print(casp_id)
          print(data[casp_id]['msacovembedding'].shape)
          raise RuntimeError
      except:
        n = len(data[casp_id]['aa'])
        val_msa_embeds.append(np.zeros((n,n,1), dtype=np.float64))

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

  if msa_flag:
    validation_data = (val_aa_seq, val_ss_seq, val_casp_id, val_pssms,
                       val_dcalphas, val_psis, val_phis, val_coords,
                       val_msa_onehot, val_msa_embeds)
  else:
    validation_data = (val_aa_seq, val_ss_seq, val_casp_id, val_pssms,
                       val_dcalphas, val_psis, val_phis, val_coords)

  # Preprocess the loaded data
  print('Preprocessing the validation protein data...')
  X_val, Y_val = preprocess_data(validation_data, 'dcalpha', msa_flag=msa_flag)

  return X_val, Y_val

# Loads the returns the test data
def load_test_data(msa_flag=False):
  test_aa_seq = [] # AA (Amino Acid) Sequence
  test_ss_seq = [] # Q8 Secondary Structure Sequence
  test_casp_ids = [] # Protein CASP ID
  test_pssms = [] # Position-Specific Scoring Matrix (PSSM)
  test_msa_onehot = [] # MSA Covariance One Hot
  test_msa_embeds = [] # MSA Covariance Embeddings

  print('Loading the test protein data...')

  with open(TEST_PATH, 'rb') as fh:
    data = pickle.load(fh)

  casp_ids = list(data.keys())

  for casp_id in tqdm(casp_ids):
    test_aa_seq.append(data[casp_id]['aa'])
    test_ss_seq.append(data[casp_id]['ss'])
    test_casp_ids.append(casp_id)
    test_pssms.append(np.array(data[casp_id]['pssm']))
    if msa_flag:
      # TODO: Fix catch
      try:
        test_msa_onehot.append(np.array(data[casp_id]['msacovonehot'])[:,:,np.newaxis]) # Expand dimension
        if data[casp_id]['msacovonehot'].shape[0] != len(data[casp_id]['aa']):
          print(casp_id)
          print(data[casp_id]['msacovonehot'].shape)
          raise RuntimeError
      except:
        n = len(data[casp_id]['aa'])
        test_msa_onehot.append(np.zeros((n,n,1), dtype=np.float64))
      try:
        test_msa_embeds.append(np.array(data[casp_id]['msacovembedding'])[:,:,np.newaxis]) # Expand dimension
        if data[casp_id]['msacovembedding'].shape[0] != len(data[casp_id]['aa']):
          print(casp_id)
          print(data[casp_id]['msacovembedding'].shape)
          raise RuntimeError
      except:
        n = len(data[casp_id]['aa'])
        test_msa_embeds.append(np.zeros((n,n,1), dtype=np.float64))
  
  # Preprocess the loaded data
  print('Preprocessing the test protein data...')
  n_proteins = len(test_aa_seq)

  # PSSM Processing
  test_pssms_processed = []
  for i in range(n_proteins):
    test_pssm_i = np.array(test_pssms[i])
    test_pssm_i = np.transpose(test_pssm_i)
    test_pssms_processed.append(test_pssm_i)

  # X_test Processing
  X_test = []
  for i in range(n_proteins):
    aa_encodings = int_encoding(test_aa_seq[i], AA_CODES)
    ss_encodings = int_encoding(test_ss_seq[i], Q8_CODES)
    pssm_processed = test_pssms_processed[i]
    if msa_flag:
      msa_onehot = test_msa_onehot[i]
      msa_embed = test_msa_embeds[i]
      X_test.append([aa_encodings, ss_encodings, pssm_processed,
                     msa_onehot, msa_embed])
    else:
      X_test.append([aa_encodings, ss_encodings, pssm_processed])

  return X_test, test_casp_ids

def evaluate_model(**kwargs):
  model_id = int(kwargs['model_id'])
  msa_flag = kwargs['msa']
  no_q8_flag = kwargs['no_q8']
  log_flag = kwargs['log']
  frobenius_flag = kwargs['frobenius']
  pred_flag = kwargs['predict']

  # Check if given options are valid
  if model_id < 1 or model_id > 3:
    raise RuntimeError('[Error] Invalid Model ID provided.')

  if model_id != 3:
    if log_flag or frobenius_flag:
      raise RuntimeError('[Error] Invalid option for Model {}'.format(model_id))
    elif msa_flag and no_q8_flag:
      raise RuntimeError('[Error] Conflicting options: "msa" and "no_q8."')

  if model_id == 3:
    if msa_flag or no_q8_flag:
      raise RuntimeError('[Error] Invalid option for Model {}'.format(model_id))
    elif log_flag and frobenius_flag:
      raise RuntimeError('[Error] Conflicting options: "log" and "frobenius."')

  # Check path to save training results
  if not os.path.exists(RESULTS_PATH):
    try:
      mkdir(RESULTS_PATH)
    except:
      raise RuntimeError('[Error] Results path setting is incorrect.')

  start_time = datetime.now()

  # Load and preprocess validation data
  X_val, Y_val = load_val_data(msa_flag)

  # Load and preprocess test data
  X_test, test_casp_ids = load_test_data(msa_flag)

  # Load trained model
  model_path = ID_TO_MODEL[model_id]
  custom_objects = {'tf': tf, 'reduce_l2_norm': reduce_l2_norm,
                    'dihedral_to_point': dihedral_to_point,
                    'point_to_coordinate': point_to_coordinate,
                    'pairwise_distance': pairwise_distance,
                    'theta': BOND_ANGLES,
                    'NUM_DIHEDRALS': NUM_DIHEDRALS,
                    'NUM_DIMENSIONS': NUM_DIMENSIONS,
                    'collections': collections}
  
  model = load_model(model_path, custom_objects)

  # Evaluate dRMSD using validation data
  print('Evaluating model using validation data...')
  
  val_dcalpha_preds = []
  for i in tqdm(range(len(X_val))):
    # Models 1 and 2
    if model_id != 3:
      if msa_flag:
        x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], 
             X_val[i][2][np.newaxis,:], np.array(X_val[i][3])[np.newaxis,:], 
             np.array(X_val[i][4])[np.newaxis,:]]
      elif no_q8_flag:
        x = [np.array(X_val[i][0])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
      else:
        x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], 
             X_val[i][2][np.newaxis,:]]

      y = Y_val[i][np.newaxis,:]

      if model_id == 1:
        dcalpha = model.predict(x)
        dcalpha = np.squeeze(dcalpha)
        val_dcalpha_preds.append(dcalpha)

      elif model_id == 2:
        _, dcalpha = model.predict(x)
        dcalpha = np.squeeze(dcalpha)
        val_dcalpha_preds.append(dcalpha)

    elif model_id == 3:
      x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], 
           X_val[i][2][np.newaxis,:]]

      if log_flag:
        #m = combined[i][np.newaxis,:] + EPSILON
        #n = np.log(m)
        #y = [n,n]
        y = [combined[i][np.newaxis,:], combined[i][np.newaxis,:]]
      elif frobenius_flag:
        #m = combined[i][np.newaxis,:]
        #y = [m,m]
        y = [combined[i][np.newaxis,:], combined[i][np.newaxis,:]]
      else:
        y = combined[i][np.newaxis,:]

      dcalpha = model.predict(x)
      dcalpha = np.squeeze(dcalpha)

      if log_flag:
        base = 2 * np.ones(shape=dcalpha.shape)
        dcalpha = np.power(base, dcalpha)
      
      val_dcalpha_preds.append(dcalpha)

  val_dcalpha_preds = np.array(val_dcalpha_preds)
  drmsd = compute_rmsd(val_dcalpha_preds, val_dcalphas, arjun_flag=True)
  print('Validation dRMSD: {}'.format(drmsd))

  # Generate predictions for test data
  print('Generating predictions for test data...')

  test_dcalpha_preds = []
  for i in tqdm(range(len(X_test))):
    # Models 1 and 2
    if model_id != 3:
      if msa_flag:
        x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], 
             X_test[i][2][np.newaxis,:], np.array(X_test[i][3])[np.newaxis,:], 
             np.array(X_test[i][4])[np.newaxis,:]]
      elif no_q8_flag:
        x = [np.array(X_test[i][0])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
      else:
        x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], 
             X_test[i][2][np.newaxis,:]]

      if model_id == 1:
        dcalpha = model.predict(x)
        dcalpha = np.squeeze(dcalpha)

      elif model_id == 2:
        _, dcalpha = model.predict(x)
        dcalpha = np.squeeze(dcalpha)

    elif model_id == 3:
      x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], 
           X_test[i][2][np.newaxis,:]]

      dcalpha = model.predict(x)
      dcalpha = np.squeeze(dcalpha)

      if log_flag:
        base = 2 * np.ones(shape=dcalpha.shape)
        dcalpha = np.power(base, dcalpha)

    # Save distance matrix as .png
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=dcalpha.min(), vmax=dcalpha.max())
    image = cmap(norm(dcalpha))
    plt.imsave(os.path.join(RESULTS_PATH), '{}_dcalpha.png'.format(test_casp_ids[i]))
      
    test_dcalpha_preds.append([test_casp_ids[i], dcalpha])

  preds_to_pkl(test_dcalpha_preds)
  print('Test predictions saved.')
  
  end_time = datetime.now()
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Training Time: {}h {}m {}s'.format(hours, minutes, seconds))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('model_id', help='Model ID')
  parser.add_argument('--msa', help='Indicate model requires MSA as input', action='store_true',
                      default=False)
  parser.add_argument('--no_q8', help='Indicate model does not take Q8 as input', 
                      action='store_true', default=False)
  parser.add_argument('--log', help='Indicate model makes log-dcalpha predictions', action='store_true',
                      default=False)
  parser.add_argument('--frobenius', help='Indicate model used Frobenius-norm-based loss', action='store_true',
                      default=False)
  parser.add_argument('--predict', help='Option to predict using test data', action='store_true',
                      default=False)
  args = parser.parse_args()

  evaluate_model(**vars(args))
  
