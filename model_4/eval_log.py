"""
  load_eval_sublime.py
  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University
  This script loads the trained reconstruction networks, computes
  the dRMSD using the validation set, and makes predictions using the
  test set.
"""

import numpy as np
from numpy import array
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
from keras.models import Sequential, save_model, model_from_json
from keras.backend.tensorflow_backend import set_session
config = tf.ConfigProto()
config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
config.log_device_placement = True  # to log device placement (on which device the operation ran)
sess = tf.Session(config=config)
set_session(sess)
from core import *
from core_harvard import *

# CONSTANTS
##EDIT WITH PATHS TO MODEL FOR LOG FILE
#MJ_FILE = '/home/arjunsrivatsa/results_old/rnnlog.json'
#MW_FILE = '/home/arjunsrivatsa/results_old/weightslog.h5'
MJ_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/rnnlog.json'
MW_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/weightslog.h5'
#MJ_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/100.json'
#MW_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/100w.h5'
#MJ_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/altrnn.json'
#MW_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/altweights.h5'
#MJ_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/normalized.json'
#MW_FILE = '/home/darshan/cu-tsp-dev/model_4_branch/results/weightsnormalized.h5'

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

# Loads a dataset given the appropriate flag
def load_val_test(data_flag):
  if not data_flag:
    print('[Error] No dataset flag provided.')
    return None

  if data_flag == 'validation':
    print('Loading and preprocessing validation protein data...')

    validation_data = load_val_data(coords_flag=True)
    X_val, val_coords = preprocess_data(validation_data, 'coords')
    _, val_dcalphas = preprocess_data(validation_data, 'dcalpha')
    _, val_angles = preprocess_data(validation_data, 'torsion')

    val_psis, val_phis = [], []
    for i in range(len(val_angles)):
      val_psis += list(val_angles[i][:,0])
      val_phis += list(val_angles[i][:,1])
    val_psis = np.array(val_psis)
    val_phis = np.array(val_phis)

    val_angles = np.transpose(np.vstack((val_psis, val_phis)))

    return (X_val, val_coords, val_dcalphas, val_angles)

  elif data_flag == 'test':
    print('Loading and preprocessing test protein data...')
    test_aa_seq = [] # AA (Amino Acid) Sequence
    test_ss_seq = [] # Q8 Secondary Structure Sequence
    test_casp_id = [] # Protein CASP ID
    test_pssms = [] # Position-Specific Scoring Matrix (PSSM)
    test_dcalphas = [] # Ground Truth Distance Matrix
    test_psis = [] # Torsion Angle (psi)
    test_phis = [] # Torsion Angle (phi)
    test_coords = [] # Ground Truth 3D Coordinates

    with open(TEST_PATH, 'rb') as fh:
      test_data = pickle.load(fh)

    casp_ids = list(test_data.keys())

    for casp_id in casp_ids:
      test_aa_seq.append(test_data[casp_id]['aa'])
      test_ss_seq.append(test_data[casp_id]['ss'])
      test_casp_id.append(casp_id)
      test_pssms.append(np.array(test_data[casp_id]['pssm']))
      test_dcalphas.append(test_data[casp_id]['dcalpha'])
      psi = test_data[casp_id]['psi']
      for i in range(len(psi)):
        if psi[i] == None:
          psi[i] = 0.
      test_psis.append(np.array(psi))
      phi = test_data[casp_id]['phi']
      for i in range(len(phi)):
        if phi[i] == None:
          phi[i] = 0.
      test_phis.append(np.array(phi))

      # 3D coordinate standardization
      coords = test_data[casp_id]['coords']

      offset = coords[0]
      orig_coords = coords - offset # Linear translation
      unit = orig_coords[1] / (np.sqrt(np.sum(orig_coords[1]**2)) + EPSILON) # Unit vector to second amino acid
      R1 = get_rot_matrix(unit, np.array([0,0,1])) # Rotation matrix to [0,0,1]
      R2 = get_rot_matrix(unit, np.array([0,1,0])) # Rotation matrix to [0,1,0]

      R1_coords = np.array([np.matmul(R1,p) for p in orig_coords]) # R1 rotation applied
      R2_coords = np.array([np.matmul(R2,p) for p in R1_coords]) # R2 rotation applied

      test_coords.append(R2_coords)

    n_proteins = len(test_aa_seq)
    # PSSM Processing
    test_pssms_processed = []
    for i in range(n_proteins):
      test_pssm_i = np.array(test_pssms[i])
      test_pssm_i = np.transpose(test_pssm_i)
      test_pssms_processed.append(test_pssm_i)

    # X_train Processing
    X_test = []
    for i in range(n_proteins):
      aa_encodings = int_encoding(test_aa_seq[i], AA_CODES)
      ss_encodings = int_encoding(test_ss_seq[i], Q8_CODES)
      pssm_processed = test_pssms_processed[i]

      X_test.append([aa_encodings, ss_encodings, pssm_processed])

    return (X_test, casp_ids)

# Saves predicted 3D coordinates, dcalpha, torsion angles in .csv format
def preds_to_csv(test_coords_preds, test_dcalpha_preds, test_angle_preds):
  csv_fh = open(CSV_OUT_PATH, 'w', newline='')
  csv_writer = csv.writer(csv_fh)

  csv_writer.writerow(['CASP ID', 'Predicted']) # Header

  for i in range(len(test_coords_preds)):
    casp_id = test_coords_preds[i][0]
    coords = test_coords_preds[i][1]
    dcalpha = test_dcalpha_preds[i][1]
    angles = test_angle_preds[i][1]

    # Write 3D coordinate predictions
    idx_to_axis = ['x','y','z']
    for j in range(coords.shape[0]):
      for k in range(coords.shape[1]):
        csv_entry = CSV_COORDS_FORMAT.format(casp_id, j+1, idx_to_axis[k])
        csv_writer.writerow([csv_entry, coords[j][k]])

    # Write distance predictions
    for j in range(dcalpha.shape[0]):
      for k in range(dcalpha.shape[1]):
        csv_entry = CSV_DCALPHA_FORMAT.format(casp_id, j+1, k+1)
        csv_writer.writerow([csv_entry, dcalpha[j][k]])

    # Write psi predictions
    for j in range(angles.shape[0]):
      csv_entry = CSV_PSI_FORMAT.format(casp_id, j+1)
      csv_writer.writerow([csv_entry, angles[j][0]])

    # Write phi predictions
    for j in range(angles.shape[0]):
      csv_entry = CSV_PHI_FORMAT.format(casp_id, j+1)
      csv_writer.writerow([csv_entry, angles[j][1]])

  csv_fh.close()

# Saves the dcalpha, torsion angle predictions in .pkl format
def preds_to_pkl(test_dcalpha_preds):
  # Format: {casp_id: {'dcalpha': np.array, 'psi': np.array, 'phi': np.array}}
  casp_to_preds = {}

  for i in range(len(test_dcalpha_preds)):
    casp_id = test_dcalpha_preds[i][0]
    dcalpha = test_dcalpha_preds[i][1]

    casp_to_preds[casp_id] = {'dcalpha': dcalpha}

  with open(PKL_OUT_PATH, 'wb') as fh:
    pickle.dump(casp_to_preds, fh)

# Evaluates the MDS-based distance matrix RNN model
def evaluate_mds_dcalpha():
  # Load test and validation data
  X_val, val_coords, val_dcalphas, val_angles = load_val_test('validation')
  X_test, test_casp_ids = load_val_test('test')

  data = (X_val, val_coords, val_dcalphas, val_angles, X_test, test_casp_ids)

  X_val = data[0]
  val_coords = data[1]
  val_dcalphas = data[2]
  val_angles = data[3]
  X_test = data[4]
  test_casp_ids = data[5]

  start_time = datetime.now()

  model = load_model()

  print('Evaluating MDS dcalpha model using validation data...')

  # Compute dRMSD using validation set
  val_dcalpha_preds = []
  for i in tqdm(range(len(X_val))):
    x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
    dcalpha, _ = model.predict(x)
    #dcalpha = model.predict(x)
    dcalpha = np.squeeze(dcalpha)
    dcalpha = np.exp(dcalpha)
    val_dcalpha_preds.append(dcalpha)

  val_dcalpha_preds = np.array(val_dcalpha_preds)
  drmsd = compute_rmsd(val_dcalpha_preds, val_dcalphas)
  print('Distance Matrix Validation dRMSD: {}'.format(drmsd))

  # Predict using test data and save as .npy and .png
  test_dcalpha_preds = []
  for i in tqdm(range(len(X_test))):
    x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
    dcalpha, _ = model.predict(x)
    #dcalpha = model.predict(x)
    dcalpha = np.squeeze(dcalpha)
    dcalpha = np.exp(dcalpha)
    test_dcalpha_preds.append([test_casp_ids[i], dcalpha])

    # Save distance matrix as .npy
    np.save(os.path.join(MDS_DCALPHA_PATH, '{}_mds_dcalpha.npy'.format(test_casp_ids[i])), dcalpha)

    # Save distance matrix as .png
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=dcalpha.min(), vmax=dcalpha.max())
    image = cmap(norm(dcalpha))
    plt.imsave(os.path.join(MDS_DCALPHA_PATH, '{}_mds_dcalpha.png'.format(test_casp_ids[i])), image)

  print('Distance matrix test predictions saved.')

  preds_to_pkl(test_dcalpha_preds)
  print('Results saved in .pkl format.')

  end_time = datetime.now()
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Evaluation Time: {}h {}m {}s'.format(hours, minutes, seconds))

if __name__ == '__main__':
  #parser = argparse.ArgumentParser()

  #parser.add_argument('--coords', help='Evaluate 3D coordinates model', action='store_true')
  #parser.add_argument('--coords_dcalpha', help='Evaluate coords_dcalpha model', action='store_true')
  #parser.add_argument('--torsion', help='Evaluate torsion angle model', action='store_true')
  #parser.add_argument('--torsion_dcalpha', help='Evaluate torsion_dcalpha model', action='store_true')
  #args = parser.parse_args()
  
  #evaluate_models(**vars(args))

  evaluate_mds_dcalpha()
