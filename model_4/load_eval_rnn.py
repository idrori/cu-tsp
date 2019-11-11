"""
  load_eval_rnn.py

  Daniel Jeong, Arjun Srivatsa
  Department of Computer Science
  Columbia University

  This script loads the trained RNN encoder-decoder networks, computes
  the RMSD using the validation set, and makes predictions using the
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

from core import *
sys.path.append('..')
from utils.rgn_utils import *

# Loads a trained model given the flag
def load_model(model_flag):
  if not model_flag:
    print('[Error] No model flag provided.')
    return None

  if model_flag == 'coords':
    print('Loading coords model...')
    mj_file = COORDS_MJ
    mw_file = COORDS_MW

  elif model_flag == 'coords_dcalpha':
    print('Loading coords_dcalpha model...')
    mj_file = COORDS_DCALPHA_MJ
    mw_file = COORDS_DCALPHA_MW

  elif model_flag == 'torsion':
    print('Loading torsion model...')
    mj_file = TORSION_MJ
    mw_file = TORSION_MW

  elif model_flag == 'torsion_dcalpha':
    print('Loading torsion_dcalpha model...')
    mj_file = TORSION_DCALPHA_MJ
    mw_file = TORSION_DCALPHA_MW

  with open(mj_file, 'r') as json_file:
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
  loaded_model.load_weights(mw_file)

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
def preds_to_pkl(test_coords_preds, test_dcalpha_preds, test_angle_preds):
  # Format: {casp_id: {'dcalpha': np.array, 'psi': np.array, 'phi': np.array}}
  casp_to_preds = {}

  for i in range(len(test_dcalpha_preds)):
    casp_id = test_dcalpha_preds[i][0]
    if not test_coords_preds == None:
      x = test_coords_preds[i][1][:,0]
      y = test_coords_preds[i][1][:,1]
      z = test_coords_preds[i][1][:,2]
    dcalpha = test_dcalpha_preds[i][1]
    psi = test_angle_preds[i][1][:,0]
    phi = test_angle_preds[i][1][:,1]

    casp_to_preds[casp_id] = {'dcalpha': dcalpha,
                              'psi': psi,
                              'phi': phi}
    if not test_coords_preds == None:
      casp_to_preds[casp_id]['x'] = x
      casp_to_preds[casp_id]['y'] = y
      casp_to_preds[casp_id]['z'] = z

  with open(PKL_OUT_PATH, 'wb') as fh:
    pickle.dump(casp_to_preds, fh)

# Evaluates the 3D coordinate RNN model
def evaluate_coords(data):
  X_val = data[0]
  val_coords = data[1]
  val_dcalphas = data[2]
  val_angles = data[3]
  X_test = data[4]
  test_casp_ids = data[5]

  model_coords = load_model('coords')

  print('Evaluating coords model using validation data...')

  # Compute RMSD using validation set
  val_coords_preds = []
  val_dcalpha_preds = []
  for i in tqdm(range(len(X_val))):
    x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
    coords = model_coords.predict(x)

    coords = np.squeeze(coords)
    val_coords_preds.append(coords)
    val_dcalpha_preds.append(coords_to_distmat(coords))

  val_coords_preds = np.array(val_coords_preds)
  crmsd = compute_rmsd(val_coords_preds, val_coords)
  val_dcalpha_preds = np.array(val_dcalpha_preds)
  drmsd = compute_rmsd(val_dcalpha_preds, val_dcalphas)
  print('3D Coordinates Validation cRMSD: {}'.format(crmsd))
  print('3D Coordinates Validation dRMSD: {}'.format(drmsd))

  # Predict using test data and save as .npy and .png
  test_coords_preds = []
  for i in tqdm(range(len(X_test))):
    x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
    coords = model_coords.predict(x)

    coords = np.squeeze(coords)
    test_coords_preds.append([test_casp_ids[i], coords])

    coords_dcalpha = coords_to_distmat(coords)

    # Save 3D coordinates as .npy
    #np.save(os.path.join(COORDS_PATH, '{}_coords.npy'.format(test_casp_ids[i], coords)))

    # Save distance matrix as .npy
    #np.save(os.path.join(COORDS_PATH, '{}_coords_to_dcalpha.npy'.format(test_casp_ids[i])), coords_dcalpha)

    # Save distance matrix as .png
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=coords_dcalpha.min(), vmax=coords_dcalpha.max())
    image = cmap(norm(coords_dcalpha))
    plt.imsave(os.path.join(COORDS_PATH, '{}_coords_to_dcalpha.png'.format(test_casp_ids[i])), image)

  print('3D coordinate test predictions saved.')

  return np.array(test_coords_preds)

# Evaluates the distance matrix RNN model
def evaluate_coords_dcalpha(data):
  X_val = data[0]
  val_coords = data[1]
  val_dcalphas = data[2]
  val_angles = data[3]
  X_test = data[4]
  test_casp_ids = data[5]

  model_coords_dcalpha = load_model('coords_dcalpha')

  print('Evaluating coords_dcalpha model using validation data...')

  # Compute dRMSD using validation set
  val_dcalpha_preds = []
  for i in tqdm(range(len(X_val))):
    x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
    dcalpha = model_coords_dcalpha.predict(x)

    dcalpha = np.squeeze(dcalpha)
    val_dcalpha_preds.append(dcalpha)

  val_dcalpha_preds = np.array(val_dcalpha_preds)
  drmsd = compute_rmsd(val_dcalpha_preds, val_dcalphas)
  print('Distance Matrix Validation dRMSD: {}'.format(drmsd))

  # Predict using test data and save as .npy and .png
  test_dcalpha_preds = []
  for i in tqdm(range(len(X_test))):
    x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
    dcalpha = model_coords_dcalpha.predict(x)

    dcalpha = np.squeeze(dcalpha)
    test_dcalpha_preds.append([test_casp_ids[i], dcalpha])

    # Save distance matrix as .npy
    #np.save(os.path.join(DCALPHA_PATH, '{}_coords_dcalpha.npy'.format(test_casp_ids[i])), dcalpha)

    # Save distance matrix as .png
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=dcalpha.min(), vmax=dcalpha.max())
    image = cmap(norm(dcalpha))
    plt.imsave(os.path.join(COORDS_DCALPHA_PATH, '{}_coords_dcalpha.png'.format(test_casp_ids[i])), image)

  print('Distance matrix test predictions saved.')

  return (np.array(test_dcalpha_preds), drmsd)

# Evaluates the torsion angle RNN model
def evaluate_torsion(data):
  X_val = data[0]
  val_coords = data[1]
  val_dcalphas = data[2]
  val_angles = data[3]
  X_test = data[4]
  test_casp_ids = data[5]

  model_torsion = load_model('torsion')

  print('Evaluating torsion model using validation data...')

  # Compute torsion angle RMSD using validation set
  val_angle_preds = []
  for i in tqdm(range(len(X_val))):
    x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
    angles = model_torsion.predict(x)

    angles = np.squeeze(angles)
    for j in range(angles.shape[0]):
      val_angle_preds.append(angles[j])

  val_angle_preds = np.array(val_angle_preds)
  rmsd = compute_rmsd(val_angle_preds, val_angles, torsion_flag=True)
  print('Torsion Angle Validation RMSD: psi = {}, phi = {}'.format(rmsd[0], rmsd[1]))

  # Predict using test data and save as .npy
  test_angle_preds = []
  for i in tqdm(range(len(X_test))):
    x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
    angles = model_torsion.predict(x)

    angles = np.squeeze(angles)
    test_angle_preds.append([test_casp_ids[i], angles])

    # Save angles as .npy
    #np.save(os.path.join(TORSION_PATH, '{}_angles.npy'.format(test_casp_ids[i])), angles)

  print('Torsion angle test predictions saved.')

  return (np.array(test_angle_preds), rmsd)

# Evaluates the torsion angle and distance matrix RNN model
def evaluate_torsion_dcalpha(data):
  X_val = data[0]
  val_coords = data[1]
  val_dcalphas = data[2]
  val_angles = data[3]
  X_test = data[4]
  test_casp_ids = data[5]

  model_torsion_dcalpha = load_model('torsion_dcalpha')

  print('Evaluating torsion_dcalpha model using validation data...')

  # Compute torsion angle and distance matrix RMSD using validation set
  val_angle_preds = []
  val_dcalpha_preds = []
  for i in tqdm(range(len(X_val))):
    x = [np.array(X_val[i][0])[np.newaxis,:], np.array(X_val[i][1])[np.newaxis,:], X_val[i][2][np.newaxis,:]]
    # Note: Angle order is switched (phi, psi) instead of (psi, phi)
    # Note: Angles are in radians
    angles, dcalpha = model_torsion_dcalpha.predict(x)

    angles = np.squeeze(angles)
    dcalpha = np.squeeze(dcalpha)
    for j in range(angles.shape[0]):
      val_angle_preds.append(np.flip(angles[j] * DEG))

    val_dcalpha_preds.append(dcalpha)

  val_angle_preds = np.array(val_angle_preds)
  val_dcalpha_preds = np.array(val_dcalpha_preds)
  torsion_rmsd = compute_rmsd(val_angle_preds, val_angles, torsion_flag=True)
  drmsd = compute_rmsd(val_dcalpha_preds, val_dcalphas)
  print('Torsion Angle Validation RMSD: psi = {}, phi = {}'.format(torsion_rmsd[0], torsion_rmsd[1]))
  print('Distance Matrix Validation dRMSD: {}'.format(drmsd))

  # Predict using test data and save as .npy and .png
  test_angle_preds = []
  test_dcalpha_preds = []
  for i in tqdm(range(len(X_test))):
    x = [np.array(X_test[i][0])[np.newaxis,:], np.array(X_test[i][1])[np.newaxis,:], X_test[i][2][np.newaxis,:]]
    # Note: Angle order is switched (phi, psi) instead of (psi, phi)
    # Note: Angles are in radians
    angles, dcalpha = model_torsion_dcalpha.predict(x)

    angles = np.squeeze(angles)
    test_angle_preds.append([test_casp_ids[i], np.flip(angles * DEG, axis=1)])
    test_dcalpha_preds.append([test_casp_ids[i], dcalpha])

    # Save angles as .npy
    #np.save(os.path.join(TORSION_DCALPHA_PATH, '{}_torsion_angles.npy'.format(test_casp_ids[i])), angles)

    # Save distance matrix as .npy
    #np.save(os.path.join(TORSION_DCALPHA_PATH, '{}_torsion_dcalpha.npy'.format(test_casp_ids[i])), dcalpha)

    # Save distance matrix as .png
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=dcalpha.min(), vmax=dcalpha.max())
    image = cmap(norm(dcalpha))
    plt.imsave(os.path.join(TORSION_DCALPHA_PATH, '{}_torsion_dcalpha.png'.format(test_casp_ids[i])), image)

  print('Torsion angle test predictions saved.')

  return (np.array(test_angle_preds), np.array(test_dcalpha_preds), torsion_rmsd, drmsd)

# Evaluates the distance matrix and torsion angle RNN models
def evaluate_models(**kwargs):
  coords_flag = kwargs['coords']
  coords_dcalpha_flag = kwargs['coords_dcalpha']
  torsion_flag = kwargs['torsion']
  torsion_dcalpha_flag = kwargs['torsion_dcalpha']

  start_time = datetime.now()

  # Load test and validation data
  X_val, val_coords, val_dcalphas, val_angles = load_val_test('validation')
  X_test, test_casp_ids = load_val_test('test')

  data = (X_val, val_coords, val_dcalphas, val_angles, X_test, test_casp_ids)

  # Evaluate 3D coordinate model
  if coords_flag:
    test_coords_preds = evaluate_coords(data)
  else:
    test_coords_preds = None

  # Evaluate 3D coordinate to distance matrix model
  if coords_dcalpha_flag:
    test_coords_dcalpha_preds, coords_dcalpha_drmsd = evaluate_coords_dcalpha(data)

  # Evaluate torsion angle model
  if torsion_flag:
    test_angle_preds, torsion_rmsd = evaluate_torsion(data)

  # Evaluate torsion angle to distance matrix model
  if torsion_dcalpha_flag:
    torsion_dcalpha_results = evaluate_torsion_dcalpha(data)
    test_torsion_dcalpha_angle_preds = torsion_dcalpha_results[0]
    test_torsion_dcalpha_dcalpha_preds = torsion_dcalpha_results[1]
    torsion_dcalpha_rmsd = torsion_dcalpha_results[2]
    torsion_dcalpha_drmsd = torsion_dcalpha_results[3]

  # Save results as .pkl
  best_coords_preds = test_coords_preds

  # Case 1
  if (coords_dcalpha_flag and torsion_flag) and torsion_dcalpha_flag:
    # Compare dcalpha predictions
    if coords_dcalpha_drmsd <= torsion_dcalpha_drmsd:
      best_dcalpha_preds = test_coords_dcalpha_preds
    else:
      best_dcalpha_preds = test_torsion_dcalpha_dcalpha_preds

    # Compare angle predictions
    if torsion_rmsd <= torsion_dcalpha_rmsd:
      best_angle_preds = test_angle_preds
    else:
      best_torsion_preds = test_torsion_dcalpha_angle_preds

    preds_to_pkl(best_coords_preds, best_dcalpha_preds, best_angle_preds)
    print('Results saved in .pkl format.')

  # Case 2
  elif (coords_dcalpha_flag and torsion_flag) and not torsion_dcalpha_flag:
    best_dcalpha_preds = test_coords_dcalpha_preds
    best_angle_preds = test_angle_preds

    preds_to_pkl(best_coords_preds, best_dcalpha_preds, best_angle_preds)
    print('Results saved in .pkl format.')

  # Case 3
  elif not (coords_dcalpha_flag and torsion_flag) and torsion_dcalpha_flag:
    best_dcalpha_preds = test_torsion_dcalpha_dcalpha_preds
    best_angle_preds = test_torsion_dcalpha_angle_preds

    preds_to_pkl(best_coords_preds, best_dcalpha_preds, best_angle_preds)
    print('Results saved in .pkl format.')

  # Case 4
  elif (not coords_dcalpha_flag) and torsion_flag and torsion_dcalpha_flag:
    best_dcalpha_preds = test_torsion_dcalpha_dcalpha_preds

    # Compare angle predictions
    if torsion_rmsd <= torsion_dcalpha_rmsd:
      best_angle_preds = test_angle_preds
    else:
      best_torsion_preds = test_torsion_dcalpha_angle_preds

    preds_to_pkl(best_coords_preds, best_dcalpha_preds, best_angle_preds)
    print('Results saved in .pkl format.')

  # Case 5
  elif coords_dcalpha_flag and (not torsion_flag) and torsion_dcalpha_flag:
    # Compare dcalpha predictions
    if coords_dcalpha_drmsd <= torsion_dcalpha_drmsd:
      best_dcalpha_preds = test_coords_dcalpha_preds
    else:
      best_dcalpha_preds = test_torsion_dcalpha_dcalpha_preds

    best_angle_preds = test_torsion_dcalpha_angle_preds

    preds_to_pkl(best_coords_preds, best_dcalpha_preds, best_angle_preds)
    print('Results saved in .pkl format.')

  end_time = datetime.now()
  elapsed_time = end_time - start_time
  hours, rem = divmod(elapsed_time.seconds, 3600)
  minutes, seconds = divmod(rem, 60)
  print('Total Training Time: {}h {}m {}s'.format(hours, minutes, seconds))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--coords', help='Evaluate 3D coordinates model', action='store_true')
  parser.add_argument('--coords_dcalpha', help='Evaluate coords_dcalpha model', action='store_true')
  parser.add_argument('--torsion', help='Evaluate torsion angle model', action='store_true')
  parser.add_argument('--torsion_dcalpha', help='Evaluate torsion_dcalpha model', action='store_true')
  args = parser.parse_args()

  evaluate_models(**vars(args))
