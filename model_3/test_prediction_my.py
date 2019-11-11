import os, sys, time 
import pickle
import numpy as np
import tensorflow as tf
from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.layers import *
from model_utils import *

print(tf.__version__)

########## load data and model
data_dir = "./data/"
model_dir = sys.argv[1]  # or your model dir


data_fields = [[] for _ in range(6)]
with open(os.path.join(data_dir, 'test.pkl'), 'rb') as f:
    data = pickle.load(f)
    for j in range(len(data_fields)):
        data_fields[j] = data[j]

indices, pdbs, length_aas, pdb_aas, q8s, msas = data_fields
msas = [np.stack(x).transpose().astype(np.float32) for x in msas]

maxlen_seq = 384

# Input data 
test_input_seqs = pdb_aas
test_input_grams = seq2ngrams(test_input_seqs, n=1)
test_q8s = q8s
test_q8s_grams = seq2ngrams(test_q8s, n=1)

# Tokenizer encoder for the input sequence
with open(os.path.join(model_dir, 'tokenizer_encoder.pickle'), 'rb') as handle:
    tokenizer_encoder = pickle.load(handle)
with open(os.path.join(model_dir, 'tokenizer_encoder_q8s.pickle'), 'rb') as handle:
    tokenizer_encoder_q8s = pickle.load(handle)

# Tokenize the input sequences for use in testing
test_input_data = tokenizer_encoder.texts_to_sequences(test_input_grams)
test_input_data = sequence.pad_sequences(test_input_data, maxlen = maxlen_seq, padding = 'post', truncating='post')
test_q8s_data = tokenizer_encoder_q8s.texts_to_sequences(test_q8s_grams)
test_q8s_data = sequence.pad_sequences(test_q8s_data, maxlen = maxlen_seq, padding = 'post', truncating='post')

# extra MSA features 
msa_dim = msas[0].shape[1]
msas_padded = np.zeros([len(msas), maxlen_seq, msa_dim], dtype=np.float32)
for i in range(len(msas)):
    msas_padded[i, :msas[i].shape[0], :] = msas[i]

X_test = tf.convert_to_tensor(test_input_data, dtype=tf.float32)
X_test_q8s = tf.convert_to_tensor(test_q8s_data, dtype=tf.float32)
X_test_msa = tf.convert_to_tensor(msas_padded, dtype=tf.float32)
X_test_seqlen = length_aas



# prediction
angle_scale = 180 / np.pi 

imported = tf.saved_model.load(model_dir)
# loading a model from tf.saved_model.save() gives a forward pass function
model_func = imported.signatures['serving_default']
model_out = model_func(input_1=X_test, input_2=X_test_q8s, input_3=X_test_msa)
torsion_angles = model_out['torsion_angles']
metric = model_out['time_distributed']
#metric = model_out['metric']
phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]  
phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
#dist_matrix, coordinates = DistanceMatrix()(torsion_angles)
dist_matrix = coordinates_to_dist_matrix(metric)

phi_scaled, psi_scaled = phi_scaled.numpy(), psi_scaled.numpy()
dist_matrix = dist_matrix.numpy()

#phi_output, psi_output, dist_matrix_output = [], [], []
output = {}
for i in range(X_test.shape[0]):
    #phi_output.append(phi_scaled[i, :X_test_seqlen[i]])
    #psi_output.append(psi_scaled[i, :X_test_seqlen[i]])
    #dist_matrix_output.append(dist_matrix[i, :X_test_seqlen[i], :X_test_seqlen[i]])
    output[pdbs[i]] = {"dcalpha": dist_matrix[i, :X_test_seqlen[i], :X_test_seqlen[i]], "psi": psi_scaled[i, :X_test_seqlen[i]], "phi": phi_scaled[i, :X_test_seqlen[i]]}

#output = [dist_matrix_output, psi_output, phi_output]
with open(os.path.join(model_dir, 'predictions.pkl'), 'wb') as handle:
    pickle.dump(output, handle)

