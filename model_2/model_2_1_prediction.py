import os, sys, time 
import pickle
import numpy as np
import tensorflow as tf
from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.layers import *

#!pip install tf-nightly-gpu-2.0-preview
print("TensorFlow version:", tf.__version__)

# fetch the util.py file for future import
print('fetching the util.py file...')
sys.path.append('..')
from utils.rgn_utils import *

print('prediction on the testing set...')

# define data loading function
test_data_path = "cuprotein/test/cuprotein_testing.pkl"
data = pickle.load(open(test_data_path, 'rb'))
test_data_list = list(data.items())
# get the data fields...
pdbs = [xx[0] for xx in test_data_list]
lengths = [xx[1]["length"] for xx in test_data_list]
aa = [xx[1]["aa"] for xx in test_data_list]
ss = [xx[1]["ss"] for xx in test_data_list]
chain = [xx[1]["chain"] for xx in test_data_list]
dcalphas = [xx[1]["dcalpha"] for xx in test_data_list]
coords = [xx[1]["coords"] for xx in test_data_list]
psis = [xx[1]['psi'] for xx in test_data_list]
phis = [xx[1]['phi'] for xx in test_data_list]
msacovonehot = [np.array(xx[1]['msacovonehot'], dtype=np.float32) for xx in test_data_list]
msacovembedding = [np.array(xx[1]['msacovembedding'], dtype=np.float32) for xx in test_data_list]

def load_test_data(test_data_dir):
    pdbs = []; lengths = []; aa = []; ss = []; dcalphas = []; coords = []; psis = []; phis = []; pssms = []
    for fold in [1]:
        # load current pickle file
        dict_curr_file = pickle.load(open(test_data_path), 'rb'))
        # read current file data and put them into data arrays
        items_curr_file = list(dict_curr_file.items())
        pdbs.extend([item[0] for item in items_curr_file])
        lengths.extend([item[1]['length'] for item in items_curr_file])
        aa.extend([item[1]['aa'] for item in items_curr_file])
        ss.extend([item[1]['ss'] for item in items_curr_file])
        dcalphas.extend([item[1]['dcalpha'] for item in items_curr_file])
        coords.extend([item[1]['coords'] for item in items_curr_file])
        psis.extend([item[1]['psi'] for item in items_curr_file])
        phis.extend([item[1]['phi'] for item in items_curr_file])
        pssms.extend([item[1]['pssm'] for item in items_curr_file])
    # return the 9 feature data arrays
    return [np.array(pdbs), np.array(lengths), np.array(aa), np.array(ss), np.array(dcalphas), np.array(coords), np.array(psis), np.array(phis), np.array(pssms)]

# load test data
model_name = 'bigru-skip-11-june'
#data_dir = '/gdrive/My Drive/CU-TSP/data/100k/train/'
#log_dir = '/gdrive/My Drive/CU-TSP/models_code/wf2232-yg2541/' + 'logs/{}'.format(model_name)
test_data_dir = 'cuprotein/test'
model_dir = 'log/{}'.format(model_name)  # by default, fetch the model from the log_dir
pdbs_test, lengths_test, aa_test, ss_test, dcalphas_test, coords_test, psis_test, phis_test, pssms_test = load_test_data(test_data_dir)

pssms_test = [np.stack(x).transpose().astype(np.float32) for x in pssms_test] # pssms_test[i].shape = (dimension lengths_test[i],21)

# maxlen_seq = np.max(lengths_test) 
maxlen_seq = 384

# Input data 
test_input_seqs = aa_test
test_input_grams = seq2ngrams(test_input_seqs, n=1)
test_q8s = ss_test
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

# extra PSSM features 
pssms_dim = pssms_test[0].shape[1] # 21
pssms_test_padded = np.zeros([len(pssms_test), maxlen_seq, pssms_dim], dtype=np.float32)
for i in range(len(pssms_test)):
    truncated_length = np.min([maxlen_seq, lengths_test[i]])
    pssms_test_padded[i, :truncated_length, :] = pssms_test[i][:truncated_length, :]
    
msacovonehot_pad = np.zeros((len(msacovonehot), maxlen_seq, maxlen_seq), dtype=np.float32)
for i in range(len(msacovonehot)):
    protein_length = msacovonehot[i].shape[0]
    msacovonehot_pad[i, :protein_length, :protein_length] = msacovonehot[i]
    
msacovembedding_pad = np.zeros((len(msacovembedding), maxlen_seq, maxlen_seq), dtype=np.float32)
for i in range(len(msacovembedding)):
    protein_length = msacovembedding[i].shape[0]
    msacovembedding_pad[i, :protein_length, :protein_length] = msacovembedding[i]

msacovonehot_pad = msacovonehot_pad / 10000
msacovembedding_pad = msacovembedding_pad / 260000

X_test = tf.convert_to_tensor(test_input_data, dtype=tf.float32)
X_test_q8s = tf.convert_to_tensor(test_q8s_data, dtype=tf.float32)
X_test_pssm = tf.convert_to_tensor(pssms_test_padded, dtype=tf.float32)
X_test_oh = tf.convert_to_tensor(msacovonehot_pad, dtype=tf.float32)
X_test_eb = tf.convert_to_tensor(msacovembedding_pad, dtype=tf.float32)
X_test_seqlen = lengths_test

# prediction
angle_scale = 180 / np.pi 

imported = tf.saved_model.load(model_dir)
# loading a model from tf.saved_model.save() gives a forward pass function
model_func = imported.signatures['serving_default']
model_out = model_func(input_1=X_test, input_2=X_test_q8s, input_3=X_test_pssm, input_4=X_test_oh, input_5=X_test_eb )
# model_out = model_func(input_28 = X_test, input_29 = X_test_q8s, input_30 = X_test_msa)
torsion_angles = model_out['torsion_angles']
coordinates = model_out['coordinates']
coordinates = coordinates.numpy()

phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]  
phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
dist_matrix = pairwise_distance_self(coordinates)

phi_scaled, psi_scaled = phi_scaled.numpy(), psi_scaled.numpy()
dist_matrix = dist_matrix.numpy()

output_dict = {}
for i in range(len(pdbs_test)):
    protein_id = pdbs_test[i]
    protein_length = lengths_test[i]
    phi = phi_scaled[i, :protein_length]
    psi = psi_scaled[i, :protein_length]
    coords = coordinates[i, :protein_length, :]
    dcalpha = dist_matrix[i, :protein_length, :protein_length]
    #print(coords.shape, ', ', sep = '', end='')
    #print(protein_id, protein_length, coords.shape, dcalpha.shape)
    output_dict[protein_id] = {"length": protein_length, "phi": phi, "psi": psi, "coords": coords, "dcalpha": dcalpha}

input_dict = pickle.load(open(os.path.join(test_data_dir, 'testing.pkl'), 'rb'))
#print(list(input_dict.keys()))
    
result_dir = 'results/model_2_9_jul_output_all.pkl'
pickle.dump(output_dict, open(result_dir, "wb"))
