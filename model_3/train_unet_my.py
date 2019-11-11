import os, sys, time 
import pickle
import numpy as np
import tensorflow as tf
from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.layers import *
from model_utils import *
from datetime import datetime

#######################################
# require tensorflow 2.0 
# CPU-only: pip install tf-nightly-2.0-preview
# GPU: pip install tf-nightly-2.0-gpu-preview
print(tf.__version__)
#######################################

# python train_unet.py

model_name = os.path.basename(sys.argv[0]).split(".")[0]   # or model name
if len(sys.argv) > 1:
    model_name = model_name + sys.argv[1]
model_name = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + model_name

#######################################
# Load data 
#######################################

# replace this with your data directory
'''import os
os.environ['CUDA_VISIBLE_DEVICES'] = '2' 
data_dir = "./data"

folds_use = [1,2,3,4,5,6,7,8,9,10]
data_fields = [[] for _ in range(9)]
for i in folds_use:
    with open(os.path.join(data_dir, 'train_fold_{}.pkl'.format(i)), 'rb') as f:
        data = pickle.load(f)
        for j in range(len(data_fields)):
            data_fields[j].append(data[j])
'''
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '3'      
train_num = 28000
data_dir = "./data100"

folds_use = [i for i in range(1,105)]
data_fields = [[] for _ in range(9)]
for i in folds_use:
    with open(os.path.join(data_dir, 'training_100_{}.pkl'.format(i)), 'rb') as f:
        data = pickle.load(f)
        for j in range(len(data_fields)):
            data_fields[j].append(data[j])

for j in range(len(data_fields)):
    data_fields[j] = np.concatenate(data_fields[j])
indices, pdbs, length_aas, pdb_aas, q8s, dcalphas, psis, phis, msas = data_fields
print("Total number of protein sequences:", len(pdb_aas))

# filter out sequences that are too long 
maxlen_seq = 384
minlen_seq = 10
lenth_mask = (length_aas < maxlen_seq) & (length_aas > minlen_seq)
for j in range(len(data_fields)):
    data_fields[j] = data_fields[j][lenth_mask]
indices, pdbs, length_aas, pdb_aas, q8s, dcalphas, psis, phis, msas = data_fields
msas = [np.stack(x).transpose().astype(np.float32) for x in msas]
phis = [np.array(x) for x in phis]
psis = [np.array(x) for x in psis]
print("Number of protein sequences shorter than {}: {}".format(maxlen_seq, len(pdb_aas)))

#############################################
# Prepare dataset for traininng 
#############################################

# Pad target distance matrix to the same size 
# Mask the padded regions later so that they don't contribute to the loss
dcalphas_pad = np.zeros((len(dcalphas), maxlen_seq, maxlen_seq), dtype=np.float32)
for i in range(len(dcalphas)):
    length = dcalphas[i].shape[0]
    dcalphas_pad[i, :length, :length] = dcalphas[i]

# Pad target torsion angles to the same size 
train_target_phis = np.zeros([len(phis), maxlen_seq], dtype=np.float32)
for i in range(len(phis)):
    train_target_phis[i, :phis[i].shape[0]] = phis[i]
train_target_psis = np.zeros([len(psis), maxlen_seq], dtype=np.float32)
for i in range(len(psis)):
    train_target_psis[i, :psis[i].shape[0]] = psis[i]

# Input data 
train_input_seqs = pdb_aas
train_input_grams = seq2ngrams(train_input_seqs, n=1)
train_q8s = q8s
train_q8s_grams = seq2ngrams(train_q8s, n=1)

# Define tokenizer encoder the input sequence
tokenizer_encoder = text.Tokenizer()
tokenizer_encoder.fit_on_texts(train_input_grams)
tokenizer_encoder_q8s = text.Tokenizer()
tokenizer_encoder_q8s.fit_on_texts(train_q8s_grams)

# Tokenize the input sequences for use in training
train_input_data = tokenizer_encoder.texts_to_sequences(train_input_grams)
train_input_data = sequence.pad_sequences(train_input_data, maxlen = maxlen_seq, padding = 'post', truncating='post')
train_q8s_data = tokenizer_encoder_q8s.texts_to_sequences(train_q8s_grams)
train_q8s_data = sequence.pad_sequences(train_q8s_data, maxlen = maxlen_seq, padding = 'post', truncating='post')

# The number of words and tags to be passed as parameters to the model
n_words = len(tokenizer_encoder.word_index) + 1
n_words_q8s = len(tokenizer_encoder_q8s.word_index) + 1

# extra MSA features 
msa_dim = msas[0].shape[1]
msas_padded = np.zeros([len(msas), maxlen_seq, msa_dim], dtype=np.float32)
for i in range(len(msas)):
    msas_padded[i, :msas[i].shape[0], :] = msas[i]

# Train/validation set split 
training_idx = np.arange(train_num)
validation_idx = np.arange(train_num, len(pdb_aas))

X_train = train_input_data[training_idx]
X_val = train_input_data[validation_idx]

X_train_q8s = train_q8s_data[training_idx]
X_val_q8s = train_q8s_data[validation_idx]

X_train_msa = msas_padded[training_idx]
X_val_msa = msas_padded[validation_idx]

X_train_seqlen = length_aas[training_idx]
X_val_seqlen = length_aas[validation_idx]

y_train_dist_matrix = dcalphas_pad[training_idx]
y_val_dist_matrix = dcalphas_pad[validation_idx]

y_train_phis = train_target_phis[training_idx]
y_val_phis = train_target_phis[validation_idx]
y_train_psis = train_target_psis[training_idx]
y_val_psis = train_target_psis[validation_idx]

pdb_names_train = pdbs[training_idx]
pdb_names_val = pdbs[validation_idx]


#############################################
# Create a directory to save model checkpoints 
#############################################

if not os.path.exists("logs"):
    os.mkdir("logs")
log_dir = 'logs/{}'.format(model_name)
os.mkdir(log_dir)

with open(os.path.join(log_dir, 'tokenizer_encoder.pickle'), 'wb') as handle:
    pickle.dump(tokenizer_encoder, handle)

with open(os.path.join(log_dir, 'tokenizer_encoder_q8s.pickle'), 'wb') as handle:
    pickle.dump(tokenizer_encoder_q8s, handle)


#############################################
# build a neural network model
#############################################

# Here we use layer functions to build neural network graph (as opposed to defining a Model class)

# Input layer: (None, ) shape allows the network to be flexible to input sequence length 
# inputs: [protein_sequence, q8_sequence, MSA features]
inputs = [Input(shape = (None, )), Input(shape = (None, )), Input(shape = (None, 21))]

# Word (amino acid) embedding layers 
embedding1 = Embedding(input_dim = n_words, output_dim = 128, input_length = None, input_shape=(None,))
embedding2 = Embedding(input_dim = n_words_q8s, output_dim = 64, input_length = None, input_shape=(None,))
#embedding2 = Embedding(input_dim = n_words, output_dim = 64, input_length = None, input_shape=(None,))
embed_input = embedding1(inputs[0]) # raw seq
embed_q8s = embedding2(inputs[1]) # q8

# Concatenate different types of inputs 
merged_input = concatenate([embed_input, embed_q8s, inputs[2]], axis = 2)

#### bidirectional RNN
# rnn = Bidirectional(GRU(units=64, return_sequences=True, recurrent_dropout=0.1, recurrent_activation='relu'))
# y = rnn(merged_input)

#### A seq-to-seq 1-D convolutional U-Net (https://arxiv.org/abs/1505.04597)
y = Conv_UNet(merged_input, droprate=0.4)
y2 = Conv_UNet(merged_input, droprate=0.2)
#y3 = Conv_UNet(merged_input, droprate=0.4)

#### Torsion angle prediction
# The model output 3 angles (phi, psi, omega) for computing distance matrix. 
# But we only use phi, psi for evaluating angles since omega angle is usually fixed in the protein.
metric = TimeDistributed(Dense(64, activation = "linear"))(y2)
#m2 = TimeDistributed(Dense(64, activation = "linear"))(y2)
#metric = tf.multiply(metric + 0.1 * m2, 1., name = 'metric')
#ya = TimeDistributed(Dense(2, activation = "tanh"))(y)
#yc = TimeDistributed(Dense(2, activation = "tanh"))(y)
y = TimeDistributed(Dense(2, activation = "tanh"))(y)
#yb = TimeDistributed(Dense(2, activation = "tanh"))(y3)
#yd = TimeDistributed(Dense(2, activation = "tanh"))(y3)
#y3 = TimeDistributed(Dense(2, activation = "tanh"))(y3)
#angles = tf.multiply(y + ya * 0.5 + 0.1 * y3 + yb * 0.05 + yc * 0.01 + yd * 0.005, np.pi, name="torsion_angles")
#angles2 = tf.multiply(y3, np.pi, name="torsion_angles2")
#y = TimeDistributed(Dense(360, activation = "tanh"))(y)
#y3 = TimeDistributed(Dense(360, activation = "tanh"))(y3)
angles = tf.multiply(y, np.pi, name="torsion_angles")
#angles2 = tf.multiply(y3, np.pi, name="torsion_angles2")


# The RGN paper predict torsion angles using softmax probabilities over a learned alphabet/mixture of angles (https://www.biorxiv.org/content/biorxiv/early/2018/08/29/265231.full-text.pdf)
# y = TimeDistributed(Dense(50, activation = "softmax"))(y)
# angles = TorsionAngles(alphabet_size=50)(y)

#### Build Model. 
# Define the computational graph using model input and model output
model = tf.keras.Model(inputs, [angles, metric])
model.summary()

##################################################
# Alternatively, try soemthing else (e.g. Attention, Transformer, etc)
##################################################


##################################################
# Train the model
##################################################

#batch_size = 32
batch_size = 128
angle_scale = 180 / np.pi  # convert angles from [-pi, pi] to [-180, 180]
n_iter = int(X_train.shape[0] / batch_size)

#### For RNN (slower)
# n_epochs = 8
# lr_decay_iters = [n_iter * 3.0, n_iter * 6.0]  # lr decay at n_iters * epoch_steps
# lr_steps = [0.001, 0.0005, 0.0001]  # learning rate decay steps

#### For ConvNet
n_epochs = 20
lr_decay_iters = [n_iter * 5.0, n_iter * 8.0]  # lr decay at n_iters * epoch_steps
lr_steps = [0.0001, 0.00005, 0.00002]  # learning rate decay steps

learning_rate_fn = tf.keras.optimizers.schedules.PiecewiseConstantDecay(lr_decay_iters, lr_steps)
optimizer = tf.keras.optimizers.RMSprop(learning_rate=learning_rate_fn)

global_step = tf.Variable(0, trainable=False)
train_loss_history = []
val_loss_history = []

torsion_loss_weight = 1.0 / 50
dist_loss_weight = 1.0 

for epoch in range(1, n_epochs+1):
    
    ################################### 
    # Train
    ################################### 
    idx_shuffle = np.random.permutation(X_train.shape[0])

    train_rmsd_dist = tf.metrics.Mean()
    train_rmsd_dist_norm = tf.metrics.Mean()
    train_rmsd_angle = tf.metrics.Mean()
    train_rmsd_all = tf.metrics.Mean()

    # Iterate through mini-batchs
    for it in range(n_iter):
        # Shuffle (alternatively, use tf.dataset)
        idx_batch = idx_shuffle[it * batch_size : (it+1) * batch_size]
        batch_seq = tf.convert_to_tensor(X_train[idx_batch])
        batch_q8s = tf.convert_to_tensor(X_train_q8s[idx_batch])
        batch_seqlen = tf.convert_to_tensor(X_train_seqlen[idx_batch], dtype=tf.float32)
        batch_dcalphas = tf.convert_to_tensor(y_train_dist_matrix[idx_batch])
        batch_msa = tf.convert_to_tensor(X_train_msa[idx_batch])

        # Compute loss 
        with tf.GradientTape() as tape:
            torsion_angles, metric = model([batch_seq, batch_q8s, batch_msa])
            phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]
            #torsion_angles, torsion_angles2, metric = model([batch_seq, batch_q8s, batch_msa])
            #phi, psi = torsion_angles[:, :, 0], torsion_angles2[:, :, 0]
            phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
            loss_phi_batch = rmsd_torsion_angle(phi_scaled, y_train_phis[idx_batch], batch_seqlen)
            loss_psi_batch = rmsd_torsion_angle(psi_scaled, y_train_psis[idx_batch], batch_seqlen)
            loss_phi = tf.reduce_mean(loss_phi_batch)  
            loss_psi = tf.reduce_mean(loss_psi_batch)

            #dist_matrix, coordinates = DistanceMatrix()(torsion_angles)
            dist_matrix = coordinates_to_dist_matrix(metric)
            loss_drmsd_batch = drmsd_dist_matrix(dist_matrix, batch_dcalphas, batch_seqlen)
            loss_drmsd = tf.reduce_mean(loss_drmsd_batch)  # drmsd metric
            loss_drmsd_normalized = tf.reduce_mean(loss_drmsd_batch / tf.sqrt(batch_seqlen))  # longer proteins have larger distance
            
            # the optimization objective can be 1. loss_drmsd or 2. (loss_phi + loss_psi), or both
            # optimizing distance mattrix (drmsd) and tortion angles separately probably gives better rerults. 
            #loss_all = loss_drmsd   # distance matrix loss
            # loss_all = loss_phi + loss_psi  # torsion angle loss
            loss_all = loss_drmsd * dist_loss_weight + (loss_phi + loss_psi) * torsion_loss_weight
            

        # Compute gradient 
        grads = tape.gradient(loss_all, model.trainable_variables)  # loss includ both distance matrix and torsion angles 
        # grads = tape.gradient(loss_all, model.trainable_variables)
        grads, global_norm = tf.clip_by_global_norm(grads, 5.0)
        
        # Backprop
        optimizer.apply_gradients(zip(grads, model.trainable_variables), global_step)

        # Record metrics
        train_rmsd_dist_norm(loss_drmsd_normalized)
        train_rmsd_dist(loss_drmsd)
        train_rmsd_angle((loss_phi + loss_psi)/2)
        train_rmsd_all(loss_all)

        if it % 10 == 0:
            print("Epoch {:04d} Batch {:03d}/{:03d}: Loss: {:.5g}, dist RMSD: {:.5g}, (length-normalized): {:.5g}, angle RMSD: {:.5g}".format(
                epoch, it, n_iter, loss_all, loss_drmsd, loss_drmsd_normalized, (loss_phi + loss_psi)/2))
        
    
    ###################################  
    # Validation 
    ################################### 
    idx_val = np.arange(X_val.shape[0])
    n_iter_val = int(len(idx_val) / batch_size) + 1  # use all

    val_rmsd_dist = tf.metrics.Mean()
    val_rmsd_dist_norm = tf.metrics.Mean()
    val_rmsd_angle = tf.metrics.Mean()
    val_rmsd_all = tf.metrics.Mean()

    for it in range(n_iter_val):
        idx_batch = idx_val[it * batch_size : (it+1) * batch_size]
        batch_seq = tf.convert_to_tensor(X_val[idx_batch])
        batch_q8s = tf.convert_to_tensor(X_val_q8s[idx_batch])
        batch_seqlen = tf.convert_to_tensor(X_val_seqlen[idx_batch], dtype=tf.float32)
        batch_dcalphas = tf.convert_to_tensor(y_val_dist_matrix[idx_batch])
        batch_msa = tf.convert_to_tensor(X_val_msa[idx_batch])

        with tf.GradientTape() as tape:
            torsion_angles, metric = model([batch_seq, batch_q8s, batch_msa])
            phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]
            #torsion_angles, torsion_angles2, metric = model([batch_seq, batch_q8s, batch_msa])
            #phi, psi = torsion_angles[:, :, 0], torsion_angles2[:, :, 0]
            phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
            loss_phi_batch = rmsd_torsion_angle(phi_scaled, y_val_phis[idx_batch], batch_seqlen)
            loss_psi_batch = rmsd_torsion_angle(psi_scaled, y_val_psis[idx_batch], batch_seqlen)
            loss_phi = tf.reduce_mean(loss_phi_batch)  
            loss_psi = tf.reduce_mean(loss_psi_batch)

            # torsion_angles are on the scale of [-3.14, 3.14], the target labels are on the scale of [-180, 180]. Need to convert after prediction
            #dist_matrix, coordinates = DistanceMatrix()(torsion_angles)
            dist_matrix = coordinates_to_dist_matrix(metric)
            loss_drmsd_batch = drmsd_dist_matrix(dist_matrix, batch_dcalphas, batch_seqlen)
            loss_drmsd = tf.reduce_mean(loss_drmsd_batch)
            loss_drmsd_normalized = tf.reduce_mean(loss_drmsd_batch / tf.sqrt(batch_seqlen))

            #loss_all = loss_drmsd
            # loss_all = loss_phi + loss_psi
            loss_all = loss_drmsd * dist_loss_weight + (loss_phi + loss_psi) * torsion_loss_weight

        # no gradient descent during validation
        val_rmsd_dist_norm(loss_drmsd_normalized)
        val_rmsd_dist(loss_drmsd)
        val_rmsd_angle((loss_phi + loss_psi)/2)
        val_rmsd_all(loss_all)

        if it == 0:
            plot_pred = dist_matrix[:8].numpy()
            plot_gt = batch_dcalphas[:8].numpy()
            plot_names = pdb_names_val[idx_batch][:8]
            plot_len = batch_seqlen[:8].numpy().astype(int)
            plot_rmds = loss_drmsd_batch[:8].numpy()
            plot_dist_matrix(plot_pred, plot_gt, plot_names, plot_len, plot_rmds,
                os.path.join(log_dir, "plot_distance_matrix_epoch-{}.pdf".format(epoch)))

    ################################### 

    # plot history if not using tensorboard
    train_loss_history.append([train_rmsd_all.result(), train_rmsd_dist.result(), train_rmsd_dist_norm.result(), train_rmsd_angle.result()])
    val_loss_history.append([val_rmsd_all.result(), val_rmsd_dist.result(), val_rmsd_dist_norm.result(), val_rmsd_angle.result()])

    print("\n===========================================\n")
    print("Epoch {:03d}: Train Loss: {:.5g}, dist RMSD: {:.5g}, (length-normalized): {:.5g}, angle RMSD: {:.5g} \n ".format(
        epoch, *train_loss_history[-1]) + 
        "             Val Loss: {:.5g}, dist RMSD: {:.5g}, (length-normalized): {:.5g}, angle RMSD: {:.5g}".format(
        *val_loss_history[-1]))
    print("\n===========================================\n")

    # save model checkpoints (this overwrites previous checkpoints)
    tf.saved_model.save(model, log_dir+str(epoch))

    # plot
    train_loss_arr = np.array(train_loss_history).transpose()
    val_loss_arr = np.array(val_loss_history).transpose()
    plot_train_val(train_loss_arr[0], val_loss_arr[0], 
            title="RMSD (dist + angle) loss", savepath=os.path.join(log_dir, "rmsd_all.pdf"))
    plot_train_val(train_loss_arr[1], val_loss_arr[1], 
            title="RMSD (distance matrix) loss", savepath=os.path.join(log_dir, "rmsd_distance_matrix.pdf"))
    plot_train_val(train_loss_arr[2], val_loss_arr[2], 
            title="RMSD (distance matrix, lengh-normalized) loss", savepath=os.path.join(log_dir, "rmsd_distance_matrix_normalized.pdf"))
    plot_train_val(train_loss_arr[3], val_loss_arr[3], 
            title="RMSD (torsion angles) loss", savepath=os.path.join(log_dir, "rmsd_torsion_angles.pdf"))
    



