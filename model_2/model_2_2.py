# run on 100k dataset with and without q8s

#######################################
# REQUIRE tensorflow 2.0 (nightly GPU version)
# CPU-only: pip install tf-nightly-2.0-preview
# GPU: pip install tf-nightly-2.0-gpu-preview
#######################################

#!pip install tf-nightly-gpu-2.0-preview==2.0.0-dev20190227
import tensorflow as tf
print("TensorFlow version:", tf.__version__)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# set model name, which determines the log directory (to store checkpoints and trained model instances)
import os, sys, time 
import pickle
import numpy as np
import tensorflow as tf
from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.layers import *
from datetime import datetime
sys.path.append('..')
from utils.rgn_utils import *

# model 2 (among models 1-4), officially
model_name = "bigru-skip-9-jul"

# may need to modify the paths accordingly
train_data_dir = "cuprotein/train/"
val_data_dir = "cuprotein/validation/"
log_dir = os.path.join("./","logs/{}".format(model_name))

# split data into super-folds
super_folds = [list(range(1, 22)), list(range(22, 43)), list(range(43, 64)), list(range(64, 85)), list(range(85,105))]
starting_epochs = [1, 16, 31, 46, 51]
n_epochs = 10 # per sub-epoch

for super_idx in range(5):
    folds = super_folds[super_idx]
    print("============ training data folds_used ============\n", folds)

    # load training data from pickled dictionary
    data_list = []
    maxlen_seq, minlen_seq = 384, 30
    ff = lambda xx: (xx[1]['length'] <= maxlen_seq and xx[1]['length'] >= minlen_seq)
    for fold in folds:
        data = pickle.load(open(os.path.join(train_data_dir, 'cuprotein_{}.pkl'.format(fold)), 'rb'))
        data_list += [xx for xx in data.items() if ff(xx)]

    # count...
    num_train = len(data_list)

    # load validation data and append it to the training data
    data_val = pickle.load(open(os.path.join(val_data_dir, 'cuprotein_validation.pkl'), 'rb'))
    data_val_list = [xx for xx in data_val.items() if ff(xx)]
    num_val = len(data_val_list)

    # put them together for padding, reshaping, type casting, etc. 
    data_list += data_val_list
    print("min and max length: {}, {}".format(minlen_seq, maxlen_seq))
    print("total number of training and validation used: {}, {}".format(num_train, num_val))

    # put data into lists and cast into numpy arrays whenever necessary
    pdbs = np.array([xx[0] for xx in data_list])
    aa = [xx[1]['aa'] for xx in data_list]
    ss = [xx[1]['ss'] for xx in data_list]
    lengths = np.array([xx[1]['length'] for xx in data_list])
    coords = [np.array(xx[1]['coords'], dtype=np.float32) for xx in data_list]
    dcalphas = [np.array(xx[1]['dcalpha'], dtype=np.float32) for xx in data_list]
    msacovonehot = [np.array(xx[1]['msacovonehot'], dtype=np.float32) for xx in data_list]
    msacovembedding = [np.array(xx[1]['msacovembedding'], dtype=np.float32) for xx in data_list]
    phis = [np.array(xx[1]['phi'], dtype=np.float32) for xx in data_list]
    psis = [np.array(xx[1]['psi'], dtype=np.float32) for xx in data_list]
    pssms = [np.stack(xx[1]['pssm']).transpose().astype(np.float32) for xx in data_list] # pssm[i].shape = (21, len[i])

    # pad dcalpha
    dcalphas_pad = np.zeros((len(dcalphas), maxlen_seq, maxlen_seq), dtype=np.float32)
    for i in range(len(dcalphas)):
        protein_length = dcalphas[i].shape[0]
        dcalphas_pad[i, :protein_length, :protein_length] = dcalphas[i]
        
    # pad pssm
    pssm_dim = pssms[0].shape[1] # 21
    pssms_padded = np.zeros([len(pssms), maxlen_seq, pssm_dim], dtype=np.float32)
    for i in range(len(pssms)):
        pssms_padded[i, :pssms[i].shape[0], :] = pssms[i]
        
    # pad torsion angles
    train_target_phis = np.zeros([len(phis), maxlen_seq], dtype=np.float32)
    train_target_psis = np.zeros([len(psis), maxlen_seq], dtype=np.float32)
    for i in range(len(phis)):
        train_target_phis[i, :len(phis[i])] = phis[i]
        train_target_psis[i, :len(psis[i])] = psis[i]
        
    # input data 
    train_input_seqs = aa
    train_input_grams = seq2ngrams(train_input_seqs, n=1)
    train_q8s = ss
    train_q8s_grams = seq2ngrams(train_q8s, n=1)

    # define tokenizer encoder the input sequence
    tokenizer_encoder = text.Tokenizer()
    tokenizer_encoder.fit_on_texts(train_input_grams)
    tokenizer_encoder_q8s = text.Tokenizer()
    tokenizer_encoder_q8s.fit_on_texts(train_q8s_grams)

    # tokenize the input sequences for use in training
    train_input_data = tokenizer_encoder.texts_to_sequences(train_input_grams)
    train_input_data = sequence.pad_sequences(train_input_data, maxlen = maxlen_seq, padding = 'post', truncating='post')
    train_q8s_data = tokenizer_encoder_q8s.texts_to_sequences(train_q8s_grams)
    train_q8s_data = sequence.pad_sequences(train_q8s_data, maxlen = maxlen_seq, padding = 'post', truncating='post')

    # number of words and tags to be passed as parameters to the model
    n_words = len(tokenizer_encoder.word_index) + 1
    n_words_q8s = len(tokenizer_encoder_q8s.word_index) + 1

    # make dir for saving training logs...create tokenizer_encoder pickles when it's the first time
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)
    with open(os.path.join(log_dir, 'tokenizer_encoder.pickle'), 'wb') as handle:
        pickle.dump(tokenizer_encoder, handle)

    with open(os.path.join(log_dir, 'tokenizer_encoder_q8s.pickle'), 'wb') as handle:
        pickle.dump(tokenizer_encoder_q8s, handle)

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

    ###### split train and val, prepare for training ############
    # split train and val data

    training_idx = np.arange(num_train)
    validation_idx = np.arange(num_train, num_train+num_val)
    training_idx[0], training_idx[-1], validation_idx[0], validation_idx[-1], num_train+num_val

    X_train = train_input_data[training_idx]
    X_val = train_input_data[validation_idx]

    X_train_q8s = train_q8s_data[training_idx]
    X_val_q8s = train_q8s_data[validation_idx]

    X_train_pssm = pssms_padded[training_idx] # padded
    X_val_pssm = pssms_padded[validation_idx]

    X_train_seqlen = lengths[training_idx]
    X_val_seqlen = lengths[validation_idx]

    y_train_dist_matrix = dcalphas_pad[training_idx]
    y_val_dist_matrix = dcalphas_pad[validation_idx]

    y_train_phis = train_target_phis[training_idx]
    y_val_phis = train_target_phis[validation_idx]
    y_train_psis = train_target_psis[training_idx]
    y_val_psis = train_target_psis[validation_idx]

    pdb_names_train = pdbs[training_idx]
    pdb_names_val = pdbs[validation_idx]

    X_train_msa_oh = msacovonehot_pad[training_idx] # padded
    X_val_msa_oh = msacovonehot_pad[validation_idx]

    X_train_msa_eb = msacovembedding_pad[training_idx] # padded
    X_val_msa_eb = msacovembedding_pad[validation_idx]

    print('train and validate sizes: {}, {}'.format(len(y_train_dist_matrix), len(y_val_dist_matrix)))


    ############################################
    # build a neural network model
    #############################################

    # Here we use layer functions to build neural network graph (as opposed to defining a Model class)

    # Input layer: (None, ) shape allows the network to be flexible to input sequence length 
    # inputs: [protein_sequence, q8_sequence, MSA features]
    inputs = [Input(shape = (None, )), Input(shape = (None, )), Input(shape = (None, 21)), Input(shape = (None, maxlen_seq)), Input(shape = (None, maxlen_seq))]

    # Word (amino acid) embedding layers 
    embedding1 = Embedding(input_dim = n_words, output_dim = 128, input_length = None, input_shape=(None,))
    embedding2 = Embedding(input_dim = n_words_q8s, output_dim = 64, input_length = None, input_shape=(None,))
    embed_input = embedding1(inputs[0]) # raw seq
    embed_q8s = embedding2(inputs[1]) # q8

    vec_oh = TimeDistributed(Dense(64, activation = "relu"))(inputs[3])
    vec_eb = TimeDistributed(Dense(64, activation = "relu"))(inputs[4])

    # Concatenate different types of inputs 
    merged_input = concatenate([embed_input, embed_q8s, inputs[2], vec_oh, vec_eb], axis = 2)

    mid1 = Bidirectional(tf.compat.v1.keras.layers.CuDNNGRU(units=512, return_sequences=True))(merged_input)
    mid2 = Bidirectional(tf.compat.v1.keras.layers.CuDNNGRU(units=256, return_sequences=True))(mid1)
    mid3 = Bidirectional(tf.compat.v1.keras.layers.CuDNNGRU(units=128, return_sequences=True))(mid2)
    h1 = TimeDistributed(Dense(64, activation = "relu"))(mid3)
    mid = concatenate([mid1, mid2, mid3, h1], axis = -1)

    y = TimeDistributed(Dense(3, activation = "tanh"))(mid)
    angles = tf.multiply(y, np.pi, name="torsion_angles") 

    y2 = TimeDistributed(Dense(3))(mid)
    coords = tf.multiply(y2, 1, name="coordinates") 

    #### Build Model. 
    # Define the computational graph using model input and model output
    model = tf.keras.Model(inputs, [angles, coords])
    model.summary()

    global_step = tf.Variable(0, trainable=False)


    #################### Training #######################


    ##################################################
    # Train the model
    ##################################################

    batch_size = 32
    angle_scale = 180 / np.pi  # convert angles from [-pi, pi] to [-180, 180]
    n_iter = int(X_train.shape[0] / batch_size)

    print('Total number of epochs = {}'.format(n_epochs))
    print('Number of iterations per epoch = {}'.format(n_iter))
    lr_decay_iters = [n_iter * 20.0]  # lr decay at n_iters * epoch_steps
    lr_steps = [0.001, 0.0005]  # learning rate decay steps

    learning_rate_fn = tf.keras.optimizers.schedules.PiecewiseConstantDecay(lr_decay_iters, lr_steps)
    optimizer = tf.keras.optimizers.RMSprop(learning_rate=learning_rate_fn)

    train_loss_history = []
    val_loss_history = []

    checkpoint_directory = log_dir
    checkpoint_prefix = os.path.join(checkpoint_directory, "ckpt")

    checkpoint = tf.train.Checkpoint(optimizer=optimizer, model=model)
    status = checkpoint.restore(tf.train.latest_checkpoint(checkpoint_directory))

    starting_epoch = starting_epochs[super_idx] # 80 epochs in total
    for epoch in range(starting_epoch, n_epochs+starting_epoch):
        
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
            batch_pssm = tf.convert_to_tensor(X_train_pssm[idx_batch])
            batch_msa_oh = tf.convert_to_tensor(X_train_msa_oh[idx_batch])
            batch_msa_eb = tf.convert_to_tensor(X_train_msa_eb[idx_batch])


            # Compute loss 
            with tf.GradientTape() as tape:
                torsion_angles, coordinates = model([batch_seq, batch_q8s, batch_pssm, batch_msa_oh, batch_msa_eb])
                phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]  
                phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
                loss_phi_batch = rmsd_torsion_angle2(phi_scaled, y_train_phis[idx_batch], batch_seqlen)
                loss_psi_batch = rmsd_torsion_angle2(psi_scaled, y_train_psis[idx_batch], batch_seqlen)
                loss_phi = tf.reduce_sum(loss_phi_batch)  
                loss_psi = tf.reduce_sum(loss_psi_batch)

                dist_matrix = pairwise_distance_self(coordinates)
                width = min(10 * max(epoch, 1), maxlen_seq) # only calculating losses of entries that are within width off-diagonal
                loss_drmsd_batch = drmsd_dist_matrix2(dist_matrix, batch_dcalphas, batch_seqlen, width)
                loss_drmsd = tf.reduce_sum(loss_drmsd_batch)
                loss_drmsd_batch_full = drmsd_dist_matrix2(dist_matrix, batch_dcalphas, batch_seqlen)
                loss_drmsd_full = tf.reduce_sum(loss_drmsd_batch_full)
                
                batch_numbers = tf.reduce_sum(batch_seqlen)
                sq_batch_numbers = tf.reduce_sum(batch_seqlen * batch_seqlen)
                dist_numbers = tf. reduce_sum(batch_seqlen * (batch_seqlen-1) - tf.maximum(batch_seqlen-width,0) * tf.maximum(batch_seqlen-width-1,0))
        
                loss_phi_avg = tf.sqrt(loss_phi / batch_numbers)
                loss_psi_avg = tf.sqrt(loss_psi / batch_numbers)
                loss_main = 1/(2*batch_numbers+sq_batch_numbers)*(loss_phi+loss_psi) + sq_batch_numbers/dist_numbers/(2*batch_numbers+sq_batch_numbers)*loss_drmsd     
                loss_drmsd_avg = tf.sqrt(loss_drmsd_full / sq_batch_numbers) 
                loss_all = 1/(2*batch_numbers+sq_batch_numbers)*(loss_phi+loss_psi) + 1/(2*batch_numbers+sq_batch_numbers)*loss_drmsd_full     
                loss_root = tf.sqrt(loss_all)
                        

            # Compute gradient 
            grads = tape.gradient(loss_main, model.trainable_variables)
            grads, global_norm = tf.clip_by_global_norm(grads, 5.0)
            
            # Backprop
            optimizer.apply_gradients(zip(grads, model.trainable_variables), global_step)

            # Record metrics
    #         train_rmsd_dist_norm(loss_drmsd_normalized)
            train_rmsd_dist(loss_drmsd_avg)
            train_rmsd_angle((loss_phi_avg + loss_psi_avg)/2)
            train_rmsd_all(loss_root)

            if it % 20 == 0:
                print("Epoch {:04d} Batch {:03d}/{:03d}: Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g}".format(
                    epoch, it, n_iter, loss_root, loss_drmsd_avg, (loss_phi_avg + loss_psi_avg)/2))
                with open(os.path.join(log_dir, 'log.txt'), 'a') as f:
                    f.write("Epoch {:04d} Batch {:03d}/{:03d}: Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g}\n".format(
                        epoch, it, n_iter, loss_root, loss_drmsd_avg, (loss_phi_avg + loss_psi_avg)/2))
            
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
            batch_pssm = tf.convert_to_tensor(X_val_pssm[idx_batch])
            batch_msa_oh = tf.convert_to_tensor(X_val_msa_oh[idx_batch])
            batch_msa_eb = tf.convert_to_tensor(X_val_msa_eb[idx_batch])

            with tf.GradientTape() as tape:
                torsion_angles, coordinates = model([batch_seq, batch_q8s, batch_pssm, batch_msa_oh, batch_msa_eb])
                phi, psi = torsion_angles[:, :, 0], torsion_angles[:, :, 1]
                phi_scaled, psi_scaled = phi * angle_scale, psi * angle_scale  # from [-pi, pi] to [-180, 180]
                loss_phi_batch = rmsd_torsion_angle2(phi_scaled, y_val_phis[idx_batch], batch_seqlen)
                loss_psi_batch = rmsd_torsion_angle2(psi_scaled, y_val_psis[idx_batch], batch_seqlen)
                loss_phi = tf.reduce_sum(loss_phi_batch)  
                loss_psi = tf.reduce_sum(loss_psi_batch)

                dist_matrix = pairwise_distance_self(coordinates)
                loss_drmsd_batch = drmsd_dist_matrix2(dist_matrix, batch_dcalphas, batch_seqlen)
                loss_drmsd = tf.reduce_sum(loss_drmsd_batch)

                batch_numbers = tf.reduce_sum(batch_seqlen)
                sq_batch_numbers = tf.reduce_sum(batch_seqlen * batch_seqlen)
        
                loss_phi_avg = tf.sqrt(loss_phi / batch_numbers)
                loss_psi_avg = tf.sqrt(loss_psi / batch_numbers)
                loss_drmsd_avg = tf.sqrt(loss_drmsd / sq_batch_numbers) 
                loss_all = 1/(2*batch_numbers+sq_batch_numbers)*(loss_phi+loss_psi) + 1/(2*batch_numbers+sq_batch_numbers)*loss_drmsd     
                loss_root = tf.sqrt(loss_all)
                

            # no gradient descent during validation
            # val_rmsd_dist_norm(loss_drmsd_normalized)
            val_rmsd_dist(loss_drmsd_avg)
            val_rmsd_angle((loss_phi_avg + loss_psi_avg)/2)
            val_rmsd_all(loss_root)

            if it == 0:
                plot_pred = dist_matrix[:8].numpy()
                plot_gt = batch_dcalphas[:8].numpy()
                plot_names = pdb_names_val[idx_batch][:8]
                plot_len = batch_seqlen[:8].numpy().astype(int)
                plot_rmds = loss_drmsd_batch[:8].numpy()
                plot_dist_matrix(plot_pred, plot_gt, plot_names, plot_len, tf.sqrt(plot_rmds) / plot_len,
                    os.path.join(log_dir, "plot_distance_matrix_epoch-{}.pdf".format(epoch)))

        ################################### 

        # plot history if not using tensorboard
        train_loss_history.append([train_rmsd_all.result(), train_rmsd_dist.result(), train_rmsd_angle.result()])
        val_loss_history.append([val_rmsd_all.result(), val_rmsd_dist.result(), val_rmsd_angle.result()])

        print("\n===========================================\n")
        print("Epoch {:03d}: Train Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g} \n ".format(
            epoch, *train_loss_history[-1]) + 
            "             Val Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g}".format(
            *val_loss_history[-1]))
        print("\n===========================================\n")
        with open(os.path.join(log_dir, 'log.txt'), 'a') as f:
            f.write("\n===========================================\n")
            f.write("Epoch {:03d}: Train Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g} \n ".format(
                epoch, *train_loss_history[-1]) + 
                "             Val Loss: {:.5g}, dist RMSD: {:.5g}, angle RMSD: {:.5g}".format(
                *val_loss_history[-1]))
            f.write("\n===========================================\n")

        # save model checkpoints (this overwrites previous checkpoints)
        tf.saved_model.save(model, log_dir)
        checkpoint.save(file_prefix=checkpoint_prefix)

        # plot
        train_loss_arr = np.array(train_loss_history).transpose()
        val_loss_arr = np.array(val_loss_history).transpose()
        plot_train_val(train_loss_arr[0], val_loss_arr[0],
                title="RMSD (dist + angle) loss", savepath=os.path.join(log_dir, "rmsd_all.pdf"))
        plot_train_val(train_loss_arr[1], val_loss_arr[1], 
                title="RMSD (distance matrix) loss", savepath=os.path.join(log_dir, "rmsd_distance_matrix.pdf"))
        plot_train_val(train_loss_arr[2], val_loss_arr[2], 
                title="RMSD (torsion angles) loss", savepath=os.path.join(log_dir, "rmsd_torsion_angles.pdf"))

