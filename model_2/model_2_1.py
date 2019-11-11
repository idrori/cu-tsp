## Load and prepare data

# set model name, which determines the log directory (to store checkpoints and trained model instances)
model_name = 'model-2-1'
data_dir = '/home/yg2541/train/'
log_dir = '/home/yg2541/logs/{}'.format(model_name)

# split data into subsets
list_of_folds_use = [list(range(1, 27)), list(range(27,53)), list(range(53,80)), list(range(80, 106))]
starting_epochs = [1, 11, 21,31]

# loading modules and util.py

print('Importing modules...')
import os, sys, time 
import pickle
import numpy as np
import tensorflow as tf
print("TensorFlow version:", tf.__version__)
from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.layers import *
from datetime import datetime
sys.path.append("..")
from utils.rgn_utils import *
print('Done importing...')

def load_data(data_dir, folds_use = [1,2,3]):
    pdbs = []; lengths = []; aa = []; ss = []; dcalphas = []; coords = []; psis = []; phis = []; pssms = []
    for fold in folds_use:
        # load current pickle file
        dict_curr_file = pickle.load(open(os.path.join(data_dir, 'training_100_{}.pkl'.format(fold)), 'rb'))
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

for super_epoch_idx in range(3):
    folds_use = list_of_folds_use[super_epoch_idx]
    starting_epoch = starting_epochs[super_epoch_idx] 
    # folds_use = list(range(1, 5)) 
    data_fields = load_data(data_dir, folds_use) # data_fields = [pdbs, lengths, aa, ss, dcalphas, coords, psis, phis, pssms]
    print('Total number of proteins loaded from data: {}'.format(len(data_fields[0])))

    ################
    # preprocessing
    ################

    # remove proteins of extreme lengths
    maxlen_seq, minlen_seq = 320, 30
    #maxlen_seq, minlen_seq = 384, 30
    length_mask = (data_fields[1] < maxlen_seq) & (data_fields[1] > minlen_seq)
    for j in range(len(data_fields)): # in-place modification of list data_fields
        data_fields[j] = data_fields[j][length_mask]
    # unroll data fields
    pdbs, lengths, aa, ss, dcalphas, coords, psis, phis, pssms = data_fields
    print("Number of protein sequences between max and min lengths: {}".format(len(pdbs)))

    # reshape some matrix-valued data
    pssms = [np.stack(x).transpose().astype(np.float32) for x in pssms]
    phis = [np.array(x) for x in phis]
    psis = [np.array(x) for x in psis]

    # remove None
    for idx in range(len(phis)):
        phis[idx] = np.array([x if x!=None else 0 for x in phis[idx]])
        psis[idx] = np.array([x if x!=None else 0 for x in psis[idx]])
                
    # pad target distance matrix to the same size
    dcalphas_pad = np.zeros((len(dcalphas), maxlen_seq, maxlen_seq), dtype=np.float32)
    for i in range(len(dcalphas)):
        length = dcalphas[i].shape[0]
        dcalphas_pad[i, :length, :length] = dcalphas[i]

    # pad target torsion angles to the same size with 0
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

    # extra PSSM features
    pssm_dim = pssms[0].shape[1]
    pssms_padded = np.zeros([len(pssms), maxlen_seq, pssm_dim], dtype=np.float32)
    for i in range(len(pssms)):
        pssms_padded[i, :pssms[i].shape[0], :] = pssms[i]

    # make dir for saving training logs...create tokenizer_encoder pickles when it's the first time
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)
    with open(os.path.join(log_dir, 'tokenizer_encoder.pickle'), 'wb') as handle:
        pickle.dump(tokenizer_encoder, handle)

    with open(os.path.join(log_dir, 'tokenizer_encoder_q8s.pickle'), 'wb') as handle:
        pickle.dump(tokenizer_encoder_q8s, handle)

    # train/val split
    permu = np.random.permutation(len(pdbs))

    validation_ratio = 0.05

    cut_off = int(len(pdbs)*(1-validation_ratio))
    training_idx = permu[:cut_off]
    validation_idx = permu[cut_off:]

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

    train_target_psis[:,0] = 0
    y_train_phis = train_target_phis[training_idx]
    y_val_phis = train_target_phis[validation_idx]
    y_train_psis = train_target_psis[training_idx]
    y_val_psis = train_target_psis[validation_idx]

    pdb_names_train = pdbs[training_idx]
    pdb_names_val = pdbs[validation_idx]

    print('train and validate sizes: {}, {}'.format(len(y_train_dist_matrix), len(y_val_dist_matrix)))

    ## Build and train the model
    ############################################
    # build a neural network model
    #############################################

    # Here we use layer functions to build neural network graph (as opposed to defining a Model class)

    # Input layer: (None, ) shape allows the network to be flexible to input sequence length 
    # inputs: [protein_sequence, q8_sequence, MSA features]
    inputs = [Input(shape = (None, )), Input(shape = (None, )), Input(shape = (None, 21))]

    # Word (amino acid) embedding layers 
    embedding1 = Embedding(input_dim = n_words, output_dim = 128, input_length = None, input_shape=(None,))
    embedding2 = Embedding(input_dim = n_words_q8s, output_dim = 64, input_length = None, input_shape=(None,))
    embed_input = embedding1(inputs[0]) # raw seq
    embed_q8s = embedding2(inputs[1]) # q8

    # Concatenate different types of inputs 
    merged_input = concatenate([embed_input, embed_q8s, inputs[2]], axis = 2)

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

    ##################################################
    # Train the model
    ##################################################

    batch_size = 32
    angle_scale = 180 / np.pi  # convert angles from [-pi, pi] to [-180, 180]
    n_iter = int(X_train.shape[0] / batch_size)

    n_epochs = 30
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

    # starting_epoch = 1
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

            # Compute loss 
            with tf.GradientTape() as tape:
                torsion_angles, coordinates = model([batch_seq, batch_q8s, batch_pssm])
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

            with tf.GradientTape() as tape:
                torsion_angles, coordinates = model([batch_seq, batch_q8s, batch_pssm])
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
    #         val_rmsd_dist_norm(loss_drmsd_normalized)
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
