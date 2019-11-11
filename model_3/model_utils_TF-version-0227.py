import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.layers import *
from geom_ops import *


# U-Net layers 
def conv_block(x, n_channels, droprate = 0.25):
    """ for UNet """
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = Conv1D(n_channels, 3, padding = 'same', kernel_initializer = 'he_normal')(x) 
    x = Dropout(droprate)(x)
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = Conv1D(n_channels, 3, padding = 'same', kernel_initializer = 'he_normal')(x)
    return x 

def up_block(x, n_channels):
    """ for UNet """
    x = BatchNormalization()(x)
    x = ReLU()(x)
    x = UpSampling1D(size = 2)(x)
    x = Conv1D(n_channels, 2, padding = 'same', kernel_initializer = 'he_normal')(x)
    return x

def Conv_UNet(x, droprate=0.25):
    """ 1-D Convolutional UNet https://arxiv.org/abs/1505.04597 """

    conv0 = Conv1D(192, 3, padding = 'same', kernel_initializer = 'he_normal')(x) 

    conv1 = conv_block(conv0, 128, droprate)
    pool1 = MaxPooling1D(pool_size=2)(conv1)

    conv2 = conv_block(pool1, 192, droprate)
    pool2 = MaxPooling1D(pool_size=2)(conv2)

    conv3 = conv_block(pool2, 384, droprate)
    pool3 = MaxPooling1D(pool_size=2)(conv3)

    conv4 = conv_block(pool3, 512, droprate)

    pool4 = MaxPooling1D(pool_size=2)(conv4)
    conv5 = conv_block(pool4, 1024, droprate)
    up5 = conv5

    up4 = up_block(up5, 512)
    up4 = concatenate([conv4,up4], axis = 2)
    up4 = conv_block(up4, 512, droprate)

    up4 = conv4

    up3 = up_block(up4, 384)
    up3 = concatenate([conv3,up3], axis = 2)
    up3 = conv_block(up3, 384, droprate)

    up2 = up_block(up3, 192)
    up2 = concatenate([conv2,up2], axis = 2)
    up2 = conv_block(up2, 192, droprate)

    up1 = up_block(up2, 128)
    up1 = concatenate([conv1,up1], axis = 2)
    up1 = conv_block(up1, 128, droprate)

    up1 = BatchNormalization()(up1)
    up1 = ReLU()(up1)

    return up1 

# some functions below are adapted from https://github.com/aqlaboratory/rgn

class DistanceMatrix(tf.keras.layers.Layer):
    """ Convert torsion angles to distance matrix 
    using differentiable geometric transformation. """
    def __init__(self):
        super(DistanceMatrix, self).__init__()
    
    def call(self, torsion_angles): 
        coordinates = torsion_angles_to_coordinates(torsion_angles)
        dist = coordinates_to_dist_matrix(coordinates)
        return dist, coordinates

class TorsionAngles(tf.keras.layers.Layer):
    """ computes torsion angles using softmax probabilities 
    and a learned alphabet of angles. (as an alternative to directly predictin angles) """
    def __init__(self, alphabet_size=50):
        super(TorsionAngles, self).__init__()
        self.alphabet = create_alphabet_mixtures(alphabet_size=alphabet_size)
    
    def call(self, probs): 
        torsion_angles = alphabet_mixtures_to_torsion_angles(probs, self.alphabet)
        return torsion_angles

def create_alphabet_mixtures(alphabet_size=50):
    """ Creates alphabet for alphabetized dihedral prediction. """
    init_range = np.pi 
    alphabet_initializer = tf.keras.initializers.RandomUniform(-init_range, init_range)
    alphabet_init = alphabet_initializer(shape=[alphabet_size, NUM_DIHEDRALS], dtype=tf.float32)
    alphabet = tf.Variable(name='alphabet', initial_value=alphabet_init, trainable=True)
    return alphabet  # [alphabet_size, NUM_DIHEDRALS]

def alphabet_mixtures_to_torsion_angles(probs, alphabet):
    """ Converts softmax probabilties + learned mixture components (alphabets) 
        into dihedral angles. 
    """
    torsion_angles = reduce_mean_angle(probs, alphabet)
    return torsion_angles  # [BATCH_SIZE, MAX_LEN, NUM_DIHEDRALS]


def torsion_angles_to_coordinates(torsion_angles, c_alpha_only=True):
    """ Converts dihedrals into full 3D structures. """
    original_shape = torsion_angles.shape
    torsion_angles = tf.transpose(torsion_angles, [1,0,2])
    # converts dihedrals to points ready for reconstruction.

    # torsion_angles: [MAX_LEN=768, BATCH_SIZE=32, NUM_DIHEDRALS=3]
    points = dihedral_to_point(torsion_angles) 
    # points: [MAX_LEN x NUM_DIHEDRALS, BATCH_SIZE, NUM_DIMENSIONS]
             
    # converts points to final 3D coordinates.
    coordinates = point_to_coordinate(points, num_fragments=6, parallel_iterations=4) 
    # [MAX_LEN x NUM_DIHEDRALS, BATCH_SIZE, NUM_DIMENSIONS]
    if c_alpha_only:
        coordinates = coordinates[1::NUM_DIHEDRALS]  # calpha starts from 1
        # [MAX_LEN x 1, BATCH_SIZE, NUM_DIMENSIONS]
    coordinates = tf.transpose(coordinates, [1,0,2])  # do not use reshape
    return coordinates

def coordinates_to_dist_matrix(u, name=None):
    """ Computes the pairwise distance (l2 norm) between all vectors in the tensor.
        Vectors are assumed to be in the third dimension. Op is done element-wise over batch.
    Args:
        u: [MAX_LEN, BATCH_SIZE, NUM_DIMENSIONS]
    Returns:
           [BATCH_SIZE, MAX_LEN, MAX_LEN]
    """
    with tf.name_scope(name, 'pairwise_distance', [u]) as scope:
        u = tf.convert_to_tensor(u, name='u')
        u = tf.transpose(u, [1,0,2])
        
        diffs = u - tf.expand_dims(u, 1)                                 # [MAX_LEN, MAX_LEN, BATCH_SIZE, NUM_DIMENSIONS]
        norms = reduce_l2_norm(diffs, reduction_indices=[3], name=scope) # [MAX_LEN, MAX_LEN, BATCH_SIZE]
        norms = tf.transpose(norms, [2,0,1])
        return norms

def drmsd_dist_matrix(mat1, mat2, batch_seqlen, name=None):
    """
    mat1, mat2: [BATCH_SIZE, MAX_LEN, MAX_LEN]
    batch_seqlen: [BATCH_SIZE,]
    """

    weights = np.zeros(shape=mat1.shape, dtype=np.float32)
    for i, length in enumerate(batch_seqlen):
        weights[i, :length, :length] = 1

    with tf.name_scope(name, 'dRMSD', [mat1, mat2, weights]) as scope:
        mat1 = tf.convert_to_tensor(mat1, name='mat1')
        mat2 = tf.convert_to_tensor(mat2, name='mat2')
        weights = tf.convert_to_tensor(weights, name='weights')
        diffs = mat1 - mat2                      # [BATCH_SIZE, MAX_LEN, MAX_LEN]         
        #diffs = tf.transpose(diffs, [1,2,0])      # [MAX_LEN, MAX_LEN, BATCH_SIZE]
        #weights = tf.transpose(weights, [1,2,0])

        norms = reduce_l2_norm(diffs, reduction_indices=[1, 2], weights=weights, name=scope) # [BATCH_SIZE]
        drmsd = norms / batch_seqlen 
        return drmsd  # [BATCH_SIZE,]

def rmsd_torsion_angle(angles1, angles2, batch_seqlen, name=None):
    """
    angles1, angles2: [BATCH_SIZE, MAX_LEN]
    batch_seqlen: [BATCH_SIZE,]
    """

    weights = np.zeros(shape=angles1.shape, dtype=np.float32)
    for i, length in enumerate(batch_seqlen):
        weights[i, :length] = 1.0

    with tf.name_scope(name, 'RMSD_torsion', [angles1, angles2, weights]) as scope:
        angles1 = tf.convert_to_tensor(angles1, name='angles1')
        angles2 = tf.convert_to_tensor(angles2, name='angles2')
        weights = tf.convert_to_tensor(weights, name='weights')
        diffs = angles1 - angles2                      # [BATCH_SIZE, MAX_LEN]         

        norms = reduce_l2_norm(diffs, reduction_indices=[1], weights=weights, name=scope) # [BATCH_SIZE]
        drmsd = norms / tf.sqrt(batch_seqlen)
        return drmsd  # [BATCH_SIZE,]

def seq2ngrams(seqs, n = 1):
    return np.array([[seq[i : i + n] for i in range(len(seq))] for seq in seqs])

def plot_train_val(train, val, title=None, savepath=None):
    fig, ax = plt.subplots(1, 1, figsize=(6,4))
    ax.plot(train, c='g', label='train')
    ax.plot(val, c='b', label='val')
    ax.legend()
    if not title is None:
        ax.set_title(title)
    fig.savefig(savepath)
    plt.close()

def plot_dist_matrix(pred, gt, protein_names, lengths, scores, savepath=None):
    assert(pred.shape[0] == gt.shape[0])
    fig, axes = plt.subplots(2, pred.shape[0], figsize=(4 * pred.shape[0],8))
    for i, pname in enumerate(protein_names):
        axes[0, i].imshow(pred[i, :lengths[i], :lengths[i]])
        axes[0, i].set_title(pname + " prediction ({:.4g})".format(scores[i]))
    for i, pname in enumerate(protein_names):
        axes[1, i].imshow(gt[i, :lengths[i], :lengths[i]])
        axes[1, i].set_title(pname + " ground truth")
    plt.savefig(savepath)
    plt.close()

def rmsd_kaggle(rmsd_batch, seqlen_batch):
    """ rmsd across the entire batch """
    norm = tf.reduce_sum(tf.multiply(tf.square(rmsd_batch), seqlen_batch))
    rmsd = tf.sqrt(norm / tf.reduce_sum(seqlen_batch))
    return rmsd 

