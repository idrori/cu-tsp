import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras import backend as k
from tensorflow.keras.models import Model


def make_generator_model_128():
  """Makes the generator network"""
  
  in_seq = layers.Input(shape=(128, 51, 1)) # one hot encoded
  seq = layers.Flatten()(in_seq)
  seq = layers.Dense(4*4*512, use_bias=False)(seq) # same size as noise
  seq = layers.Reshape((4, 4, 512))(seq)
  
  # input noise
  in_noise = layers.Input(shape=(100,))
  noise = layers.Dense(4*4*512, use_bias=False)(in_noise)
  noise = layers.Reshape((4, 4, 512))(noise)
  merge = layers.Concatenate()([noise, seq])
  
  gen = layers.Conv2DTranspose(512, (4, 4), strides=(1, 1), padding='same', use_bias=False)(merge)
  gen = layers.BatchNormalization()(gen)
  gen = layers.LeakyReLU()(gen)
  
  gen = layers.Conv2DTranspose(256, (4, 4), strides=(2, 2), padding='same', use_bias=False)(gen)
  gen = layers.BatchNormalization()(gen)
  gen = layers.LeakyReLU()(gen)
  
  gen = layers.Conv2DTranspose(128, (4, 4), strides=(4, 4), padding='same', use_bias=False)(gen)
  gen = layers.BatchNormalization()(gen)
  gen = layers.LeakyReLU()(gen)
  
  gen = layers.Conv2DTranspose(64, (4, 4), strides=(2, 2), padding='same', use_bias=False)(gen)
  gen = layers.BatchNormalization()(gen)
  gen = layers.LeakyReLU()(gen)
  
  output = layers.Conv2DTranspose(1, (4, 4), strides=(2, 2), padding='same', use_bias=False, activation='tanh')(gen)
  model = Model([in_noise, in_seq], output)

  # add optimizer and compile
  opt = tf.keras.optimizers.Adam(lr=1e-4, beta_1=0.65, beta_2=0.999)
  model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
  
  return model
  
def make_discriminator_model_128():
  """Makes the discriminator network"""

  in_seq = layers.Input(shape=(128, 51, 1)) # one hot encoded
  seq = layers.Flatten()(in_seq)
  
  seq = layers.Dense(128*128)(seq)
  seq = layers.Reshape((128, 128, 1))(seq)  # same shape as input matrix
  
  in_matrix = layers.Input(shape=(128,128,1))
  
  # concat the sequence as a channel
  merge = layers.Concatenate()([in_matrix, seq])
  
  d = layers.Conv2D(64, (4,4), strides=(4, 4), padding='same')(merge)
  d = layers.LeakyReLU(alpha=0.2)(d)
  d = layers.Dropout(0.1)(d)
  
  d = layers.Conv2D(128, (4,4), strides=(2, 2), padding='same')(d)
  d = layers.BatchNormalization()(d)
  d = layers.LeakyReLU(alpha=0.2)(d)
  d = layers.Dropout(0.1)(d)
  
  d = layers.Conv2D(256, (4,4), strides=(4, 4), padding='same')(d)
  d = layers.BatchNormalization()(d)
  d = layers.LeakyReLU(alpha=0.2)(d)
  d = layers.Dropout(0.1)(d)
  
  d = layers.Conv2D(512, (4,4), strides=(2, 2), padding='same')(d)
  d = layers.BatchNormalization()(d)
  d = layers.LeakyReLU(alpha=0.2)(d)
  d = layers.Dropout(0.1)(d)
  
  d = layers.Conv2D(1, (4, 4), strides=(1, 1), padding='same')(d)
  flat = layers.Flatten()(d)
  output = layers.Dense(1, activation='sigmoid')(flat)
  model = Model(inputs=[in_matrix, in_seq], outputs=[output])
  
  # add optimizer and compile
  opt = tf.keras.optimizers.Adam(lr=1e-4, beta_1=0.5, beta_2=0.999)
  model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
  
  return model

# define the combined generator and discriminator model, for updating the generator
def make_gan(g_model, d_model):
  # make weights in the discriminator not trainable
  d_model.trainable = False
  # get noise and label inputs from generator model
  gen_noise, gen_label = g_model.input
  # get image output from the generator model
  gen_output = g_model.output

  # Clamp generated values to > 0
  gen_output = tf.clip_by_value(gen_output,0,1)

  # Enforce symmetry for generated_maps
  gen_output = tf.math.divide(tf.math.add(gen_output, tf.transpose(gen_output, perm=[0, 2, 1, 3] )), 2)
  
  # connect image output and label input from generator as inputs to discriminator
  gan_output = d_model([gen_output, gen_label])
  # define gan model as taking noise and label and outputting a classification
  model = Model([gen_noise, gen_label], gan_output)
  # compile model
  opt = tf.keras.optimizers.Adam(learning_rate=1e-4, beta_1=0.5, beta_2=0.999)
  model.compile(loss='binary_crossentropy', optimizer=opt)
  return model