import os
import time
import pickle
import numpy as np

from numpy import zeros
from numpy import ones
from numpy.random import randn
from numpy.random import randint

OUTPUT_DIR = 'cgan_output'

def normalize_predictions(generated_maps):
  # Clamp generated values to > 0
  generated_maps = np.clip(generated_maps, 0, 1)
  # Enforce symmetry for generated_maps
  generated_maps = (generated_maps + np.transpose(generated_maps, (0,2,1,3)))/2
  return generated_maps

def grab_real_data(dataset, n_samples):
  maps, conds = dataset
  # choose random instances
  ix = randint(0, maps.shape[0], n_samples)
  # select dcalpha maps and conditional data
  map_samples, cond_samples = maps[ix], conds[ix]
  return map_samples, cond_samples

# select real samples
def generate_real_samples(dataset, n_samples):
  # select images and labels
  x, conds = grab_real_data(dataset, n_samples)
  # generate soft T/F labels for reals
  y = np.random.uniform(low=0.9, high=1, size=(n_samples,1))
  # add a few mislabels to help train
  y[:1] = 1-y[:1]
  np.random.shuffle(y)
  return [x, conds], y
 
# generate points in latent space as input for the generator
def generate_latent_points(latent_dim, n_samples, dataset):
  maps, cond = dataset
  # generate points in the latent space
  x_input = randn(latent_dim * n_samples)
  # reshape into a batch of inputs for the network
  z_input = x_input.reshape(n_samples, latent_dim)
  # grab random conditioning data to condition on
  x_real, labels = grab_real_data(dataset, n_samples)
  labels = labels.reshape(labels.shape[0], 128, 51, 1)
  return [z_input, labels, x_real]
 
# use the generator to generate n fake examples
def generate_fake_samples(generator, latent_dim, n_samples, dataset):
  # generate points in latent space along with conditioning data and corresponding real maps
  z_input, labels_input, real_maps = generate_latent_points(latent_dim, n_samples, dataset)
  # predict outputs
  generated_maps = generator.predict([z_input, labels_input])
  # clamp and enforce symmetry
  generated_maps = normalize_predictions(generated_maps)
  # create soft T/F labels for fakes
  y = np.random.uniform(low=0, high=0.1, size=(n_samples,1))
  return [generated_maps, labels_input, real_maps], y

def frobenius(real, fake):
  return np.linalg.norm(real-fake)/np.linalg.norm(real)

def nearest_neighbor(dataset, gen_im):
  argmin = np.argmin(np.linalg.norm(dataset-gen_im, axis=(1,2))/np.linalg.norm(dataset, axis=(1,2)))
  print(argmin)
  nearest = dataset[argmin, :, :]
  frobenius = np.linalg.norm(nearest-gen_im)/np.linalg.norm(nearest)
  return frobenius, nearest

def plot_training(real, fake):
    display.clear_output(wait=True)
    fig = plt.figure(figsize=(8,8))
    plt.subplot(1, 2, 1)
    plt.imshow(real[0,:,:,0],cmap='winter')
    plt.subplot(1, 2, 2)
    plt.imshow(fake[0,:,:,0],cmap='winter')
    plt.show()

def train(g_model, d_model, gan_model, dataset, latent_dim, n_epochs=20, n_batch=15, plot=False):
  # track variables
  # TODO helper data structure for this
  g_losses = []
  d_losses1 = []
  d_losses2 = []
  g_frobs = []

  bat_per_epo = int(dataset[0].shape[0] / n_batch)
  half_batch = int(n_batch / 2)
  # manually enumerate epochs
  for i in range(n_epochs):
    # enumerate batches over the training set
    for j in range(bat_per_epo):
      # get randomly selected 'real' samples
      [X_real, labels_real], y_real = generate_real_samples(dataset, half_batch)
      # update discriminator model weights
      d_loss1, _ = d_model.train_on_batch([X_real, labels_real], y_real)
      # generate 'fake' examples aalongside their corresponding real maps
      [X_fake, labels, X_corresponding], y_fake = generate_fake_samples(g_model, latent_dim, half_batch, dataset)
      # update discriminator model weights
      d_loss2, _ = d_model.train_on_batch([X_fake, labels], y_fake)
      
      if (plot and j==0):
        plot_training(X_real, X_fake)
      
      # prepare points in latent space as input for the generator
      [z_input, labels_input, _] = generate_latent_points(latent_dim, n_batch, dataset)
      # create inverted labels for the fake samples
      y_gan = ones((n_batch, 1))
      # update the generator via the discriminator's error
      g_loss = gan_model.train_on_batch([z_input, labels_input], y_gan)
			# summarize loss on this batch
      print('>%d, %d/%d, d1=%.5f, d2=%.5f g=%.5f' %
				(i+1, j+1, bat_per_epo, d_loss1, d_loss2, g_loss))
      
      # track variables 
      # TODO helper data structure for this
      g_losses.append(g_loss)
      d_losses1.append(d_loss1)
      d_losses2.append(d_loss2)
      g_frobs.append(frobenius(X_corresponding, X_fake))
    
    # save the generator every 5 epochs
    if (i % 5 == 0):
      g_model.save('cgan_generator_128_epoch_{0}.h5'.format(i))

  # save the losses and frobenius norms
  print(os.getcwd())
  if not os.path.exists(os.path.join(os.getcwd(), OUTPUT_DIR)):
    os.makedirs(os.path.join(os.getcwd(), OUTPUT_DIR))
  with open(os.path.join(os.getcwd(),OUTPUT_DIR,'training_stats_'+time.strftime("%d%m%y_%H")+'.pkl'), 'wb') as f:
    pickle.dump((d_losses1,d_losses2,g_losses,g_frobs), f)
  # save the generator model
  g_model.save('cgan_generator_128.h5')
  d_model.save('cgan_discriminator_128.h5')