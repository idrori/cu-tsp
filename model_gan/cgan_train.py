from cgan_model import *
from cgan_data_loader import *
from cgan_train_utils import *


# load dcalphas and conditioning data (seq+pssm+ss)
training_maps, cond, pdbs = load_training_data(128, 105)
# downscale by factor of 100 as in Anand et al. 2018
training_maps = training_maps/100 

generator = make_generator_model_128()
discriminator = make_discriminator_model_128()
full_gan = make_gan(generator, discriminator)

# training loop variables
EPOCHS = 20
BATCH_SIZE = 15
NOISE_DIM = 100

train(generator, discriminator, full_gan, (training_maps, cond), NOISE_DIM, n_epochs=EPOCHS, n_batch=BATCH_SIZE)