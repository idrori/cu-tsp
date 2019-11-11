import numpy as np
import tensorflow as tf
import pickle

from tensorflow.keras.utils import to_categorical
from tensorflow.keras.preprocessing.text import Tokenizer

from cgan_data_loader import load_test_data

GENERATOR_FILE = 'cgan_generator_128.h5'
TEST_PATH = '../100k/100k/train/testing.pkl'

def nearest_neighbor(dataset, gen_im):
  argmin = np.argmin(np.linalg.norm(dataset-gen_im, axis=(1,2))/np.linalg.norm(dataset, axis=(1,2)))
  print(argmin)
  nearest = dataset[argmin, :, :]
  frobenius = np.linalg.norm(nearest-gen_im)/np.linalg.norm(nearest)
  return frobenius, nearest

def frobenius(real, fake):
  return np.linalg.norm(real-fake)/np.linalg.norm(real)

def test(generator_file):
  # load data
  dcalphas, cond, _ = load_test_data()
  # load generator
  generator = tf.keras.models.load_model(generator_file)
  # generate predictions
  outputs = []
  for ix in range(cond.shape[0]):
    noise = tf.random.normal([1, 100])
    generated = generator.predict([noise, (cond[ix,:,:,:]).reshape(1,128,51,1)], batch_size=1, steps=1)
    generated = np.clip(generated, 0, 1)
    generated = (generated + np.transpose(generated, (0,2,1,3)))/2
    
    frob = frobenius(dcalphas[ix,:,:,0], generated[0,:,:,0])
    print("frobenius norm: {0}".format(frob))
    outputs.append(generated)
  # save predictions as pkl
  if not os.path.exists(os.path.join(os.getcwd(), PREDICTION_DIR)):
    os.makedirs(os.path.join(os.getcwd(), PREDICTION_DIR))
  with open(os.path.join(os.getcwd(),PREDICTION_DIR,'predictions_'+time.strftime("%d%m%y_%H")+'.pkl'), 'wb') as f:
    pickle.dump(outputs, f)

test(GENERATOR_FILE)