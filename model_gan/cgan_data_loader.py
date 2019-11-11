import numpy as np
import pickle

from tensorflow.keras.utils import to_categorical
from tensorflow.keras.preprocessing.text import Tokenizer

def load_training_data(frag_size, n_folds):
  """
  Loads tuple of dcalpha matrices, ss arrays, and pssms into numpy nd arrays
  Output dcalpha ndarray of shape (num_examples, frag_size, frag_size, 1)
  """
  
  dcalphas = np.array([])
  seqs = np.array([]) # nx1
  pssms = np.array([]) # nx21
  secondary = np.array([]) #nx1
  m_frags = []
  s_frags = []
  p_frags = []
  ss_frags = []
  all_pdbs = []
  for fold in range(1, n_folds): # Full Dataset is 75 folds
    path = 'train/training_100_{}.pkl'.format(fold)
    data = pickle.load(open(path, 'rb'))
    pdbs = list(data.keys())
    all_pdbs += pdbs

    for p in pdbs:
      full_matrix = np.array(data[p]['dcalpha'])
      full_seq = np.array(data[p]['aa']) # amino acid sequence
      full_pssm = np.array(data[p]['pssm']) # pssm
      full_ss = np.array(data[p]['ss']) # secondary structure

      num_bonds = full_matrix.shape[0]
      
      # get non-overlapping fragments of frag_size length
      if (full_matrix.shape[0] >= frag_size): 
        idx = 0
        for i in range(frag_size, num_bonds, frag_size):
          matrix_frag = full_matrix[idx:i, idx:i]
          seq_frag = full_seq[()][idx:i]
          pssm_frag = full_pssm[:, idx:i]
          ss_frag = full_ss[()][idx:i]
          
          idx = i # update start index
          m_frags.append(matrix_frag)
          s_frags.append(seq_frag)
          p_frags.append(pssm_frag)
          ss_frags.append(ss_frag)


      num_loaded = len(m_frags)
    print("fold {} complete, so far loaded {} pairwise matrices of size {}"
          .format(fold, num_loaded, frag_size))

  dcalphas = np.stack(m_frags)
  
  # one hot encode the sequences
  tokenizer = Tokenizer(char_level=True)
  tokenizer.fit_on_texts(s_frags)
  int_seq = tokenizer.texts_to_sequences(s_frags)
  s_frags = to_categorical(int_seq)
  
  # one hot encode the sequences
  tokenizer2 = Tokenizer(char_level=True)
  tokenizer2.fit_on_texts(ss_frags)
  int_seq = tokenizer2.texts_to_sequences(ss_frags)
  ss_frags = to_categorical(int_seq)
  
  seqs = np.stack(s_frags)
  pssms = np.stack(p_frags)
  secondary = np.stack(ss_frags)

  dcalphas = dcalphas.reshape(dcalphas.shape[0], frag_size, frag_size, 1).astype('float32')
  pssms = pssms.reshape(pssms.shape[0], 128, 21).astype('float32')
  
  cond = np.array([])
  cond = np.concatenate((seqs,pssms,secondary), axis=2) # thing we are conditioning on nx51
  cond = cond.reshape(cond.shape[0], 128, 51, 1).astype('float32')
  print(dcalphas.shape)
  print(cond.shape)
  print(len(all_pdbs))

  return dcalphas, cond, all_pdbs

def load_test_data(test_path):
  """
  Loads tuple of dcalpha matrices, ss arrays, and pssms into numpy nd arrays
  Output dcalpha ndarray of shape (num_examples, frag_size, frag_size, 1)
  """
  frag_size = 128 
  dcalphas = np.array([])
  seqs = np.array([]) # nx1
  pssms = np.array([]) # nx21
  secondary = np.array([]) #nx1
  m_frags = []
  s_frags = []
  p_frags = []
  ss_frags = []
  all_pdbs = []
  
  data = pickle.load(open(test_path, 'rb'))
  pdbs = list(data.keys())
  all_pdbs += pdbs

  for p in pdbs:
    full_matrix = np.array(data[p]['dcalpha'])
    full_seq = np.array(data[p]['aa']) # amino acid sequence
    full_pssm = np.array(data[p]['pssm']) # pssm
    full_ss = np.array(data[p]['ss']) # secondary structure

    num_bonds = full_matrix.shape[0]
    
    # get non-overlapping fragments of frag_size length
    if (full_matrix.shape[0] >= frag_size): 
      idx = 0
      for i in range(frag_size, num_bonds, frag_size):
        matrix_frag = full_matrix[idx:i, idx:i]
        seq_frag = full_seq[()][idx:i]
        pssm_frag = full_pssm[:, idx:i]
        ss_frag = full_ss[()][idx:i]
        
        idx = i # update start index
        m_frags.append(matrix_frag)
        s_frags.append(seq_frag)
        p_frags.append(pssm_frag)
        ss_frags.append(ss_frag)

  dcalphas = np.stack(m_frags)
  
  # one hot encode the sequences
  tokenizer = Tokenizer(char_level=True)
  tokenizer.fit_on_texts(s_frags)
  int_seq = tokenizer.texts_to_sequences(s_frags)
  s_frags = to_categorical(int_seq)
  
  # one hot encode the sequences
  tokenizer2 = Tokenizer(char_level=True)
  tokenizer2.fit_on_texts(ss_frags)
  int_seq = tokenizer2.texts_to_sequences(ss_frags)
  ss_frags = to_categorical(int_seq)
  
  seqs = np.stack(s_frags)
  pssms = np.stack(p_frags)
  secondary = np.stack(ss_frags)

  dcalphas = dcalphas.reshape(dcalphas.shape[0], frag_size, frag_size, 1).astype('float32')
  pssms = pssms.reshape(pssms.shape[0], 128, 21).astype('float32')
  
  cond = np.array([])
  cond = np.concatenate((seqs,pssms,secondary), axis=2) # thing we are conditioning on nx51
  cond = cond.reshape(cond.shape[0], 128, 51, 1).astype('float32')
  print(dcalphas.shape)
  print(cond.shape)
  print(len(all_pdbs))
  