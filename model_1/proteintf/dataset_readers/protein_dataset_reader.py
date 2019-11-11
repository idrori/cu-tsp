import os
import pickle
import numpy as np
import pandas as pd

from proteintf.dataset_readers.datasets import SingleSeqDataset, SampledSingleSeqDataset

_AAS_VOCAB = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
_Q8S_VOCAB = ('-', 'B', 'E', 'G', 'H', 'I', 'S', 'T')
_TRAIN_COLUMNS = ('index', 'pdb', 'length_aas', 'pdb_aas', 'q8s', 'dcalphas', 'psis', 'phis', 'msas')
_TEST_COLUMNS = ('index', 'pdb', 'length_aas', 'pdb_aas', 'q8s', 'msas')
_MSAS_DIM = 21


def ensure_args(inputs, targets, window_size):
    assert len(inputs) == len(set(inputs)), 'No duplicated inputs'
    assert len(targets) == len(set(targets)), 'No duplicated targets'
    assert len(set(inputs)) > 0, 'At least one input'
    assert len(set(targets)) > 0, 'At least one target'
    assert len(set(inputs) - {'pdb_aas', 'q8s', 'msas'}) == 0, \
        'Inputs must be a subset of [pdb_aas, q8s, msas]'
    assert len(set(targets) - {'dcalphas', 'psis', 'phis'}) == 0, \
        'Targets must be a subset of [dcalphas, psis, phis]'
    assert window_size is None or \
        set(targets) == {'dcalphas'} or 'dcalphas' not in targets, \
        'When sampling data, cannot predict dcalphas and (psis or phis) at the same time'


def preprocessor(tb, aa_to_idx, q8_to_idx):
    tb['pdb_aas'] = [np.array([aa_to_idx[aa] for aa in aas], dtype=np.int64) for aas in tb['pdb_aas']]
    tb['q8s'] = [np.array([q8_to_idx[q8] for q8 in q8s], dtype=np.int64) for q8s in tb['q8s']]
    tb['msas'] = [np.array(msas, dtype=np.float32).transpose() for msas in tb['msas']]

    if 'dcalphas' in tb:
        for col in ['dcalphas', 'phis', 'psis']:
            tb[col] = [np.array(dcalphas, dtype=np.float32) for dcalphas in tb[col]]
    return tb


class ProteinDatasetReader:
    def __init__(self,
                 root='',
                 train_file_prefix='train_fold_',
                 total_folds=10,
                 aas_vocab=_AAS_VOCAB,
                 q8s_vocab=_Q8S_VOCAB,
                 train_columns=_TRAIN_COLUMNS,
                 test_columns=_TEST_COLUMNS):
        aa_to_idx = {aa: idx + 1 for idx, aa in enumerate(aas_vocab)}  # 0 reserved for padding
        q8_to_idx = {q8: idx + 1 for idx, q8 in enumerate(q8s_vocab)}  # 0 reserved for padding

        self.train_data = {}
        for i in range(1, total_folds + 1):
            with open(os.path.join(root, '{}{}.pkl'.format(train_file_prefix, i)), 'rb') as fin:
                tb = pd.DataFrame(list(zip(*pickle.load(fin))), columns=train_columns)
                self.train_data[i] = preprocessor(tb, aa_to_idx, q8_to_idx)
        with open(os.path.join(root, 'test.pkl'), 'rb') as fin:
            tb = pd.DataFrame(list(zip(*pickle.load(fin))), columns=test_columns)
            self.test_data = preprocessor(tb, aa_to_idx, q8_to_idx)

        self.inputs = None
        self.targets = None
        self.position_emb_dim = None
        self.window_size = None
        self.msas_dim = None
        self.target_stats = None

    def set_params(self,
                   inputs=('pdb_aas',),
                   targets=('dcalphas', 'psis', 'phis',),
                   position_emb_dim=None,  # None if not using position encodings
                   window_size=None,  # None if not sampling data
                   msas_dim=_MSAS_DIM):
        ensure_args(inputs, targets, window_size)
        self.inputs = inputs
        self.targets = targets
        self.position_emb_dim = position_emb_dim
        self.window_size = window_size
        self.msas_dim = msas_dim

    def make_train_dataset(self,
                           folds,
                           batch_size,
                           max_length=None,
                           normalize_target=False,
                           shuffle=True):
        assert self.inputs is not None, 'Call set_params before making datasets'
        assert max_length is None or self.window_size is None, \
            'Cannot set max_length when sampling data'

        train_data = pd.concat([self.train_data[i] for i in folds], ignore_index=True)
        if normalize_target:
            self.target_stats = {}
            for t in self.targets:
                target_data = np.concatenate([item.reshape(-1) for item in train_data[t]])
                mean, std = target_data.mean(), target_data.std()
                train_data[t] = (train_data[t] - mean) / std
                self.target_stats[t] = {'mean': mean, 'std': std}
        return self._make_dataset(train_data, batch_size, max_length, shuffle)

    def make_valid_dataset(self,
                           folds,
                           batch_size,
                           max_length=None,
                           normalize_target=False,
                           shuffle=True):
        assert self.inputs is not None, 'Call set_params before making datasets'
        assert max_length is None or self.window_size is None, \
            'Cannot set max_length when sampling data'
        assert not normalize_target or self.target_stats is not None, \
            'Normalzing targets requires target_stats, which is generated by ' \
            'calling make_train_dataset with normalize_target=True'

        valid_data = pd.concat([self.train_data[i] for i in folds], ignore_index=True)
        if normalize_target:
            for t in self.targets:
                mean, std = self.target_stats[t]['mean'], self.target_stats[t]['std']
                valid_data[t] = (valid_data[t] - mean) / std
                self.target_stats[t] = {'mean': mean, 'std': std}
        return self._make_dataset(valid_data, batch_size, max_length, shuffle)

    def make_test_dataset(self, batch_size):
        return self._make_dataset(self.test_data, batch_size, test=True)

    def _make_dataset(self, data, batch_size, max_length=None, shuffle=True, test=False):
        if self.window_size is None:
            return SingleSeqDataset(data=data,
                                    batch_size=batch_size,
                                    inputs=self.inputs,
                                    targets=self.targets,
                                    target_stats=self.target_stats,
                                    position_emb_dim=self.position_emb_dim,
                                    max_length=max_length,
                                    shuffle=shuffle,
                                    msas_dim=self.msas_dim,
                                    test=test)
        elif 'dcalphas' not in self.targets:
            return SampledSingleSeqDataset(data=data,
                                           batch_size=batch_size,
                                           inputs=self.inputs,
                                           targets=self.targets,
                                           target_stats=self.target_stats,
                                           position_emb_dim=self.position_emb_dim,
                                           window_size=self.window_size,
                                           shuffle=shuffle,
                                           msas_dim=self.msas_dim,
                                           test=test)
        else:
            raise NotImplementedError('SampledTwoSeqDataset has not been implemented')

    def get_target_stats(self):
        assert self.target_stats is not None, \
            'Target stats do not exist; Generate by ' \
            'calling make_train_dataset with normalize_target=True'
        return self.target_stats
