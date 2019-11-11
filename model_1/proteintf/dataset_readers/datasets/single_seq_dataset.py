import math
import numpy as np
import tensorflow as tf

from proteintf.dataset_readers.datasets.utils import get_position_encoding, get_output_info


class SingleSeqDataset:
    def __init__(self, data, batch_size, inputs, targets, target_stats,
                 position_emb_dim, max_length, shuffle, msas_dim, test):
        output_types, output_shapes, padded_shapes = get_output_info(inputs,
                                                                     targets,
                                                                     msas_dim,
                                                                     position_emb_dim,
                                                                     test=test)
        if test:
            max_length = None
            shuffle = False
            batch_size = len(data)

        def _generator():
            for i in range(len(data)):
                length = data['length_aas'][i]
                if max_length is not None and length > max_length:
                    length = max_length

                # The first input is the mask
                output = [np.ones(length, dtype=np.int64)]
                if position_emb_dim is not None:
                    output += [get_position_encoding(length, position_emb_dim)]

                output += [data[col][i][:length] for col in inputs]
                if not test:
                    tgt = [data[col][i][:length, :length] if col == 'dcalphas'
                            else data[col][i][:length] for col in targets]
                    output = (tuple(output), tuple(tgt))
                else:
                    output = tuple(output)

                yield output

        dataset = tf.data.Dataset.from_generator(_generator, output_types, output_shapes)
        if shuffle:
            dataset = dataset.shuffle(buffer_size=min(len(data), 10000))

        self.dataset = dataset.repeat().padded_batch(batch_size, padded_shapes)
        self._len = math.ceil(len(data) / batch_size)

        self.length_aas = data['length_aas']
        self.targets = targets
        self.target_stats = target_stats

    def __len__(self):
        return self._len

    def format_predictions(self, predictions):
        if isinstance(predictions, np.ndarray):
            predictions = [predictions]
        outputs = {}
        for i, tgt in enumerate(self.targets):
            pred = predictions[i]
            if self.target_stats is not None:
                mean, std = self.target_stats[tgt]['mean'], self.target_stats[tgt]['std']
                pred = pred * std + mean
            outputs[tgt] = []
            for l, p in zip(self.length_aas, pred):
                if tgt == 'dcalphas':
                    p = p[:l, :l]
                    p = np.triu(p, 1) + np.tril(p.transpose(), -1)
                elif tgt == 'phis':
                    p = p[:l]
                    p[-1] = 0.
                else:
                    p = p[:l]
                    p[0] = 0.
                outputs[tgt] += [p]
        return outputs
