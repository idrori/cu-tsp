import math
import numpy as np
import tensorflow as tf

from proteintf.dataset_readers.datasets.utils import get_position_encoding, get_output_info


class SampledSingleSeqDataset:
    def __init__(self, data, batch_size, inputs, targets, target_stats,
                 position_emb_dim, window_size, shuffle, msas_dim, test):
        assert not test, 'Test functionalities have not been implemented for SampledSingleSeqDataset'

        output_types, output_shapes, padded_shapes = get_output_info(inputs,
                                                                     targets,
                                                                     msas_dim,
                                                                     position_emb_dim,
                                                                     window_size,
                                                                     test=test)
        # Each instance is sampled using a fixed window size
        # Record the original index, and the start position of each instance
        self.instance_id = []
        self.window_size = window_size

        def _generator():
            for i in range(len(data)):
                length = data['length_aas'][i]
                mask = np.ones(length, dtype=np.int64)
                if position_emb_dim is not None:
                    position_encoding = get_position_encoding(length, position_emb_dim)
                else:
                    position_encoding = None
                for start in range(max(length - window_size + 1, 1)):
                    self.instance_id += [(i, start)]
                    end = start + window_size

                    output = [mask[start: end]]
                    if position_encoding is not None:
                        output += [position_encoding[start: end]]
                    output += [data[col][i][start: end] for col in inputs]
                    if not test:
                        tgt = [data[col][i][start: end] for col in targets]
                        output = (tuple(output), tuple(tgt))
                    else:
                        output = tuple(output)
                    yield output

        dataset = tf.data.Dataset.from_generator(_generator, output_types, output_shapes)
        if shuffle:
            dataset = dataset.shuffle(buffer_size=5000)

        self.dataset = dataset.repeat().padded_batch(batch_size, padded_shapes)
        self._len = math.ceil(sum([max(l - window_size + 1, 1) for l in data['length_aas']]) / batch_size)

        self.length_aas = data['length_aas']
        self.target_stats = target_stats

    def __len__(self):
        return self._len

    def format_predictions(self, predictions):
        raise NotImplementedError('format_predictions has not been implemented for SampledSingleSeqDataset')
