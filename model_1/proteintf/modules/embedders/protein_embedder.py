import tensorflow as tf
from tensorflow.keras import layers

from proteintf.modules.embedders.embedder import Embedder


class ProteinEmbedder:
    @classmethod
    def from_params(cls, params):
        input_names = params['input_names']
        embedders = params['embedders']
        merge_mode = params.get('merge_mode', 'concat')
        project_dim = params.get('project_dim', None)
        dropout = params.get('dropout', 0.)

        embedder_layers = {}
        output_dim = 0
        for n in input_names:
            embedder_type = embedders[n]['type']
            embedder, emb_dim = Embedder.by_name(embedder_type).from_params(embedders[n])
            embedder_layers[n] = embedder
            if merge_mode == 'concat':
                output_dim += emb_dim
            elif merge_mode == 'sum' or merge_mode == 'ave':
                output_dim = emb_dim
            else:
                raise ValueError('Unrecognized merge mode {}'.format(merge_mode))
        output_dim = emb_dim if project_dim is None else project_dim

        def _embedder(inputs):
            assert len(inputs) == len(input_names), \
                f'Expected {len(input_names)} input(s) {input_names}, ' \
                    f'but got {len(inputs)}'

            x = []
            for n, inp in zip(input_names, inputs):
                for layer in embedder_layers[n]:
                    inp = layer(inp)
                x += [inp]

            if len(x) > 1:
                if merge_mode == 'concat':
                    x = layers.concatenate(x)
                elif merge_mode == 'sum':
                    x = layers.add(x)
                elif merge_mode == 'ave':
                    x = layers.average(x)
                else:
                    raise ValueError('Unrecognized merge mode {}'.format(merge_mode))
            else:
                x = x[0]

            if project_dim is not None:
                x = layers.Dense(project_dim, input_shape=(emb_dim,))(x)
            return layers.Dropout(rate=dropout)(x)

        return _embedder, output_dim
