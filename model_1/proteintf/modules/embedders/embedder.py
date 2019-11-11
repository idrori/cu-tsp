import tensorflow as tf
from tensorflow.keras import layers

from proteintf.common import Registrable


class Embedder(Registrable):
    pass


@Embedder.register('embedding')
class Embedding(Embedder):
    @classmethod
    def from_params(cls, params):
        emb_dim = params['emb_dim']
        vocab_size = params['vocab_size']
        mask_zero = False
        if params['mask_zero']:
            mask_zero = True
            vocab_size += 1
        embedder = [layers.Embedding(input_dim=vocab_size,
                                     output_dim=emb_dim,
                                     mask_zero=mask_zero,
                                     embeddings_initializer='glorot_normal')]
        if 'project_dim' in params and params['project_dim'] is not None:
            project_dim = params['project_dim']
            embedder += [layers.Dense(project_dim, input_shape=(emb_dim,))]
            emb_dim = project_dim
        return embedder, emb_dim


@Embedder.register('one_hot')
class OneHotEmbedder(Embedder):
    @classmethod
    def from_params(cls, params):
        depth = params['vocab_size']
        if params['mask_zero']:
            depth += 1
        embedder = [layers.Lambda(lambda inp: tf.one_hot(inp, depth))]
        return embedder, depth


@Embedder.register('dense')
class DenseEmbedder(Embedder):
    @classmethod
    def from_params(cls, params):
        embedder = []
        if 'batch_norm' in params and params['batch_norm']:
            embedder += [layers.BatchNormalization(axis=-1)]
        embedder += [layers.Dense(params['output_dim'], input_shape=(params['input_dim'],))]
        return embedder, params['output_dim']


@Embedder.register('pass_through')
class PassThroughEmbedder(Embedder):
    @classmethod
    def from_params(cls, params):
        if 'batch_norm' in params and params['batch_norm']:
            embedder = [layers.BatchNormalization(axis=-1)]
        else:
            embedder = [lambda inp: inp]
        return embedder, params['input_dim']
