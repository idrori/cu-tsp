import tensorflow as tf
from tensorflow.keras import layers

from proteintf.common import Registrable
from proteintf.modules.feedforward import FeedForward


class Decoder(Registrable):
    pass


@Decoder.register('pairwise_bilinear')
class PairwiseBilinearDecoder(Decoder):
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        name = params.get('name', None)

        def _decoder(x, y=None):
            if y is None:
                y = x
            y = layers.Dense(input_dim, input_shape=(input_dim,))(y)
            return layers.Lambda(lambda inp: tf.matmul(inp[0], inp[1], transpose_b=True), name=name)((x, y))
        return _decoder


@Decoder.register('pairwise_linear')
class PairwiseLinearDecoder(Decoder):
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        name = params.get('name', None)

        def _decoder(x, y=None):
            if y is None:
                y = x
            x = layers.Dense(1, input_shape=(input_dim,))(x)
            y = layers.Dense(1, input_shape=(input_dim,))(y)
            y = layers.Lambda(lambda inp: tf.expand_dims(tf.squeeze(inp, axis=-1), axis=1))(y)
            return layers.add([x, y], name=name)
        return _decoder


@Decoder.register('pairwise_dot')
class PairwiseDotDecoder(Decoder):
    @classmethod
    def from_params(cls, params):
        name = params.get('name', None)

        def _decoder(x, y=None):
            if y is None:
                y = x
            return layers.Lambda(lambda inp: tf.matmul(inp[0], inp[1], transpose_b=True), name=name)((x, y))
        return _decoder


@Decoder.register('elementwise_linear')
class ElementwiseLinearDecoder(Decoder):
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        name = params.get('name', None)
        if 'feedforward' in params:
            feedforward_config = params['feedforward']
            feedforward_config['input_dim'] = input_dim
            feedforward, output_dim = FeedForward.from_params(feedforward_config)
            input_dim = output_dim
        else:
            feedforward = lambda inp: inp

        def _decoder(x):
            x = feedforward(x)
            x = layers.Dense(1, input_shape=(input_dim,))(x)
            return layers.Lambda(lambda inp: tf.squeeze(inp, axis=-1),name=name)(x)
        return _decoder
