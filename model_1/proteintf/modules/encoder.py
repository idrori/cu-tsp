import tensorflow as tf
from tensorflow.keras import layers

from proteintf.common import Registrable
from proteintf.modules.transformer import EncoderStack, get_padding, get_padding_bias


class Encoder(Registrable):
    pass


@Encoder.register('pass_through')
class PassThroughEncoder(Encoder):
    @classmethod
    def from_params(cls, params):
        def _encoder(x, mask=None):
            return x
        return _encoder, params['input_dim']


@Encoder.register('rnn')
class RNNEncoder(Encoder):
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        hid_dim = params['hid_dim']
        num_layers = params.get('num_layers', 1)
        bidirectional = params.get('bidirectional', False)
        batch_norm = params.get('batch_norm', False)
        residual = params.get('residual', False)
        dropout = params.get('dropout', 0.)
        RNN = params.get('RNN', 'LSTM')

        assert not residual or \
           (not bidirectional and input_dim == hid_dim) or \
           (bidirectional and input_dim == hid_dim * 2), \
            'When using residual connections, make sure input_dim == hid_dim * directions, ' \
            f'got input_dim = {input_dim} and hid_dim = {hid_dim}'

        def _encoder(x, mask=None):
            mask = layers.Lambda(
                lambda inp: tf.cast(tf.expand_dims(inp, axis=-1), tf.float32)
            )(mask)
            x = layers.multiply([x, mask])
            x = layers.Masking(mask_value=0., input_shape=(None, input_dim,))(x)
            for i in range(num_layers):
                y = x
                if batch_norm:
                    x = layers.BatchNormalization(axis=-1)(x)
                if 'CuDNN' not in RNN:
                    layer = eval(f'layers.{RNN}')(hid_dim, dropout=dropout, return_sequences=True)
                else:
                    layer = eval(f'layers.{RNN}')(hid_dim, return_sequences=True)
                if bidirectional:
                    dim = input_dim if i == 0 else hid_dim * 2
                    layer = layers.Bidirectional(layer, input_shape=(None, dim,))
                x = layer(x) if 'CuDNN' not in RNN else layer(x, mask=None)
                if residual:
                    x = layers.add([x, y])
            return x

        output_dim = hid_dim * 2 if bidirectional else hid_dim
        return _encoder, output_dim


@Encoder.register('transformer')
class TransformerEncoder(Encoder):
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        num_hidden_layers = params['num_hidden_layers']
        num_heads = params['num_heads']
        filter_size = params['filter_size']
        attention_dropout = params.get('attention_dropout', 0.)
        relu_dropout = params.get('relu_dropout', 0.)
        allow_ffn_pad = params.get('allow_ffn_pad', True)
        postprocess_dropout = params.get('postprocess_dropout', 0.)

        encoder = EncoderStack(num_hidden_layers,
                               input_dim,
                               num_heads,
                               attention_dropout,
                               filter_size,
                               relu_dropout,
                               allow_ffn_pad,
                               postprocess_dropout,
                               train=True)

        def _encoder(x, mask=None):
            if mask is None:
                mask = layers.Lambda(lambda inp: tf.ones_like(inp[:, :, 0]))(x)
            bias = get_padding_bias(mask)
            padding = get_padding(mask)

            mask = layers.Lambda(
                lambda inp: tf.cast(tf.expand_dims(inp, axis=-1), tf.float32)
            )(mask)
            x = layers.multiply([x, mask])
            x = layers.Masking(mask_value=0., input_shape=(None, input_dim,))(x)
            return encoder(x, attention_bias=bias, inputs_padding=padding)

        return _encoder, input_dim
