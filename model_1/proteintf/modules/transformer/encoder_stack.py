"""
Code taken from https://github.com/tensorflow/models/blob/master/official/transformer/model/transformer.py#L291
"""

import tensorflow as tf
from tensorflow.keras import layers

from proteintf.modules.transformer.attention_layer import SelfAttention
from proteintf.modules.transformer.utils import PrePostProcessingWrapper
from proteintf.modules.transformer.feedforward_network import FeedFowardNetwork


class EncoderStack(layers.Layer):
    """Transformer encoder stack.
    The encoder stack is made up of N identical layers. Each layer is composed
    of the sublayers:
      1. Self-attention layer
      2. Feedforward network (which is 2 fully-connected layers)
    """
    def __init__(self, num_hidden_layers, hidden_size, num_heads, attention_dropout,
                 filter_size, relu_dropout, allow_ffn_pad, postprocess_dropout, train):
        super(EncoderStack, self).__init__()
        self.train = train

        self.layers = []
        for _ in range(num_hidden_layers):
            # Create sublayers for each layer.
            self_attention_layer = SelfAttention(hidden_size, num_heads, attention_dropout, train)
            feed_forward_network = FeedFowardNetwork(hidden_size, filter_size, relu_dropout, train, allow_ffn_pad)

            self.layers.append([
                PrePostProcessingWrapper(self_attention_layer, postprocess_dropout, train),
                PrePostProcessingWrapper(feed_forward_network, postprocess_dropout, train)])

        self.output_normalization = layers.BatchNormalization()

    def call(self, encoder_inputs, attention_bias, inputs_padding):
        """Return the output of the encoder layer stacks.
        Args:
          encoder_inputs: tensor with shape [batch_size, input_length, hidden_size]
          attention_bias: bias for the encoder self-attention layer.
            [batch_size, 1, 1, input_length]
          inputs_padding: P
        Returns:
          Output of encoder layer stack.
          float32 tensor with shape [batch_size, input_length, hidden_size]
        """
        for n, layer in enumerate(self.layers):
            # Run inputs through the sublayers.
            self_attention_layer = layer[0]
            feed_forward_network = layer[1]

            with tf.variable_scope("layer_%d" % n):
                with tf.variable_scope("self_attention"):
                    encoder_inputs = self_attention_layer(encoder_inputs, bias=attention_bias)
                with tf.variable_scope("ffn"):
                    encoder_inputs = feed_forward_network(encoder_inputs, inputs_padding)

            return self.output_normalization(encoder_inputs)
