import tensorflow as tf
from tensorflow.keras import layers

_NEG_INF = -1e9


def get_padding(x, padding_value=0):
    """Return float tensor representing the padding values in x.
    Code taken from https://github.com/tensorflow/models/blob/master/official/transformer/model/model_utils.py#L77

    Args:
        x: int tensor with any shape
        padding_value: int value that
    Returns:
        flaot tensor with same shape as x containing values 0 or 1.
            0 -> non-padding, 1 -> padding
    """
    with tf.name_scope("padding"):
        return tf.to_float(tf.equal(x, padding_value))


def get_padding_bias(x):
    """Calculate bias tensor from padding values in tensor.
    Code taken from https://github.com/tensorflow/models/blob/master/official/transformer/model/model_utils.py#L92

    Bias tensor that is added to the pre-softmax multi-headed attention logits,
    which has shape [batch_size, num_heads, length, length]. The tensor is zero at
    non-padding locations, and -1e9 (negative infinity) at padding locations.
    Args:
        x: int tensor with shape [batch_size, length]
    Returns:
        Attention bias tensor of shape [batch_size, 1, 1, length].
    """
    with tf.name_scope("attention_bias"):
        padding = get_padding(x)
        attention_bias = padding * _NEG_INF
        attention_bias = tf.expand_dims(
                tf.expand_dims(attention_bias, axis=1), axis=1)
    return attention_bias


class PrePostProcessingWrapper(object):
    """
    Code taken from https://github.com/tensorflow/models/blob/master/official/transformer/model/transformer.py#L267
    Wrapper class that applies layer pre-processing and post-processing.
    """

    def __init__(self, layer, layer_postprocess_dropout, train):
        self.layer = layer
        self.postprocess_dropout = layer_postprocess_dropout
        self.train = train

        # Create normalization layer
        self.layer_norm = layers.BatchNormalization()

    def __call__(self, x, *args, **kwargs):
        # Preprocessing: apply layer normalization
        y = self.layer_norm(x)

        # Get layer output
        y = self.layer(y, *args, **kwargs)

        # Postprocessing: apply dropout and residual connection
        if self.train:
            y = tf.nn.dropout(y, 1 - self.postprocess_dropout)
        return x + y
