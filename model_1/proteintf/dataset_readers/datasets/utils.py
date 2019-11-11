import math
import numpy as np
import tensorflow as tf


def get_output_info(inputs, targets, msas_dim, position_emb_dim, window_size=None, masking=True, test=False):
    output_types, output_shapes, padded_shapes = [], [], []
    if masking:
        output_types = [tf.int64]
        output_shapes = [tf.TensorShape([None])]
        padded_shapes = [(window_size,)]

    if position_emb_dim is not None:
        output_types += [tf.float32]
        output_shapes += [tf.TensorShape([None, position_emb_dim])]
        padded_shapes += [(window_size, position_emb_dim,)]

    for col in inputs:
        if col != 'msas':
            output_types += [tf.int64]
            output_shapes += [tf.TensorShape([None])]
            padded_shapes += [(window_size,)]
        else:
            output_types += [tf.float32]
            output_shapes += [tf.TensorShape([None, msas_dim])]
            padded_shapes += [(window_size, msas_dim,)]

    if not test:
        target_types, target_shapes, target_padded_shapes = [], [], []
        for col in targets:
            target_types += [tf.float32]
            if col == 'dcalphas':
                target_shapes += [tf.TensorShape([None, None])]
                target_padded_shapes += [(window_size, window_size,)]
            else:
                target_shapes += [tf.TensorShape([None])]
                target_padded_shapes += [(window_size,)]
        output_types = (tuple(output_types), tuple(target_types))
        output_shapes = (tuple(output_shapes), tuple(target_shapes))
        padded_shapes = (tuple(padded_shapes), tuple(target_padded_shapes))
    else:
        output_types = tuple(output_types)
        output_shapes = tuple(output_shapes)
        padded_shapes = tuple(padded_shapes)

    return output_types, output_shapes, padded_shapes


def get_position_encoding(
        length, hidden_size, min_timescale=1.0, max_timescale=1.0e4):
    """Return positional encoding.
    Code taken from https://github.com/tensorflow/models/blob/master/official/transformer/model/model_utils.py#L28

    Calculates the position encoding as a mix of sine and cosine functions with
    geometrically increasing wavelengths.
    Defined and formulized in Attention is All You Need, section 3.5.
    Args:
        length: Sequence length.
        hidden_size: Size of the
        min_timescale: Minimum scale that will be applied at each position
        max_timescale: Maximum scale that will be applied at each position
    Returns:
        Tensor with shape [length, hidden_size]
    """
    position = np.arange(length).astype(np.float32)
    num_timescales = hidden_size // 2
    log_timescale_increment = (
            math.log(float(max_timescale) / float(min_timescale)) /
            (float(num_timescales) - 1))
    inv_timescales = min_timescale * np.exp(
            np.arange(num_timescales).astype(np.float32) * -log_timescale_increment)
    scaled_time = position[:, np.newaxis] * inv_timescales[np.newaxis, :]
    signal = np.concatenate([np.sin(scaled_time), np.cos(scaled_time)], axis=1)
    return signal
