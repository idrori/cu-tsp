import tensorflow as tf
from proteintf.common import Registrable


class Loss(Registrable):
    pass


@Loss.register('rmse')
class RootedMeanSquaredErrorLoss(Registrable):
    @classmethod
    def from_params(cls, params):
        name = params.get('name', 'rmse')
        target_stats = params.get('target_stats', None)

        def _rmse_loss(mask=None):
            def _loss(y_true, y_pred):
                diff = y_pred - y_true
                if target_stats is not None:
                    diff *= target_stats['std']
                m = mask
                if m is not None:
                    dim = params['dim']
                    assert dim == 2 or dim == 3
                    if dim == 3:
                        m = tf.expand_dims(m, axis=-1)
                        m = tf.matmul(m, m, transpose_b=True)
                        m_tril = 1 - tf.matrix_band_part(tf.ones_like(m), -1, -0)
                        m *= m_tril
                    m = tf.cast(m, tf.float32)
                    return tf.sqrt(tf.reduce_sum((diff ** 2) * m) / tf.reduce_sum(m))
                else:
                    return tf.sqrt(tf.reduce_mean(diff ** 2))
            _loss.__name__ = name
            return _loss
        return _rmse_loss


@Loss.register('mse')
class MeanSquaredErrorLoss(Registrable):
    @classmethod
    def from_params(cls, params):
        name = params.get('name', 'mse')
        target_stats = params.get('target_stats', None)

        def _mse_loss(mask=None):
            def _loss(y_true, y_pred):
                diff = y_pred - y_true
                if target_stats is not None:
                    diff *= target_stats['std']
                m = mask
                if m is not None:
                    dim = params['dim']
                    assert dim == 2 or dim == 3
                    if dim == 3:
                        m = tf.expand_dims(m, axis=-1)
                        m = tf.matmul(m, m, transpose_b=True)
                        m_tril = 1 - tf.matrix_band_part(tf.ones_like(m), -1, -0)
                        m *= m_tril
                    m = tf.cast(m, tf.float32)
                    return tf.reduce_sum((diff ** 2) * m) / tf.reduce_sum(m)
                else:
                    return tf.reduce_mean(diff ** 2)
            _loss.__name__ = name
            return _loss
        return _mse_loss


@Loss.register('mae')
class MeanAbsoluteErrorLoss(Registrable):
    @classmethod
    def from_params(cls, params):
        name = params.get('name', 'mae')
        target_stats = params.get('target_stats', None)

        def _mae_loss(mask=None):
            def _loss(y_true, y_pred):
                diff = y_pred - y_true
                if target_stats is not None:
                    diff *= target_stats['std']
                m = mask
                if m is not None:
                    dim = params['dim']
                    assert dim == 2 or dim == 3
                    if dim == 3:
                        m = tf.expand_dims(m, axis=-1)
                        m = tf.matmul(m, m, transpose_b=True)
                        m_tril = 1 - tf.matrix_band_part(tf.ones_like(m), -1, -0)
                        m *= m_tril
                    m = tf.cast(m, tf.float32)
                    return tf.reduce_sum(tf.abs(diff) * m) / tf.reduce_sum(m)
                else:
                    return tf.reduce_mean(tf.abs(diff))
            _loss.__name__ = name
            return _loss
        return _mae_loss
