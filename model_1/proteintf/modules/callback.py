import tensorflow as tf
from tensorflow.keras import backend as K

from proteintf.common import Registrable


class Callback(Registrable):
    pass


@Callback.register('model_checkpoint')
class CheckPointCallback(Callback):
    @classmethod
    def from_params(cls, params):
        return tf.keras.callbacks.ModelCheckpoint(**params)


@Callback.register('early_stopping')
class EarlyStoppingCallback(Callback):
    @classmethod
    def from_params(cls, params):
        return tf.keras.callbacks.EarlyStopping(**params)


@Callback.register('tensorboard')
class TensorBoardCallback(Callback):
    @classmethod
    def from_params(cls, params):
        return tf.keras.callbacks.TensorBoard(**params)


@Callback.register('noam_lr_scheduler')
class NoamLearningRateScheduler(tf.keras.callbacks.Callback, Callback):
    def __init__(self, model_size, warmup_steps, factor=1., verbose=0):
        super(NoamLearningRateScheduler, self).__init__()

        self.model_size = model_size
        self.warmup_steps = warmup_steps
        self.factor = factor
        self.verbose = verbose

        self.iterations = 0
        self.current_lr = None

    def on_train_begin(self, logs):
        self.iterations = 0

    def on_batch_begin(self, batch, logs):
        self.iterations += 1
        if not hasattr(self.model.optimizer, 'lr'):
            raise ValueError('Optimizer must have a "lr" attribute.')
        lr = self.factor * (self.model_size ** (-0.5) *
                            min(self.iterations ** (-0.5),
                                self.iterations * self.warmup_steps ** (-1.5)))
        K.set_value(self.model.optimizer.lr, lr)
        self.current_lr = lr

    def on_epoch_begin(self, epoch, logs=None):
        if self.verbose > 0 and epoch != 0:
            print('NoamLRScheduler: Current learning rate: %s.' % self.current_lr)

    @classmethod
    def from_params(cls, params):
        model_size = params['model_size']
        warmup_steps = params['warmup_steps']
        factor = params.get('factor', 1.)
        verbose = params.get('verbose', 0)
        return NoamLearningRateScheduler(model_size, warmup_steps, factor, verbose)
