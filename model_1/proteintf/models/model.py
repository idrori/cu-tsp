from copy import deepcopy
import tensorflow as tf

from proteintf.common import Registrable
from proteintf.modules import ProteinEmbedder, Encoder, FeedForward, Decoder, Loss


class Model(Registrable):
    pass


@Model.register('single_seq_model')
class SingleSeqModel(Model):
    @classmethod
    def from_params(cls, params):
        params = deepcopy(params)
        input_names = params['input_names']
        target_names = params['target_names']
        embedder_config = params['embedder']
        encoder_config = params['encoder']
        decoder_config = params['decoder']
        loss_config = params['loss']
        optimizer_config = params['optimizer']
        metrics_config = params.get('metrics', {})
        target_stats = params.get('target_stats', None)
        msas_dim = params.get('msas_dim', 21)
        window_size = params.get('window_size', None)

        mask = tf.keras.Input(shape=(window_size,), name='mask', dtype=tf.int64)
        inputs = []
        for n in input_names:
            if n != 'msas':
                inputs += [tf.keras.Input(shape=(window_size,), name=n, dtype=tf.int64)]
            else:
                inputs += [tf.keras.Input(shape=(window_size, msas_dim,), name=n, dtype=tf.float32)]

        if 'position' in embedder_config['embedders']:
            input_names = ['position'] + input_names
            position_emb_dim = embedder_config['embedders']['position']['input_dim']
            inputs = [tf.keras.Input(shape=(window_size, position_emb_dim,), name='position_emb', dtype=tf.float32)] + inputs

        embedder_config['input_names'] = input_names
        embedder, output_dim = ProteinEmbedder.from_params(embedder_config)

        encoder_config['input_dim'] = output_dim
        encoder, output_dim = Encoder.by_name(encoder_config['type']).from_params(encoder_config)

        if 'feedforward' in params:
            feedforward_config = params['feedforward']
            feedforward_config['input_dim'] = output_dim
            feedforward, output_dim = FeedForward.from_params(feedforward_config)
        else:
            feedforward = lambda inp: inp

        decoders, losses, metrics = [], [], []
        for n in target_names:
            decoder_config[n]['input_dim'] = output_dim
            decoders += [Decoder.by_name(decoder_config[n]['type']).from_params(decoder_config[n])]
            losses += [Loss.by_name(loss_config[n]['type']).from_params(loss_config[n])]
            if n in metrics_config:
                if target_stats is not None:
                    metrics_config[n]['target_stats'] = target_stats[n]
                metrics += [Loss.by_name(metrics_config[n]['type']).from_params(metrics_config[n])]

        x = embedder(inputs)
        x = encoder(x, mask)
        x = feedforward(x)
        x = [decoder_model(x) for decoder_model in decoders]

        model = tf.keras.Model(inputs=[mask] + inputs, outputs=x)

        optimizer_type = optimizer_config.pop('type')
        optimizer = eval('tf.keras.optimizers.{}'.format(optimizer_type))(**optimizer_config)
        model.compile(optimizer,
                      loss={n: loss(mask=mask) for n, loss in zip(target_names, losses)},
                      metrics={n: metric(mask=mask) for n, metric in zip(target_names, metrics)})

        return model
