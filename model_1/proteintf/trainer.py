import os
import json
import shutil
import argparse

from proteintf.models import Model
from proteintf.configs import CONFIGS
from proteintf.modules import Callback
from proteintf.dataset_readers import ProteinDatasetReader
from proteintf.utils import with_fallback, prepare_environment, prepare_global_logging, cleanup_global_logging


def parse_args():
    parser = argparse.ArgumentParser(description='Train a model for protein tertiary structure prediction')
    parser.add_argument('config', type=str)
    parser.add_argument('-s', '--serialization-dir', type=str)
    parser.add_argument('-o', '--overrides', type=str)
    parser.add_argument('-d', '--seed', type=int, default=3435)
    return parser.parse_args()


def run(args):
    prepare_environment(args.seed)

    def _print_config(c, indent):
        def _print_with_indent(msg):
            print('{}{}'.format(' ' * indent, msg))
        for k, v in c.items():
            if isinstance(v, dict):
                _print_with_indent(k + ':')
                _print_config(v, indent + 2)
            else:
                _print_with_indent('{}: {}'.format(k, v))

    config = CONFIGS[args.config]
    if args.overrides is not None:
        config_overrides = json.loads(args.overrides)
        config = with_fallback(preferred=config_overrides, fallback=config)

    stdout_handler = None
    if args.serialization_dir is not None:
        output_path = args.serialization_dir
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path)
        stdout_handler = prepare_global_logging(output_path, log_prefix='train_')

        if 'model_checkpoint' in config['callbacks']:
            config['callbacks']['model_checkpoint']['save_weights_only'] = True
            config['callbacks']['model_checkpoint']['filepath'] = os.path.join(output_path, 'best.hdf5')    # model_e{epoch:02d}.hdf5
        if 'tensorboard' in config['callbacks']:
            config['callbacks']['tensorboard']['log_dir'] = os.path.join(output_path, 'log')

        with open(os.path.join(output_path, 'config.json'), 'w') as fout:
            json.dump(config, fout)

    # print('Configs: {}\n'.format(config))
    _print_config({'configs': config}, 0)
    print()

    input_names = config['inputs']
    target_names = config['targets']

    print('Making datasets ...\n')
    reader_config = config['dataset_reader']
    root = reader_config['root']
    train_file_prefix = reader_config.get('train_file_prefix', 'train_fold_')
    total_folds = reader_config.get('total_folds', 10)
    reader = ProteinDatasetReader(root, train_file_prefix, total_folds)
    reader_params = reader_config['params']
    reader.set_params(
        inputs=tuple(input_names),
        targets=tuple(target_names),
        position_emb_dim=reader_params['position_emb_dim'],
        window_size=reader_params['window_size']
    )

    train_dataset = reader.make_train_dataset(folds=reader_config['train_folds'],
                                              batch_size=reader_config['batch_size'],
                                              max_length=reader_config['max_length'],
                                              normalize_target=reader_config['normalize_target'],
                                              shuffle=reader_config['shuffle'])
    valid_dataset = reader.make_valid_dataset(folds=reader_config['valid_folds'],
                                              batch_size=reader_config['batch_size'],
                                              max_length=reader_config['max_length'],
                                              normalize_target=reader_config['normalize_target'],
                                              shuffle=reader_config['shuffle'])
    print('Done\n')

    print('Initializing model ...\n')
    model_config = config['model']
    model_config['input_names'] = input_names
    model_config['target_names'] = target_names
    if reader_config['normalize_target']:
        model_config['target_stats'] = reader.get_target_stats()
    model = Model.by_name(model_config['type']).from_params(model_config)
    print('Done\n')

    print('Training ...\n')
    model.fit(
        x=train_dataset.dataset,
        epochs=config['num_epochs'],
        steps_per_epoch=len(train_dataset),
        validation_data=valid_dataset.dataset,
        validation_steps=len(valid_dataset),
        callbacks=[Callback.by_name(n).from_params(p) for n, p in config['callbacks'].items()],
        verbose=2
    )

    if stdout_handler is not None:
        cleanup_global_logging(stdout_handler)
    print('All done\n')


if __name__ == '__main__':
    args = parse_args()
    run(args)
