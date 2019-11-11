import os
import json
import pickle
import argparse

from proteintf.models import Model
from proteintf.dataset_readers import ProteinDatasetReader
from proteintf.utils import prepare_environment, prepare_global_logging, cleanup_global_logging


def parse_args():
    parser = argparse.ArgumentParser(description='Make predictions using a trained model')
    parser.add_argument('model_path', type=str)
    parser.add_argument('-m', '--model-file', type=str, default='best.hdf5')
    parser.add_argument('-d', '--seed', type=int, default=3435)
    return parser.parse_args()


def run(args):
    prepare_environment(args.seed)
    stdout_handler = prepare_global_logging(args.model_path, log_prefix='test_')

    print('Loading config.json ...\n')
    with open(os.path.join(args.model_path, 'config.json'), 'r') as fin:
        config = json.load(fin)
    print('Done\n')

    input_names = config['inputs']
    target_names = config['targets']

    print('Making test dataset ...\n')
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

    if reader_config['normalize_target']:
        # Need this call to get target_stats
        reader.make_train_dataset(folds=reader_config['train_folds'],
                                  batch_size=reader_config['batch_size'],
                                  max_length=reader_config['max_length'],
                                  normalize_target=reader_config['normalize_target'],
                                  shuffle=reader_config['shuffle'])
    test_dataset = reader.make_test_dataset(batch_size=reader_config['batch_size'])
    print('Done\n')

    print('Loading model ...\n')
    model_config = config['model']
    model_config['input_names'] = input_names
    model_config['target_names'] = target_names
    if reader_config['normalize_target']:
        model_config['target_stats'] = reader.get_target_stats()
    model = Model.by_name(model_config['type']).from_params(model_config)
    model.load_weights(os.path.join(args.model_path, args.model_file))
    print('Done\n')

    print('Making predictions ...\n')
    iterator = test_dataset.dataset.make_one_shot_iterator()
    next_element = iterator.get_next()
    predictions = model.predict(next_element, steps=len(test_dataset))
    outputs = test_dataset.format_predictions(predictions)
    print('Done\n')

    output_path = os.path.join(args.model_path, 'predictions')
    print('Writing outputs to {}\n'.format(output_path))
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for tgt in outputs:
        with open(os.path.join(output_path, tgt + '.pkl'), 'wb') as fout:
            pickle.dump(outputs[tgt], fout)

    cleanup_global_logging(stdout_handler)
    print('All done\n')


if __name__ == '__main__':
    args = parse_args()
    run(args)
