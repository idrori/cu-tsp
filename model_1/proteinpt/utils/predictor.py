import torch
import pickle
import argparse
import numpy as np

from allennlp.common.tqdm import Tqdm
from allennlp.data import DatasetReader
from allennlp.data.iterators import BasicIterator
from allennlp.nn.util import move_to_device
from allennlp.models.archival import load_archive

from proteinpt import *


def parse_args():
    parser = argparse.ArgumentParser(description='Make predictions using a trained model')
    parser.add_argument('--model-path', type=str)
    parser.add_argument('--input-path', type=str)
    parser.add_argument('--output-path', type=str)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--device', type=str)
    return parser.parse_args()


def run(args):
    print('\nArguments:')
    for k, v in vars(args).items():
        print('{}: {}'.format(k, v))
    print()

    device = args.device
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'

    print('Loading archive ...')
    archive = load_archive(args.model_path)
    # predictor = Predictor.from_archive(archive, 'protein_predictor')
    config = archive.config.duplicate()
    dataset_reader = DatasetReader.from_params(config["dataset_reader"])
    model = archive.model.to(device).eval()

    print('Loading data ...')
    dataset_reader.lazy = False
    dataset = dataset_reader.read(args.input_path)
    iterator = BasicIterator(args.batch_size)
    iterator.index_with(model.vocab)
    num_batches = iterator.get_num_batches(dataset)
    data_generator = iterator(dataset, num_epochs=1, shuffle=False)

    print('Predicting ...')
    output_dict = {}
    with torch.no_grad():
        for batch in Tqdm.tqdm(data_generator, total=num_batches):
            batch = move_to_device(batch, model._get_prediction_device())
            outputs = model(**batch)
            predictions = outputs['predictions'].cpu().numpy()
            for pid, length, pred in zip(outputs['protein_id'], outputs['length'], predictions):
                if model.target == 'dcalpha':
                    dcalpha = pred[:length, :length]
                    dcalpha = np.triu(dcalpha, 1) + np.tril(dcalpha.transpose(), -1)
                    output_dict[pid] = {'dcalpha': dcalpha}
                elif model.target == 'angles':
                    psi, phi = pred[:length, 0], pred[:length, 1]
                    # psi[0] = 0.
                    # phi[-1] = 0.
                    output_dict[pid] = {'psi': psi, 'phi': phi}
                else:
                    coords = pred[:length]
                    output_dict[pid] = {'coords': coords}

    print('Writing to {}'.format(args.output_path))
    with open(args.output_path, 'wb') as fout:
        pickle.dump(output_dict, fout)

    print('All done.')


if __name__ == '__main__':
    args = parse_args()
    run(args)
