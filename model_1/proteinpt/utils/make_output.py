import os
import pickle
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Make output pickle')
    parser.add_argument('--dcalphas-path', type=str)
    parser.add_argument('--angles-path', type=str)
    parser.add_argument('--coords-path', type=str)
    parser.add_argument('--output-path', type=str)
    return parser.parse_args()


def run(args):
    output_dir = os.path.split(args.output_path)[0]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(args.dcalphas_path, 'rb') as fin:
        dcalphas = pickle.load(fin)
    with open(args.angles_path, 'rb') as fin:
        angles = pickle.load(fin)
    coords = None
    if args.coords_path is not None:
        with open(args.coords_path, 'rb') as fin:
            coords = pickle.load(fin)

    output = {}
    for pid in dcalphas:
        output[pid] = {'dcalphas': dcalphas[pid]['dcalpha'], 'psi': angles[pid]['psi'], 'phi': angles[pid]['phi']}
        if coords is not None:
            output[pid]['coords'] = coords[pid]['coords']

    with open(args.output_path, 'wb') as fout:
        pickle.dump(output, fout)


if __name__ == '__main__':
    args = parse_args()
    run(args)
