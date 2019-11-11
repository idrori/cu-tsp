import os
import pickle
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Make output pickle')
    parser.add_argument('--testpkl-path', type=str)
    parser.add_argument('--dcalphas-path', type=str)
    parser.add_argument('--psis-path', type=str)
    parser.add_argument('--phis-path', type=str)
    parser.add_argument('--output-dir', type=str)
    return parser.parse_args()


def run(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_path = os.path.join(args.output_dir, 'output.pkl')

    with open(args.testpkl_path, 'rb') as fin:
        pdbs = pickle.load(fin)[1]
    with open(args.dcalphas_path, 'rb') as fin:
        dcalphas = pickle.load(fin)
    with open(args.psis_path, 'rb') as fin:
        psis = pickle.load(fin)
    with open(args.phis_path, 'rb') as fin:
        phis = pickle.load(fin)

    output = {}
    for pdb, d, ps, ph in zip(pdbs, dcalphas, psis, phis):
        output[pdb] = {'dcalphas': d, 'psi': ps, 'phi': ph}

    with open(output_path, 'wb') as fout:
        pickle.dump(output, fout)


if __name__ == '__main__':
    args = parse_args()
    run(args)
