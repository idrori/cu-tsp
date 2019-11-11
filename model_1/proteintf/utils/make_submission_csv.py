import os
import pickle
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Make csv submission file for Kaggle')
    parser.add_argument('--testpkl-path', type=str)
    parser.add_argument('--dcalphas-path', type=str)
    parser.add_argument('--psis-path', type=str)
    parser.add_argument('--phis-path', type=str)
    parser.add_argument('--output-dir', type=str)
    return parser.parse_args()


def run(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_path = os.path.join(args.output_dir, 'submission.csv')

    with open(args.testpkl_path, 'rb') as fin:
        pdbs = pickle.load(fin)[1]
    with open(args.dcalphas_path, 'rb') as fin:
        dcalphas = pickle.load(fin)
    with open(args.psis_path, 'rb') as fin:
        psis = pickle.load(fin)
    with open(args.phis_path, 'rb') as fin:
        phis = pickle.load(fin)

    tb = []
    for pdb, d, ps, ph in zip(pdbs, dcalphas, psis, phis):
        for i in range(len(d)):
            for j in range(len(d[i])):
                tb += [('{}_d_{}_{}'.format(pdb, i + 1, j + 1), d[i][j])]
        for i in range(len(ps)):
            tb += [('{}_psi_{}'.format(pdb, i + 1), ps[i])]
        for i in range(len(ph)):
            tb += [('{}_phi_{}'.format(pdb, i + 1), ph[i])]
    tb = pd.DataFrame(tb, columns=['Id', 'Predicted'])
    tb.to_csv(output_path, index=False)


if __name__ == '__main__':
    args = parse_args()
    run(args)
