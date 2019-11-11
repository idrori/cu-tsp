import os
import pickle
import argparse
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description='Make plots of distance matrices')
    parser.add_argument('testpkl_path', type=str, help='Path to test.pkl')
    parser.add_argument('dcalphas_path', type=str, help='Path to dcalphas.pkl')
    parser.add_argument('output_dir', type=str, help='Path to the directory where the plots will be saved')
    return parser.parse_args()


def run(args):
    print('Reading distance matrices from {}'.format(args.dcalphas_path))
    with open(args.dcalphas_path, 'rb') as fin:
        dcalphas = pickle.load(fin)
    print('Reading test file from {}'.format(args.testpkl_path))
    with open(args.testpkl_path, 'rb') as fin:
        names = pickle.load(fin)[1]

    print('Saving plots to {}'.format(args.output_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for i in range(len(names)):
        plt.imshow(dcalphas[i], interpolation='nearest')
        plt.title(names[i])
        plt.savefig(os.path.join(args.output_dir, '{}_{}.png'.format(i+1, names[i])))

    print('All done.')


if __name__ == '__main__':
    args = parse_args()
    run(args)
