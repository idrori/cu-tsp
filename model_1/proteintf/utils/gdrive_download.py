import os
import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description='Download files from google drive using gdrive')
    parser.add_argument('file_ids_path', type=str,
                        help='Path to file where each line contains tab-separated file id and file name')
    parser.add_argument('-o', '--output-path', type=str, default='./')
    parser.add_argument('-g', '--gdrive-path', type=str, default='/home/Wayne/gdrive')
    return parser.parse_args()


def run(args):
    file_ids_path = args.file_ids_path
    output_path = args.output_path
    gdrive_path = args.gdrive_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    print('Reading file ids from {}'.format(file_ids_path))
    print('Downloading to {}'.format(output_path))

    failed = []
    with open(file_ids_path, 'r') as fin:
        total = len([_ for _ in fin])
        fin.seek(0)

        for idx, line in enumerate(fin):
            if not line:
                msg = '({}/{}) Skipping empty line'.format(idx + 1, total)
                failed.append(msg)
                print(msg)
                continue

            file_id, file_name = line.strip().split('\t')
            success = False
            while not success:
                result = subprocess.run(
                    [gdrive_path, 'download', '--path', output_path, file_id],
                    stdout=subprocess.PIPE
                )
                outputs = result.stdout.decode('utf-8')
                if not 'rateLimitExceeded' in outputs:
                    success = True
                    if 'Downloaded' in outputs:
                        msg = '({}/{}) {} downloaded'.format(idx + 1, total, file_name)
                    else:
                        msg = '({}/{}) {}'.format(idx + 1, total, ' '.join(outputs))
                        failed.append(msg)
                    print(msg)
    print('\nAll done.')
    print('{} lines failed'.format(len(failed)))
    for msg in failed:
        print(msg)
    print()


if __name__ == '__main__':
    args = parse_args()
    run(args)
