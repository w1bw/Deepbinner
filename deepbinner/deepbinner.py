"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Deepbinner/

This file is part of Deepbinner. Deepbinner is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Deepbinner is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Deepbinner.
If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import pathlib
import sys
from .help_formatter import MyParser, MyHelpFormatter
from .version import __version__


def main():
    parser = MyParser(description='Deepbinner: a deep convolutional neural network '
                                  'barcode demultiplexer for Oxford Nanopore reads',
                      formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    classify_subparser(subparsers)
    bin_subparser(subparsers)
    realtime_subparser(subparsers)
    prep_subparser(subparsers)
    balance_subparser(subparsers)
    train_subparser(subparsers)
    refine_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)
    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.subparser_name == 'classify':
        check_classify_and_realtime_arguments(args)
        from .classify import classify
        classify(args)

    if args.subparser_name == 'bin':
        from .bin import bin_reads
        bin_reads(args)

    if args.subparser_name == 'realtime':
        check_classify_and_realtime_arguments(args)
        from .realtime import realtime
        realtime(args)

    elif args.subparser_name == 'prep':
        check_prep_arguments(args)
        from .prep import prep
        prep(args)

    elif args.subparser_name == 'balance':
        check_balance_arguments(args)
        from .balance import balance_training_samples
        balance_training_samples(args)

    elif args.subparser_name == 'train':
        from .train_network import train
        train(args)

    elif args.subparser_name == 'refine':
        from .refine import refine_training_samples
        refine_training_samples(args)


def classify_subparser(subparsers):
    group = subparsers.add_parser('classify', description='Classify fast5 reads',
                                  formatter_class=MyHelpFormatter, add_help=False)

    positional_args = group.add_argument_group('Positional')
    positional_args.add_argument('input', type=str,
                                 help='One of the following: a single fast5 file, a directory of '
                                      'fast5 files (will be searched recursively) or a '
                                      'tab-delimited file of training data')
    classify_and_realtime_options(group)

    other_args = group.add_argument_group('Other')
    other_args.add_argument('--verbose', action='store_true',
                            help='Include the output probabilities for all barcodes in the '
                                 'results (default: just show the final barcode call)')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')


def classify_and_realtime_options(group):
    """
    A few options are used in both the classify and realtime command, so they are described in this
    separate function.
    """
    model_args = group.add_argument_group('Model presets')
    model_args.add_argument('--native', action='store_true',
                            help='Preset for EXP-NBD103 read start and end models')
    model_args.add_argument('--rapid', action='store_true',
                            help='Preset for SQK-RBK004 read start model')

    model_args = group.add_argument_group('Models (at least one is required if not using a preset)')
    model_args.add_argument('-s', '--start_model', type=str, required=False,
                            help='Model trained on the starts of reads')
    model_args.add_argument('-e', '--end_model', type=str, required=False,
                            help='Model trained on the ends of reads')

    barcode_args = group.add_argument_group('Barcoding')
    barcode_args.add_argument('--scan_size', type=float, required=False, default=6144,
                              help="This much of a read's start/end signal will examined for "
                                   "barcode signals")
    barcode_args.add_argument('--score_diff', type=float, required=False, default=0.5,
                              help='For a read to be classified, there must be this much '
                                   'difference between the best and second-best barcode scores')

    two_model_args = group.add_argument_group('Two model (read start and read end) behaviour')
    two_model_args.add_argument('--require_either', action='store_true',
                                help='Most lenient approach: a barcode call on either the start '
                                     'or end is sufficient to classify a read, as long as they do '
                                     'not disagree on the barcode')
    two_model_args.add_argument('--require_start', action='store_true',
                                help='Moderate approach: a start barcode is required to classify '
                                     'a read but an end barcode is optional (default behaviour)')
    two_model_args.add_argument('--require_both', action='store_true',
                                help='Most stringent approach: both start and end barcodes must be '
                                     'present and agree to classify a read')

    perf_args = group.add_argument_group('Performance')
    perf_args.add_argument('--batch_size', type=int, required=False, default=256,
                           help='Neural network batch size')
    perf_args.add_argument('--intra_op_parallelism_threads', type=int, required=False, default=12,
                           help='TensorFlow\'s intra_op_parallelism_threads config option')
    perf_args.add_argument('--inter_op_parallelism_threads', type=int, required=False, default=1,
                           help='TensorFlow\'s inter_op_parallelism_threads config option')
    perf_args.add_argument('--device_count', type=int, required=False, default=1,
                           help='TensorFlow\'s device_count config option')
    perf_args.add_argument('--omp_num_threads', type=int, required=False, default=12,
                           help='OMP_NUM_THREADS environment variable value')


def bin_subparser(subparsers):
    group = subparsers.add_parser('bin', description='Bin fasta/q reads',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required')
    required_args.add_argument('--classes', type=str, required=True,
                               help='Deepbinner classification file (made with the deepbinner '
                                    'classify command)')
    required_args.add_argument('--reads', type=str, required=True,
                               help='FASTA or FASTQ reads')
    required_args.add_argument('--type', type=str, choices=('fastq', 'fasta', 'fast5'),
                               help='type of read file (or directory of files if fast5)')
    required_args.add_argument('--out_dir', type=str, required=True,
                               help='Directory to output binned read files')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')


def realtime_subparser(subparsers):
    group = subparsers.add_parser('realtime', description='Sort fast5 files during sequencing',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required')
    required_args.add_argument('--in_dir', type=str, required=True,
                               help='Directory where sequencer deposits fast5 files')
    required_args.add_argument('--out_dir', type=str, required=True,
                               help='Directory to output binned fast5 files')

    classify_and_realtime_options(group)

    other_args = group.add_argument_group('Other')
    other_args.add_argument('--nowait', action='store_true',
                            help="After processing input directory, don't wait for more files")
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')


def prep_subparser(subparsers):
    group = subparsers.add_parser('prep', description='Prepare training data',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required')
    required_args.add_argument('--fastq', type=str, required=True,
                               help='a FASTQ file of basecalled reads')
    required_args.add_argument('--fast5_dir', type=str, required=True,
                               help='The directory containing the fast5 files (will be searched '
                                    'recursively, so can contain subdirectories)')
    required_args.add_argument('--kit', type=str, required=True,
                               choices=['EXP-NBD103_start', 'EXP-NBD103_end', 'SQK-RBK004_start'],
                               help='Which kit was used to sequence the data')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('--ref', type=str,
                            help='Reference FASTA file (required for EXP-NBD103 kit)')
    other_args.add_argument('--signal_size', type=int, required=False, default=1024,
                            help='Amount of signal (number of samples) that will be used in the '
                                 'neural network')
    other_args.add_argument('--sequencing_summary', type=str, required=False,
                            help='Basecalling sequencing_summary.txt file (if provided, will be '
                                 'used for barcode classification verification)')
    other_args.add_argument('--read_limit', type=int, required=False,
                            help='If provided, will limit the training to this many reads')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')


def balance_subparser(subparsers):
    group = subparsers.add_parser('balance', description='Select balanced training set',
                                  formatter_class=MyHelpFormatter)

    # Positional arguments
    group.add_argument('training_data', type=str, nargs='+',
                       help='Files of raw training data produced by the prep command')

    # Optional arguments
    group.add_argument('--barcodes', type=str,
                       help='A comma-delimited list of which barcodes to include (default: '
                            'include all barcodes)')
    group.add_argument('--random_signal', type=float, required=False, default=1.0,
                       help='This many random signals will be added to the output as no-barcode '
                            'training samples (expressed as a multiple of the balanced '
                            'per-barcode count, default: DEFAULT)')


def train_subparser(subparsers):
    group = subparsers.add_parser('train', description='Train the neural network',
                                  formatter_class=MyHelpFormatter)

    required_args = group.add_argument_group('Required')
    required_args.add_argument('--train', type=str, required=True,
                               help='Balanced training data produced by the balance command')
    required_args.add_argument('--val', type=str, required=True,
                               help='Validation data used to assess the training')
    required_args.add_argument('--model_out', type=str, required=True,
                               help='Filename for the trained model')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('--model_in', type=str,
                            help='An existing model to use as a starting point for training')
    other_args.add_argument('--epochs', type=int, required=False, default=100,
                            help='Number of training epochs')
    other_args.add_argument('--aug', type=float, required=False, default=2.0,
                            help='Data augmentation factor (1 = no augmentation)')
    other_args.add_argument('--batch_size', type=int, required=False, default=20,
                            help='Training batch size')
    other_args.add_argument('--batches_per_epoch', type=int, required=False, default=5000,
                            help='The number of samples per epoch will be this times the batch '
                                 'size')


def refine_subparser(subparsers):
    group = subparsers.add_parser('refine', description='Refine the training set',
                                  formatter_class=MyHelpFormatter)

    # Positional arguments
    group.add_argument('training_data', type=str,
                       help='Balanced training data produced by the balance command')
    group.add_argument('classification_data', type=str,
                       help='Training data barcode calls produced by the classify command')


def check_classify_and_realtime_arguments(args):
    if args.native and args.rapid:
        sys.exit('Error: you can only use one model preset (--native or --rapid)')
    if args.native or args.rapid:
        preset_name = 'native' if args.native else 'rapid'
        if args.start_model is not None or args.end_model is not None:
            sys.exit('Error: you cannot explicitly specify a model and '
                     'also use a model preset (--{})'.format(preset_name))
    if args.native:
        args.start_model = find_native_start_model()
        args.end_model = find_native_end_model()

    if args.rapid:
        args.start_model = find_rapid_start_model()

    model_count = (0 if args.start_model is None else 1) + (0 if args.end_model is None else 1)
    if model_count == 0:
        sys.exit('Error: you must provide at least one model')
    if args.score_diff <= 0.0 or args.score_diff > 1.0:
        sys.exit('Error: --score_diff must be in the range (0, 1] (greater than 0 and less than or '
                 'equal to 1)')

    if model_count < 2 and args.require_either:
        sys.exit('Error: --require_either can only be used with two models (start and end)')
    if model_count < 2 and args.require_start:
        sys.exit('Error: --require_start can only be used with two models (start and end)')
    if model_count < 2 and args.require_both:
        sys.exit('Error: --require_both can only be used with two models (start and end)')

    if two_model_args_used(args) > 1:
        sys.exit('Error: only one of the following options can be used: --require_either, '
                 '--require_start, --require_both')
    if not args.require_either and not args.require_start and not args.require_both:
        args.require_either = True
    assert two_model_args_used(args) == 1


def find_native_start_model():
    return find_model('EXP-NBD103_read_starts')


def find_native_end_model():
    return find_model('EXP-NBD103_read_ends')


def find_rapid_start_model():
    return find_model('SQK-RBK004_read_starts')


def find_model(model_name):
    try:
        start_model = pathlib.Path(__file__).parents[1] / 'models' / model_name
        if start_model.is_file():
            return str(start_model)
    except IndexError:
        pass
    try:
        start_model = pathlib.Path(__file__).parents[0] / 'models' / model_name
        if start_model.is_file():
            return str(start_model)
    except IndexError:
        pass
    sys.exit('Error: could not find {} - did Deepbinner install correctly?'.format(model_name))


def two_model_args_used(args):
    return sum([1 if args.require_either else 0,
                1 if args.require_start else 0,
                1 if args.require_both else 0])


def check_prep_arguments(args):
    if args.kit == 'EXP-NBD103_start' or args.kit == 'EXP-NBD103_end':
        if args.ref is None:
            sys.exit('Error: --ref is required for the EXP-NBD103 kit')


def check_balance_arguments(args):
    if args.barcodes is not None:
        args.barcodes = args.barcodes.split(',')
        try:
            _ = [int(x) for x in args.barcodes]
        except ValueError:
            sys.exit('Error: if used, --barcodes must be a comma-delimited list of numbers (no '
                     'spaces)')


if __name__ == '__main__':
    main()
