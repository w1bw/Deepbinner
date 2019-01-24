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

# Silence some warnings so the output isn't cluttered.
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
warnings.simplefilter('ignore', FutureWarning)

import os
import h5py
import sys
import gzip


def get_read_id(h5readgroup):
    return h5readgroup.attrs['read_id'].decode()

def get_signal(h5readgroup):
    return h5readgroup['Signal'][:]

def get_read_id_and_signal_from_read_group(h5readgroup):
    read_id = h5readgroup.attrs['read_id'].decode()
    signal = h5readgroup['Signal'][:]
    return get_read_id(h5readgroup), get_signal(h5readgroup)

def get_read_id_and_signal(fast5_file):
    try:
        with h5py.File(fast5_file, 'r') as hdf5_file:
            read_group = list(hdf5_file['Raw/Reads/'].values())[0]
            return get_read_id_and_signal_from_read_group(read_group)
    except OSError:
        return None, None


def find_all_fast5s(directory, verbose=False):
    if not os.path.isdir(directory):
        return [directory]
    if verbose:
        print('Looking for fast5 files in {}... '.format(directory), file=sys.stderr, end='',
              flush=True)
    fast5s = []
    for dir_name, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fast5'):
                fast5s.append(os.path.join(dir_name, filename))
    if verbose:
        noun = 'fast5' if len(fast5s) == 1 else 'fast5s'
        print('{} {} found'.format(len(fast5s), noun), file=sys.stderr)
    return fast5s


def find_fast5_reads(fast5s):
    for fast5 in fast5s:
        with h5py.File(fast5, 'r') as hdf5_file:
            if 'Raw' in hdf5_file:
                # this is a single-read fasta file
                read_group = list(hdf5_file['Raw/Reads/'].values())[0]
                yield hdf5_file, read_group, get_read_id(read_group), fast5, False
            else:
                # In multi-fast5 file, each read is in a top-level group, and the read data are one level down in "Raw"
                for read in hdf5_file.values():
                    # this is a single-read fasta file
                    read_group = read['Raw']
                    yield read, read_group, get_read_id(read_group), fast5, True


def find_fast5_read_ids_and_signals(fast5s):
    for read, group, id, fast5, multi in find_fast5_reads(fast5s):
        signal = get_signal(group)
        yield id, signal, fast5

def index_fast5_read_ids(fast5s):
    index = {}
    for fast5 in fast5s:
        with h5py.File(fast5, 'r') as hdf5_file:
            if 'Raw' in hdf5_file:
                # this is a single-read fasta file
                read_group = list(hdf5_file['Raw/Reads/'].values())[0]
                read_id = get_read_id(read_group)
                index[read_id] = (fast5, False)
            else:
                # In multi-fast5 file, each read is in a top-level group, and the read data are one level down in "Raw"
                for read in hdf5_file.values():
                    # this is a single-read fasta file
                    read_group = read['Raw']
                    read_id = get_read_id(read_group)
                    index[read_id] = (fast5, True)
    return index


class Fast5Writer:
    def __init__(self, dir, prefix, fastq=False, max_reads=4000):
        self.dir = dir
        self.prefix = prefix
        self.fastq = fastq
        self.max_reads = max_reads
        self.batch = 0
        self.batch_count = 0
        self.read_count = 0
        self.fast5file = None
        self.fastqfile = None


    def add_read(self, read):
        self.ensure_output_files()
        read.copy(read, self.fast5file)
        if self.fastq:
            self.write_fastq(read)
        self.read_count += 1
        self.batch_count += 1
        if self.batch_count >= self.max_reads:
            self.fast5file.close()
            self.fast5file = None
            self.batch += 1
            self.batch_count = 0


    def write_fastq(self, read):
        fastq = get_best_fastq_location(read)
        if fastq:
            #print(read[fastq].value, end='', file=self.fastqfile)
            self.fastqfile.write(read[fastq][()])


    def ensure_output_files(self):
        if not self.fast5file:
            filename = os.path.join(self.dir, self.prefix + '_' + str(self.batch) + '.fast5')
            self.fast5file = h5py.File(filename, 'w')
        if self.fastq and not self.fastqfile:
            # a hack to put fastq at same level as per-barcode fast5 dir
            self.fastqfile = gzip.open(self.dir + '.fastq.gz', 'wb')


    def close(self):
        if self.fast5file:
            self.fast5file.close()
        if self.fastqfile:
            self.fastqfile.close()


def get_best_fastq_location(read):
    """
    This function returns the path in the FAST5 file to the best FASTQ. If there are multiple
    basecall locations, it returns the last one (hopefully from the most recent basecalling).

    Taken from rrwick's fast5_to_fastq
    """
    names = []
    read.visit(names.append)

    basecall_locations = sorted([x for x in names if x.upper().endswith('FASTQ')])
    two_d_locations = [x for x in basecall_locations if 'BASECALLED_2D' in x.upper()]
    template_locations = [x for x in basecall_locations if 'TEMPLATE' in x.upper()]
    complement_locations = [x for x in basecall_locations if 'COMPLEMENT' in x.upper()]

    # If the read has 2D basecalling, then that's what we use.
    if two_d_locations:
        return two_d_locations[-1]

    # If the read has both template and complement basecalling, then we choose the best based on
    # mean qscore.
    elif template_locations and complement_locations:
        template_location = template_locations[-1]
        complement_location = complement_locations[-1]
        mean_template_qscore = get_mean_score(hdf5_file, template_location)
        mean_complement_qscore = get_mean_score(hdf5_file, complement_location)
        if mean_template_qscore >= mean_complement_qscore:
            return template_location
        else:
            return complement_location

    # If the read has only template basecalling (normal for 1D) or only complement, then that's
    # what we use.
    elif template_locations:
        return template_locations[-1]
    elif complement_locations:
        return complement_locations[-1]

    # If the read has none of the above, but still has a fastq value in its hdf5, that's weird, but
    # we'll consider it a 1d read and use it.
    elif basecall_locations:
        return basecall_locations[-1]

    return None
