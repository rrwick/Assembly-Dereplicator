#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Assembly-dereplicator

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import multiprocessing
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

__version__ = '0.1.0'


def get_arguments(args):
    parser = MyParser(description='Assembly dereplicator', add_help=False,
                      formatter_class=MyHelpFormatter)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('in_dir', type=str,
                               help='Directory containing all assemblies')
    required_args.add_argument('out_dir', type=str,
                               help='Directory where dereplicated assemblies will be copied')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--threshold', type=float, default=0.005,
                              help='Mash distance clustering threshold')
    setting_args.add_argument('--sketch_size', type=int, default=10000,
                              help='Mash assembly sketch size')
    setting_args.add_argument('--batch_size', type=int, default=1000000,
                              help='Dereplication iterations will occur on random batches of this '
                                   'many assemblies - smaller numbers will reduce memory usage')
    setting_args.add_argument('--threads', type=int, default=get_default_thread_count(),
                              help='Number of CPU threads for Mash')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Display more clustering information')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='Dereplicate assemblies v' + __version__,
                            help="Show program's version number and exit")

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    random.seed(0)
    all_assemblies = find_all_assemblies(args.in_dir)
    initial_count = len(all_assemblies)
    os.makedirs(args.out_dir, exist_ok=True)

    excluded_assemblies = set()
    while True:
        if len(all_assemblies) <= args.batch_size:
            break

        print('Running dereplication on a random batch of {} assemblies...'.format(args.batch_size))
        random.shuffle(all_assemblies)
        batch_assemblies = all_assemblies[:args.batch_size]
        newly_excluded = dereplicate(batch_assemblies, args.threads, args.threshold,
                                     args.sketch_size, args.verbose)

        if len(newly_excluded) == 0:
            print('  no clusters found\n')
            break
        else:
            excluded_assemblies |= newly_excluded

        all_assemblies = [x for x in all_assemblies if x not in excluded_assemblies]
        print('{:,} assemblies remain\n'.format(len(all_assemblies)))

    if initial_count <= args.batch_size:
        print('Running dereplication on all {} assemblies...'.format(len(all_assemblies)))
    else:
        print('Running a final dereplication on all {} assemblies...'.format(len(all_assemblies)))
    excluded_assemblies |= dereplicate(all_assemblies, args.threads, args.threshold,
                                       args.sketch_size, args.verbose)
    all_assemblies = [x for x in all_assemblies if x not in excluded_assemblies]

    print('\nFinal dereplication: {:,} / {:,} assemblies'.format(len(all_assemblies),
                                                                 initial_count))
    print('Copying dereplicated assemblies to {}'.format(args.out_dir))
    for a in all_assemblies:
        shutil.copy(a, args.out_dir)
    print()


def dereplicate(all_assemblies, threads, threshold, sketch_size, verbose):
    all_assemblies = sorted(all_assemblies)
    excluded_assemblies = set()

    with tempfile.TemporaryDirectory() as temp_dir:
        mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir, sketch_size)
        pairwise_distances = pairwise_mash_distances(mash_sketch, threads)
        assemblies, graph = create_graph_from_distances(pairwise_distances, threshold)
        clusters = cluster_assemblies(assemblies, graph)

        for assemblies in clusters:
            if len(assemblies) > 1:
                n50, representative = sorted([(get_assembly_n50(a), a) for a in assemblies])[-1]
                rep_name = os.path.basename(representative)
                if verbose:
                    print(os.path.basename(representative) + '*,', end='')
                    non_rep_assemblies = [os.path.basename(a) for a in assemblies
                                          if a != representative]
                    print(','.join(non_rep_assemblies))
                else:
                    print('  cluster of {} assemblies: {} (N50 = {:,})'.format(len(assemblies),
                                                                               rep_name, n50))
                excluded_assemblies |= set([x for x in assemblies if x != representative])
            elif verbose:
                assert len(assemblies) == 1
                print(os.path.basename(assemblies[0]) + '*')

    return excluded_assemblies


def build_mash_sketch(assemblies, threads, temp_dir, sketch_size):
    mash_command = ['mash', 'sketch', '-p', str(threads), '-o', temp_dir + '/mash',
                    '-s', str(sketch_size)] + assemblies
    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return temp_dir + '/mash.msh'


def pairwise_mash_distances(mash_sketch, threads):
    mash_command = ['mash', 'dist', '-p', str(threads), mash_sketch, mash_sketch]
    mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    return mash_out.splitlines()


def find_all_assemblies(in_dir):
    print('\nLooking for assembly files in {}:'.format(in_dir))
    all_assemblies = [str(x) for x in sorted(pathlib.Path(in_dir).glob('**/*'))
                      if x.is_file()]
    all_assemblies = [x for x in all_assemblies if
                      x.endswith('.fasta') or x.endswith('.fasta.gz') or
                      x.endswith('.fna') or x.endswith('.fna.gz') or
                      x.endswith('.fa') or x.endswith('.fa.gz')]
    print('found {:,} files\n'.format(len(all_assemblies)))
    return all_assemblies


def create_graph_from_distances(pairwise_distances, threshold):
    """
    Builds an undirected graph where nodes are assemblies and edges connect assemblies which have
    a pairwise Mash distance below the threshold.
    """
    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)
    for line in pairwise_distances:
        parts = line.split('\t')
        assembly_1 = parts[0]
        assembly_2 = parts[1]
        distance = float(parts[2])
        assemblies.add(assembly_1)
        assemblies.add(assembly_2)
        if assembly_1 == assembly_2:
            continue
        all_connections[assembly_1].add(assembly_2)
        all_connections[assembly_2].add(assembly_1)
        if distance < threshold:
            graph[assembly_1].add(assembly_2)
            graph[assembly_2].add(assembly_1)
    assemblies = sorted(assemblies)
    assembly_count = len(assemblies)
    for assembly in assemblies:  # sanity check: make sure we have all the connections
        assert len(all_connections[assembly]) == assembly_count - 1
    return assemblies, graph


def cluster_assemblies(assemblies, graph):
    visited = set()
    clusters = []
    for assembly in assemblies:
        if assembly in visited:
            continue
        connected = dfs(graph, assembly)
        clusters.append(sorted(connected))
        visited |= connected
    return clusters


def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


def get_assembly_n50(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            return contig_length
    return 0


def get_contig_lengths(filename):
    lengths = []
    with get_open_func(filename)(str(filename), 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if the script is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 1:  # if no arguments were given.
            self.print_help(file=sys.stderr)
            sys.exit(1)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default' not in help_text.lower():
                help_text += ' (default: {})'.format(action.default)
            elif 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT',
                                              'default: {}'.format(action.default))
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))
        elif not action_header.endswith('\n'):
            parts.append('\n')
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)


if __name__ == '__main__':
    main()
