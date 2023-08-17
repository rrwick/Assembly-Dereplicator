"""
This module contains some tests for Assembly Dereplicator. To run them, execute `python3 -m pytest`
from the root Assembly Dereplicator directory.

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

import collections
import dereplicator
import glob
import gzip
import os
import pathlib
import pytest
import tempfile


def test_get_default_thread_count():
    thread_count = dereplicator.get_default_thread_count()
    assert thread_count >= 1
    assert thread_count <= 16


def test_get_colours_from_tput():
    colours = dereplicator.get_colours_from_tput()
    assert colours >= 1


def test_compression_type_1():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'uncompressed'
    compression_type = dereplicator.get_compression_type(filename)
    assert compression_type == 'plain'
    open_func = dereplicator.get_open_func(filename)
    assert open_func == open


def test_compression_type_2():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'gzipped'
    compression_type = dereplicator.get_compression_type(filename)
    assert compression_type == 'gz'
    open_func = dereplicator.get_open_func(filename)
    assert open_func == gzip.open


def test_compression_type_3():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'bzip2ed'
    with pytest.raises(SystemExit) as exit_message:
        dereplicator.get_compression_type(filename)
    assert 'cannot use bzip2' in str(exit_message.value)


def test_compression_type_4():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'zipped'
    with pytest.raises(SystemExit) as exit_message:
        dereplicator.get_compression_type(filename)
    assert 'cannot use zip' in str(exit_message.value)


def test_get_assembly_n50_1():
    filename = pathlib.Path(__file__).resolve().parent / 'assemblies' / '1.fasta'
    n50 = dereplicator.get_assembly_n50(filename)
    assert n50 == 100000


def test_get_assembly_n50_2():
    filename = pathlib.Path(__file__).resolve().parent / 'assemblies' / '2.fasta.gz'
    n50 = dereplicator.get_assembly_n50(filename)
    assert n50 == 68315


def test_get_assembly_n50_3():
    filename = pathlib.Path(__file__).resolve().parent / 'assemblies' / '4.fna.gz'
    n50 = dereplicator.get_assembly_n50(filename)
    assert n50 == 52598


def test_get_assembly_n50_4():
    filename = pathlib.Path(__file__).resolve().parent / 'other' / 'empty_file'
    n50 = dereplicator.get_assembly_n50(filename)
    assert n50 == 0


def test_get_assembly_n50_5():
    filename = pathlib.Path(__file__).resolve().parent / 'other' / 'extra_line_breaks.fasta'
    n50 = dereplicator.get_assembly_n50(filename)
    assert n50 == 100


def test_find_all_assemblies_1():
    directory = pathlib.Path(__file__).resolve().parent / 'assemblies'
    assemblies = dereplicator.find_all_assemblies(directory)
    assemblies = [os.path.basename(a) for a in assemblies]
    assert assemblies == ['1.fasta', '2.fasta.gz', '3.fna', '4.fna.gz', '5.fa', '6.fa.gz']


def test_find_all_assemblies_2():
    directory = pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri'
    assemblies = dereplicator.find_all_assemblies(directory)
    assemblies = [os.path.basename(a) for a in assemblies]
    assert assemblies == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                          'GCF_003214135.1.fna.gz', 'GCF_003214255.1.fna.gz',
                          'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                          'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_distance_dereplication_1():
    """
    With a distance of 1%, only 2.fasta.gz should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.01', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '4.fna.gz', '5.fa', '6.fa.gz']


def test_distance_dereplication_2():
    """
    With a distance of 3.5%, 2.fasta.gz and 4.fna.gz should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.035', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '5.fa', '6.fa.gz']


def test_distance_dereplication_3():
    """
    With a distance of 10%, 2.fasta.gz, 4.fna.gz and 5.fa should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.1', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '6.fa.gz']


def test_distance_dereplication_4():
    """
    When given a very low distance, all assemblies will be returned.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.0001', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214135.1.fna.gz', 'GCF_003214255.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_distance_dereplication_5():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.001', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214135.1.fna.gz', 'GCF_003214655.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz',
                                   'GCF_003215515.1.fna.gz']


def test_distance_dereplication_6():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.002', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_distance_dereplication_7():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.003', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_distance_dereplication_8():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.004', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz']


def test_distance_dereplication_9():
    """
    When given a very high distance, the highest N50 assembly will be returned.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.1', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003215265.1.fna.gz']


def test_count_dereplication_1():
    """
    When given a count of one, the highest N50 assembly will be returned.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '1', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003215265.1.fna.gz']


def test_count_dereplication_2():
    """
    When given a count equal to the size of the input set, all assemblies will be returned.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '8', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214135.1.fna.gz', 'GCF_003214255.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_count_dereplication_3():
    """
    When given a count greater than the size of the input set, all assemblies will be returned.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '100', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214135.1.fna.gz', 'GCF_003214255.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_count_dereplication_4():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '7', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214135.1.fna.gz', 'GCF_003214655.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz',
                                   'GCF_003215515.1.fna.gz']


def test_count_dereplication_5():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '6', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_count_dereplication_6():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '5', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz',
                                   'GCF_003215515.1.fna.gz']


def test_count_dereplication_7():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '4', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz']


def test_count_dereplication_8():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '3', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz']


def test_count_dereplication_9():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '2', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz']


def test_distance_and_count_dereplication_1():
    """
    When both --distance and --count are used both conditions must be satisfied.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.0001', '--count', '2', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz']


def test_distance_and_count_dereplication_2():
    """
    When both --distance and --count are used both conditions must be satisfied.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.004', '--count', '6', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz']


def test_fraction_dereplication_1():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--fraction', '0.5', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz']


def test_fraction_dereplication_2():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--fraction', '0.75', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003214655.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz', 'GCF_003215515.1.fna.gz']


def test_fraction_dereplication_3():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--fraction', '0.25', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz']


def test_fraction_dereplication_4():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--fraction', '0.000001', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003215265.1.fna.gz']


def test_count_and_fraction_dereplication_1():
    """
    When both --count and --fraction are used, the end condition is the lower number.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '2', '--fraction', '0.9', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz']


def test_count_and_fraction_dereplication_2():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--count', '7', '--fraction', '0.5', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003213895.1.fna.gz',
                                   'GCF_003215265.1.fna.gz', 'GCF_003215395.1.fna.gz']


def test_distance_and_count_and_fraction_dereplication_1():
    """
    When both --distance, --count and --fraction are all used, both distance and count/fraction
    (whichever is lower) must be satisfied.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.000001', '--count', '2', '--fraction', '0.9',
                           in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz']


def test_distance_and_count_and_fraction_dereplication_2():
    """
    When both --distance, --count and --fraction are all used, both distance and count/fraction
    (whichever is lower) must be satisfied.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.000001', '--count', '10', '--fraction', '0.125',
                           in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003215265.1.fna.gz']


def test_distance_and_count_and_fraction_dereplication_3():
    """
    When both --distance, --count and --fraction are all used, both distance and count/fraction
    (whichever is lower) must be satisfied.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.004', '--count', '10', '--fraction', '0.9',
                           in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz']


def test_verbose_dereplication():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'sulcia_muelleri')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--distance', '0.004', '--verbose', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['GCF_003213775.1.fna.gz', 'GCF_003215265.1.fna.gz',
                                   'GCF_003215395.1.fna.gz']


def test_no_assemblies():
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'compression')
    with tempfile.TemporaryDirectory() as out_dir:
        with pytest.raises(SystemExit) as sysexit:
            dereplicator.main(['--distance', '0.01', in_dir, out_dir])


def test_help_1():
    with pytest.raises(SystemExit) as sysexit:
        dereplicator.main(['--help'])
        assert sysexit.code == 0


def test_help_2():
    with pytest.raises(SystemExit) as sysexit:
        dereplicator.main([])
        assert sysexit.code == 1


def test_check_args():
    Args = collections.namedtuple('Args', ['in_dir', 'out_dir', 'distance', 'count', 'fraction'])
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=None, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=0.0, count=None, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=-0.1, count=None, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=2.0, count=None, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=0, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=-1, fraction=None))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=None, fraction=-0.1))
    with pytest.raises(SystemExit):
        dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=None, fraction=2.0))
    dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=0.001, count=None, fraction=None))
    dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=100, fraction=None))
    dereplicator.check_args(Args(in_dir='in', out_dir='out', distance=None, count=None, fraction=0.5))
