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


def test_entire_script_1():
    """
    With a thresholds of 1%, only 2.fasta.gz should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--threshold', '0.01', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '4.fna.gz', '5.fa', '6.fa.gz']


def test_entire_script_2():
    """
    With a thresholds of 3.5%, 2.fasta.gz and 4.fna.gz should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--threshold', '0.035', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '5.fa', '6.fa.gz']


def test_entire_script_3():
    """
    With a thresholds of 10%, 2.fasta.gz, 4.fna.gz and 5.fa should be excluded.
    """
    in_dir = str(pathlib.Path(__file__).resolve().parent / 'assemblies')
    with tempfile.TemporaryDirectory() as out_dir:
        dereplicator.main(['--threshold', '0.1', in_dir, out_dir])
        derep_assembles = sorted(glob.glob(out_dir + '/*'))
        derep_assembles = [os.path.basename(a) for a in derep_assembles]
        assert derep_assembles == ['1.fasta', '3.fna', '6.fa.gz']


def test_help_1():
    with pytest.raises(SystemExit) as sysexit:
        dereplicator.main(['--help'])
        assert sysexit.code == 0


def test_help_2():
    with pytest.raises(SystemExit) as sysexit:
        dereplicator.main([])
        assert sysexit.code == 1
